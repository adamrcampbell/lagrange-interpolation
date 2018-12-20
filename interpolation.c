
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interpolation.h"

int main(int argc, char** argv) {

    // Note: only functional on even sized resolution and texture sizes, for now
    
    int resSize = 64;
    int texSize = 32;
    
    DoubleComplex *resolution = calloc(resSize * resSize, sizeof(DoubleComplex));
    DoubleComplex *texture = calloc(texSize * texSize, sizeof(DoubleComplex));
    
    double width = (double) resSize;
    double rowTaper, rowShift, colShift, taper = 0.0;
    
    for(int r = 0; r < resSize; r++)
    {   
        // calculate current shift in range [-1.0, 1.0]
        rowShift = fabs(calcSphrShift(r, width));
        // calculate 1D spheroidal sample
        rowTaper = calcSphrSample(rowShift) * (1.0 - rowShift * rowShift);
        
        for(int c = 0; c < resSize; c++)
        {
            // calculate current shift in range [-1.0, 1.0]
            colShift = fabs(calcSphrShift(c, width));
            // calculate 2D spheroidal sample
            taper = rowTaper * calcSphrSample(colShift) * (1.0 - colShift * colShift);
            resolution[r * resSize + c] = (DoubleComplex) {.real = taper, .imaginary = 0.0};
        }
    }
    
    printf("Printing resolution:\n");
    printMatrix(resolution, resSize);
    printf("\n");
    
    interpolateKernel(resolution, texture, resSize, texSize);
    
    printf("Printing texture:\n");
    printMatrix(texture, texSize);
    printf("\n");
    
    return (EXIT_SUCCESS);
}

void interpolateKernel(DoubleComplex *source, DoubleComplex *destination, 
    int sourceSupport, int destinationSupport)
{
    // Determine distance between source samples in range [-1.0, 1.0]
    double sourceShift = calcDistance(sourceSupport-1);
    // Storage for neighbours, synthesized points
    DoubleComplex n[16], p[4];
    // Neighbours shift values (rs = row shift, cs = col shift)
    double rs[16], cs[16];
    double rowShift, colShift = 0.0;
    
    for(int r = 0; r < destinationSupport; r++)
    {
        // Determine relative shift for interpolation row [-1.0, 1.0]
        rowShift = calcSphrShift((double) r, (double) destinationSupport-1);
        
        for(int c = 0; c < destinationSupport; c++)
        {
            // Determine relative shift for interpolation col [-1.0, 1.0]
            colShift = calcSphrShift((double) c, (double) destinationSupport-1);
            
            // gather 16 neighbours
            getBicubicNeighbours(rowShift, colShift, n, rs, cs, sourceSupport, source);
            // interpolate intermediate sample
            p[0] = interpolateCubicSample(n[0], n[1], n[2], n[3],
                cs[0], cs[1], cs[2], cs[3], sourceShift, colShift);
            p[1] = interpolateCubicSample(n[4], n[5], n[6], n[7],
                cs[4], cs[5], cs[6], cs[7], sourceShift, colShift);
            p[2] = interpolateCubicSample(n[8], n[9], n[10], n[11],
                cs[8], cs[9], cs[10], cs[11], sourceShift, colShift);
            p[3] = interpolateCubicSample(n[12], n[13], n[14], n[15],
                cs[12], cs[13], cs[14], cs[15], sourceShift, colShift);
            
            // interpolate final sample
            destination[r * destinationSupport + c] = interpolateCubicSample(p[0], p[1], p[2], p[3],
               rs[1], rs[5], rs[9], rs[13], sourceShift, rowShift);
        }
    }
}

void getBicubicNeighbours(double rowShift, double colShift, DoubleComplex *n, double *rs, double *cs,
        int sourceSupport, DoubleComplex *source)
{
    // determine where to start locating neighbours in source matrix
    int x = calcRelativeIndex(colShift, (double) sourceSupport);
    int y = calcRelativeIndex(rowShift, (double) sourceSupport);
    // counter for active neighbour
    int nIndex = 0;
    // define neighbour boundaries
    int rStart = (rowShift < 0.0) ? y-1 : y-2;
    int rEnd = (rowShift < 0.0) ? y+3 : y+2;
    int cStart = (colShift < 0.0) ? x-1 : x-2;
    int cEnd = (colShift < 0.0) ? x+3 : x+2;
    
    // gather 16 neighbours
    for(int r = rStart; r < rEnd; r++)
    {   
        for(int c = cStart; c < cEnd; c++)
        {
            // set row and col shifts for neighbour
            rs[nIndex] = (rowShift < 0.0) ? calcSphrShift(r-1, sourceSupport-1) : calcSphrShift(r, sourceSupport-1);
            cs[nIndex] = (colShift < 0.0) ? calcSphrShift(c-1, sourceSupport-1) : calcSphrShift(c, sourceSupport-1);            
            // neighbour falls out of bounds
            if(r < 1 || c < 1 || r >= sourceSupport || c >= sourceSupport)
                n[nIndex] = (DoubleComplex) {.real = 0.0, .imaginary = 0.0};
            // neighbour exists
            else
                n[nIndex] = source[r * sourceSupport + c];

            nIndex++;
        }
    }   
}

DoubleComplex interpolateCubicSample(DoubleComplex z0, DoubleComplex z1, 
    DoubleComplex z2, DoubleComplex z3, double x0, double x1, double x2,
    double x3, double h, double x)
{
    double hCube = pow(h, 3.0);
    double scale0 = -(x-x1)*(x-x2)*(x-x3)/(6.0*hCube);
    double scale1 = (x-x0)*(x-x2)*(x-x3)/(2.0*hCube);
    double scale2 = -(x-x0)*(x-x1)*(x-x3)/(2.0*hCube);
    double scale3 = (x-x0)*(x-x1)*(x-x2)/(6.0*hCube);
   
    DoubleComplex z = complexScale(z0, scale0);
    z = complexAdd(z, complexScale(z1, scale1));
    z = complexAdd(z, complexScale(z2, scale2));
    z = complexAdd(z, complexScale(z3, scale3));
    
    return z;
}

double calcDistance(double width)
{
    return 2.0/width;
}

double calcSphrShift(double index, double width)
{   
    return -1.0 + index * calcDistance(width);
}

double calcInterpShift(double index, double width)
{
    return -1.0 + ((2.0 * index) / width);
}

//double calcInterpShift(double index, double width)
//{
//    return -1.0 + ((2.0 * index + 1.0) / width);
//}

double calcResolutionShift(double index, double width)
{
    return -1.0 + ((index) * (2.0 / (width-2.0)));
}

double calcSphrSample(double x)
{   
    static double p[] = {0.08203343, -0.3644705, 0.627866, -0.5335581, 0.2312756,
        0.004028559, -0.03697768, 0.1021332, -0.1201436, 0.06412774};
    static double q[] = {1.0, 0.8212018, 0.2078043,
        1.0, 0.9599102, 0.2918724};
    
    int part, sp, sq;
    double xend, delta, top, bottom;
    
    if(x >= 0.0 && x < 0.75)
    {
        part = 0;   
        xend = 0.75;
    }
    else if(x >= 0.75 && x <= 1.0)
    {
        part = 1;
        xend = 1.0;
    }
    else
        return 0.0;

    delta = x * x - xend * xend;
    sp = part * 5;
    sq = part * 3;
    top = p[sp];
    bottom = q[sq];
    
    for(int i = 1; i < 5; i++)
        top += p[sp+i] * pow(delta, i);
    for(int i = 1; i < 3; i++)
        bottom += q[sq+i] * pow(delta, i);
    return (bottom == 0.0) ? 0.0 : top/bottom;
}

DoubleComplex complexScale(DoubleComplex z, double scalar)
{
    return (DoubleComplex) {.real=z.real*scalar, .imaginary=z.imaginary*scalar};
}

DoubleComplex complexAdd(DoubleComplex x, DoubleComplex y)
{
    return (DoubleComplex) {.real = x.real + y.real, .imaginary=x.imaginary + y.imaginary};
}

int calcRelativeIndex(double x, double width)
{
    int offset = (x < 0.0) ? 1 : 2;
    return ((int) floor(((x+1.0f)/2.0f) * (width-offset)))+1;
}

void printMatrix(DoubleComplex *matrix, int support)
{
    for(int r = 0; r < support; r++)
    {
        for(int c = 0; c < support; c++)
            printf("%.3f, ", matrix[r * support + c].real);
        printf("\n");
    }
    printf("\n");
}
