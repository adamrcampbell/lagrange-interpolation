
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interpolation.h"

int main(int argc, char** argv) {

    int resSize = 8;
    int texSize = 8;
    
    DoubleComplex *resolution = calloc(resSize * resSize, sizeof(DoubleComplex));
    DoubleComplex *texture = calloc(texSize * texSize, sizeof(DoubleComplex));
    
    double width = (double) resSize;
    double rowIndex, colIndex, rowShift, colShift, rowTaper, colTaper = 0.0;
    int matrixIndex = 0;
    
    for(int r = 0; r < resSize; r++)
    {
        rowIndex = (double) r;
        rowShift = fabs(calcSphrShift(rowIndex, width)); // plus/minus 1.0 to offset for odd/even spheroidals
        rowTaper = calcSphrSample(rowShift) * (1.0 - rowShift * rowShift);
        
        for(int c = 0; c < resSize; c++)
        {
            colIndex = (double) c;
            colShift = fabs(calcSphrShift(colIndex, width));
            colTaper = calcSphrSample(colShift) * (1.0 - colShift * colShift);
            matrixIndex = r * resSize + c;
            
            resolution[matrixIndex] = (DoubleComplex) {.real = rowTaper * colTaper, .imaginary = 0.0};
        }
    }
    
    printMatrix(resolution, resSize);
    
    interpolateKernel(resolution, texture, resSize, texSize);
    printf("\n\n");
    printMatrix(texture, texSize);
    
    return (EXIT_SUCCESS);
}

void interpolateKernel(DoubleComplex *source, DoubleComplex *destination, 
    int sourceSupport, int destinationSupport)
{
    double rowShift, colShift = 0.0;
    double shift = calcDistance(sourceSupport-2);
    int rowIndex, colIndex, destIndex = 0;
    
    DoubleComplex n[16];
    double s[16];
    DoubleComplex p[4];
    
    for(int r = 0; r < destinationSupport; r++)
    {
        rowShift = calcInterpShift((double) r, (double) destinationSupport);
        rowIndex = calcRelativeIndex(rowShift, sourceSupport-2)+1;
        
        for(int c = 0; c < destinationSupport; c++)
        {
            colShift = calcInterpShift((double) c, (double) destinationSupport);
            colIndex = calcRelativeIndex(colShift, sourceSupport-2)+1;
            
            getBicubicNeighbours(rowIndex, colIndex, n, s, sourceSupport, source);
            
            // Interpolate new samples
            for(int i = 0; i < 4; i++)
            {
                p[i] = interpolateCubicSample(n[i * 4], n[i * 4 + 1], n[i * 4 + 2], n[i * 4 + 3],
                        s[i * 4], s[i * 4 + 1], s[i * 4 + 2], s[i * 4 + 3], shift, colShift);
            }
            
            // Interpolate final sample
            destIndex = r * destinationSupport + c;
            destination[destIndex] = interpolateCubicSample(p[0], p[1], p[2], p[3], 
                    calcResolutionShift(rowIndex-1, sourceSupport),
                    calcResolutionShift(rowIndex, sourceSupport),
                    calcResolutionShift(rowIndex+1, sourceSupport),
                    calcResolutionShift(rowIndex+2, sourceSupport),
                    shift, rowShift);
        }        
    }
}

void getBicubicNeighbours(int rowIndex, int colIndex, DoubleComplex *n, double *s,
        int sourceSupport, DoubleComplex *source)
{
    int nIndex = 0;
    
    for(int r = rowIndex-1; r < rowIndex+3; r++)
    {
        for(int c = colIndex-1; c < colIndex+3; c++)
        {
            // Calculate each neighbours shift
            s[nIndex] = calcResolutionShift(c, sourceSupport);
            
            // printf("%f ", s[nIndex]);
            
            // Neighbour falls out of bounds
            if(r < 1 || r >= sourceSupport || c < 1 || c >= sourceSupport)
                n[nIndex] = (DoubleComplex) {.real = 0.0, .imaginary = 0.0};
            // Neighbour exists
            else
                n[nIndex] = source[r * sourceSupport + c];
            
            nIndex++;
        }
        
        // printf("\n");
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
    return -1.0 + ((2.0 * index + 1.0) / width);
}

double calcResolutionShift(double index, double width)
{
    return -1.0 + ((index-1.0) * (2.0 / (width-2.0)));
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

int calcRelativeIndex(float x, int width)
{
    int offset = (x < 0.0) ? 1 : 2;
    return ((int) floor(((x+1.0f)/2.0f) * (width-offset)))+1;
}

void printMatrix(DoubleComplex *matrix, int support)
{
    int index = 0;
    
    for(int r = 0; r < support; r++)
    {
        for(int c = 0; c < support; c++)
        {
            index = r * support + c;
            printf("%.3f,", matrix[index].real);
        }
        
        printf("\n");
    }
}
