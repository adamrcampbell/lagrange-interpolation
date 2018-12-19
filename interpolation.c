
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interpolation.h"

int main(int argc, char** argv) {

    // Note: only functional on even sized resolution and texture sizes, for now
    
    int resSize = 16;
    int texSize = 16;
    
    DoubleComplex *resolution = calloc(resSize, sizeof(DoubleComplex));
    DoubleComplex *texture = calloc(texSize, sizeof(DoubleComplex));
    
    double width = (double) resSize;
    double taper, shift = 0.0;
    
    for(int i = 0; i < resSize; i++)
    {
        // printf("%f\n", calcSphrShift(i, width));
        shift = fabs(calcSphrShift(i, width)); // plus/minus 1.0 to offset for odd/even spheroidals
        taper = calcSphrSample(shift) * (1.0 - shift * shift);
        resolution[i] = (DoubleComplex) {.real = taper, .imaginary = 0.0};
    }
    
    printf("Original Curve: ");
    printMatrix(resolution, resSize);
    
    interpolateKernel(resolution, texture, resSize, texSize);
    
    printf("Interpol Curve: ");
    printMatrix(texture, texSize);
    
    
//    DoubleComplex n0 = (DoubleComplex) {.real = 0.0, .imaginary = 0.0};
//    DoubleComplex n1 = (DoubleComplex) {.real = 0.035890, .imaginary = 0.0};
//    DoubleComplex n2 = (DoubleComplex) {.real = 0.270799, .imaginary = 0.0};
//    DoubleComplex n3 = (DoubleComplex) {.real = 0.734366, .imaginary = 0.0};
//    double s0 = 1.285714;
//    double s1 = 1.0;
//    double s2 = 0.714286;
//    double s3 = 0.428571;
//    printf("%f\n", interpolateCubicSample(n3, n2, n1, n0, s3, s2, s1, s0, calcDistance(7), 1.0).real);
    
    return (EXIT_SUCCESS);
}

void interpolateKernel(DoubleComplex *source, DoubleComplex *destination, 
    int sourceSupport, int destinationSupport)
{
    double sourceShift = calcDistance(sourceSupport-1);
    
    DoubleComplex n[4];
    double s[4];
    int y[4];
    DoubleComplex p[4];
    
    for(int i = 0; i < destinationSupport; i++)
    {
        //double shift = calcInterpShift((double) i, (double) destinationSupport-1);
        double shift = calcSphrShift((double) i, (double) destinationSupport-1);
        // printf("%d: %f\n", i, shift);
        int j = calcRelativeIndex(shift, (double) sourceSupport);

        getBicubicNeighbours(shift, n, s, y, sourceSupport, source);
        
        printf("[%d:%f] [%d:%f] [%d:%f] [%d:%f] [%d:%f]\n", y[0], s[0], y[1], s[1], j, shift, y[2], s[2], y[3], s[3]);
        
        // printf("%d: %f %f %f %f\n", i, n[0].real, n[1].real, n[2].real, n[3].real);
        
        //printf("%d %d %d %d %d\n", s[0], s[1], j, s[2], s[3]);DoubleComplex n0 = (DoubleComplex) {.real = 0.0, .imaginary = 0.0};
        
        destination[i] = interpolateCubicSample(n[0], n[1], n[2], n[3],
                    s[0], s[1], s[2], s[3], sourceShift, shift);
    }
}

void getBicubicNeighbours(double shift, DoubleComplex *n, double *s, int *y,
        int sourceSupport, DoubleComplex *source)
{
    int j = calcRelativeIndex(shift, (double) sourceSupport);
    
    int nIndex = 0;
    
    // define neighbour boundaries
    int start = (shift < 0.0) ? j-1 : j-2;
    int end = (shift < 0.0) ? j+3 : j+2;
    
    for(int i = start; i < end; i++)
    {
        // Calculate each neighbours shift
        s[nIndex] = (shift < 0.0) ? calcSphrShift(i-1, sourceSupport-1) : calcSphrShift(i, sourceSupport-1) ;
        y[nIndex] = i;
        
        // Neighbour falls out of bounds
        if(i < 1 || i >= sourceSupport)
            n[nIndex] = (DoubleComplex) {.real = 0.0, .imaginary = 0.0};
        // Neighbour exists
        else
            n[nIndex] = source[i];

        nIndex++;
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
    for(int i = 0; i < support; i++)
    {
        printf("%.3f, ", matrix[i].real);
    }
    printf("\n");
}
