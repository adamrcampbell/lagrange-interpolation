

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct DoubleComplex {
        double real;
        double imaginary;
    } DoubleComplex;

    double calcDistance(double width);
    double calcSphrShift(double index, double width);
    double calcInterpShift(double index, double width);
    double calcSphrSample(double x);
    double calcResolutionShift(double index, double width);
    int calcRelativeIndex(float x, int width);
    void interpolateKernel(DoubleComplex *source, DoubleComplex *destination, 
            int resolutionSupport, int textureSupport);
    void getBicubicNeighbours(int rowIndex, int colIndex, DoubleComplex *n, double *s,
            int sourceSupport, DoubleComplex *source);
    DoubleComplex interpolateCubicSample(DoubleComplex z0, DoubleComplex z1, 
            DoubleComplex z2, DoubleComplex z3, double x0, double x1, double x2,
            double x3, double h, double x);
    void printMatrix(DoubleComplex *matrix, int support);
    
    DoubleComplex complexScale(DoubleComplex z, double scalar);
    DoubleComplex complexAdd(DoubleComplex z0, DoubleComplex z1);
    

#ifdef __cplusplus
}
#endif

#endif /* INTERPOLATION_H */

