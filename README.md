# lagrange-interpolation
Performs largrange polynomial bicubic interpolation on two-dimensional matrices. This project is specifically tailored for the  interpolation of custom W-Projection convolutional kernels for a currently private research project.

The requirements for this interpolation process are to scale from some two-dimensional _source_ matrix, down to some two-dimensional _destination_ matrix. It assumes the dimensions of both matrices are some power of two, and that during the interpolation process, the entire first row and colum of the _source_ matrix are not included during interpolation. This is needed as the first row and column contain redundant information obtained during the FFT of the _source_ matrix during the creation of kernels.
