/**
 * @file LineMandelCalculator.h
 * @author Lukáš Plevač <xpleva07@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */

#include <BaseMandelCalculator.h>

class LineMandelCalculator : public BaseMandelCalculator
{
public:
    LineMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~LineMandelCalculator();
    int *calculateMandelbrot();

private:
    int *data;    /* 64 aligned */
    float *zReal; /* 64 aligned */
    float *zImag; /* 64 aligned */
};
