/**
 * @file BatchMandelCalculator.h
 * @author Lukáš Plevač <xpleva07@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */
#ifndef BATCHMANDELCALCULATOR_H
#define BATCHMANDELCALCULATOR_H

#define blockSize 256

// Points of main cluter in mandel
#define mainCluterIStart 0.33
#define mainCluterIEnd   0.67
#define mainCluterRStart 0.50
#define mainCluterREnd   0.73

#include <BaseMandelCalculator.h>

class BatchMandelCalculator : public BaseMandelCalculator
{
public:
    BatchMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~BatchMandelCalculator();
    int * calculateMandelbrot();

private:
    int *data;
    float *zReal;
    float *zImag;
};

#endif
