/**
 * @file LineMandelCalculator.cc
 * @author Lukáš Plevač <xpleva07@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 1.11.2022
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	data = (int *)(malloc(height * width * sizeof(int)));
}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	data = NULL;
}

int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;

	float zReal[width];
	float zImag[width];

	for (int i = 0; i < height; i++) {
		
		int rowIndex = i * width;
		float imag = y_start + i * dy;

		//init values
		#pragma omp simd
		for (int j = 0; j < width; j++) {
			zReal[j] = x_start + j * dx; // current real value
			zImag[j] = imag;             // current imaginary value
		}

		//fill empty
		#pragma omp simd
		for (int j = 0; j < width; j++) {
			pdata[rowIndex + j] = -1;
		}

		for (int iteration = 0; iteration < limit; ++iteration) {

			#pragma omp simd
			for (int j = 0; j < width; j++) {
				if (pdata[rowIndex + j] == -1) {
					float r2 = zReal[j] * zReal[j];
					float i2 = zImag[j] * zImag[j];

					if (r2 + i2 > 4.0f) {
						pdata[rowIndex + j] = iteration;
					} else {
						zImag[j] = 2.0f * zReal[j] * zImag[j] + imag;
						zReal[j] = r2 - i2 + (x_start + j * dx);
					}
				}
			}
		}
		
	}

	return data;
}
