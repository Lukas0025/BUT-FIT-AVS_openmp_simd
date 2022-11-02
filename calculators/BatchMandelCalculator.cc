/**
 * @file BatchMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data  = (int *)(aligned_alloc(64, height * width * sizeof(int)));
	zReal = (float *)(aligned_alloc(64, width * sizeof(float)));
	zImag = (float *)(aligned_alloc(64, width * sizeof(float)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	free(data);
	free(zImag);
	free(zReal);
	data = NULL;
}


int * BatchMandelCalculator::calculateMandelbrot () {

	int h2 = height / 2;

	#pragma omp parallel for
	for (int i = 0; i <= h2; i++) {
		
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
			data[rowIndex + j] = -1;
		}

		for (int block = 0; block < width / blockSize; block++) {
			
			int blockStart = blockSize * block;
			int blockEnd   = blockSize * block + blockSize;
			int sum        = 0;
			
			for (int iteration = 0; iteration < limit; ++iteration) {

				#pragma omp simd
				for (int j = blockStart; j < blockEnd; j++) {

					if (data[rowIndex + j] == -1) {
						float r2 = zReal[j] * zReal[j];
						float i2 = zImag[j] * zImag[j];

						if (r2 + i2 > 4.0f) {
							data[rowIndex + j] = iteration;
							sum += 1;
						} else {
							zImag[j] = 2.0f * zReal[j] * zImag[j] + imag;
							zReal[j] = r2 - i2 + (x_start + j * dx);
						}
					}
				}

				if (sum == blockSize) break;
			}
		}
	}

	//correct offset
	if (height % 2 == 0) {
		height -= 1;
	}

	//copy half
	for (int i = h2; i < height; i++) {

		int cpRowIndex = (height - i) * width;
		int rowIndex   = i * width;
		
		#pragma omp simd
		for (int j = 0; j < width; j++) {
			data[rowIndex + j] = data[cpRowIndex + j];
		}
	}

	return data;
}
