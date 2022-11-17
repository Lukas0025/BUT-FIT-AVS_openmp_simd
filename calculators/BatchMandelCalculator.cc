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
#include <omp.h>
#include <cmath>
#include <string.h>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data  = (int *)(_mm_malloc(height * width * sizeof(int), 64));
	zReal = (float *)(_mm_malloc(width * sizeof(float),      64));
	zImag = (float *)(_mm_malloc(width * sizeof(float),      64));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	_mm_free(data);
	_mm_free(zImag);
	_mm_free(zReal);
	data = NULL;
}


int * BatchMandelCalculator::calculateMandelbrot () {

	const int h2 = (float) height / 2;

	const int mainGroupXBlockStart = std::ceil ((float) std::ceil (mainCluterRStart * width) / blockSize) + 1;
	const int mainGroupXBlockEnd   = std::floor((float) std::floor(mainCluterREnd   * width) / blockSize) - 1;
	const int mainGroupYStart      = std::ceil (mainCluterIStart * height);
	const int mainGroupYEnd        = std::floor(mainCluterIEnd   * height);

	//#pragma omp parallel for
	for (int i = 0; i <= h2; i++) {
		
		const int rowIndex = i * width;
		const float imag = y_start + i * dy;

		const int cpRowIndex = (height - i - 1) * width;

		//init values
		#pragma omp simd simdlen(64)
		for (int j = 0; j < width; j++) {
			zReal[j] = x_start + j * dx; // current real value
			zImag[j] = imag;             // current imaginary value
		}

		//fill empty
		#pragma omp simd simdlen(64)
		for (int j = 0; j < width; j++) {
			data[rowIndex + j] = limit;
		}

		for (int block = 0; block < std::ceil(((float) width) / blockSize); block++) {

			//drop it is in bandelbrot set
			if (block >= mainGroupXBlockStart && block <= mainGroupXBlockEnd && i >= mainGroupYStart && i <= mainGroupYEnd) {
				continue;
			}
			
			const int blockStart = blockSize * block;
			int blockEnd   = blockSize * block + blockSize;
			int sum        = blockSize;

			//last iteratin
			if (blockEnd > width) {
				blockEnd = width;
			}
		
			for (int iteration = 0; iteration < limit; ++iteration) {

				#pragma omp simd reduction(-: sum) simdlen(64)
				for (int j = blockStart; j < blockEnd; j++) {

					if (data[rowIndex + j] == limit) {
						float r2 = zReal[j] * zReal[j];
						float i2 = zImag[j] * zImag[j];

						if (r2 + i2 > 4.0f) {
							data[rowIndex + j] = iteration;
							sum = sum - 1;
						} else {
							zImag[j] = 2.0f * zReal[j] * zImag[j] + imag;
							zReal[j] = r2 - i2 + (x_start + j * dx);
						}
					}
				}

				if (!sum) break;
			}
		}

		//copy
		#pragma omp simd simdlen(64)
		for (int j = 0; j < width; j++) {
			data[cpRowIndex + j] = data[rowIndex + j];
		}
	}

	return data;
}
