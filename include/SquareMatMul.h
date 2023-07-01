/*
 * SquareMatMul.h
 *
 *  Created on: Nov 6, 2022
 *      Author: cagnur
 */

#ifndef SQUAREMATMUL_H_
#define SQUAREMATMUL_H_
#include "MatrixGen.h"

class SquareMatMul {
public:
	int bsize=100;
	int col;
	int row;
	int point;
	float **C, **firstMat, **secondMat;
	SquareMatMul(MatrixGen first, MatrixGen second);
	virtual ~SquareMatMul();
	// Experiment-1
	void reset(void);
	void ijk(void);
	void jik(void);
	void jki(void);
	void kji(void);
	void kij(void);
	void ikj(void );
	void bijk(void);
	void possibleIndexingTimeAnalysis(void);
	// Experiment-2
	void unrollingTwo(void);
	void unrollingFour(void);
	void unrollingEight(void);
	void unrollingSixteen(void);
	void onlyUnrollingTimeAnalysis(void);
	// Experiment-3
	void unrollingFourWithFission(void);
	void unrollingEightWithFission(void);
	void unrollingSixteenWithFission(void);
	void unrollingWithFissionTimeAnalysis(void);
	// Experiment-4
	void changeBlockSize(int blockSize);
	void blockSizeTimeAnalysis(int start, int stop, int stepSize);
	// Experiment-5
	void printMat(void);
	void bikj(void);
	void bijk_updated(void);
	void loopBlockingTimeAnalysis(void);
	// Experiment-6
	void bikjUnrolledFour(void);
	void bikjUnrolledEight(void);
	void bikjUnrolledSixteen(void);
	void bikjUnrolledTimeAnalysis(void);
};

#endif /* SQUAREMATMUL_H_ */
