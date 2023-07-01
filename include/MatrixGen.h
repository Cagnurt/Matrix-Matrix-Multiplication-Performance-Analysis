/*
 * MatrixGen.h
 *
 *  Created on: Oct 28, 2022
 *      Author: cagnur
 */

#ifndef MATRIXGEN_H_
#define MATRIXGEN_H_

class MatrixGen {
private:
	float lb = 0;
	float ub = 1;

public:
	float **data;
	int col;
	int row;
	virtual ~MatrixGen();
	//MatrixGen(void);
	MatrixGen(int dim);
	float randBtwFloat(void);
	void fill(void);
	void print(void);
	//MatrixGen& operator=(const MatrixGen& obj);

};

#endif /* MATRIXGEN_H_ */
