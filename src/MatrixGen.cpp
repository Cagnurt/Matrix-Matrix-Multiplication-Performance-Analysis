/*
 * MatrixGen.cpp
 *
 *  Created on: Oct 28, 2022
 *      Author: cagnur
 */

#include "MatrixGen.h"
#include <cstdlib>
#include <cassert>
#include <iostream>
using namespace std;


MatrixGen::~MatrixGen() {
    for(int i = 0; i < row; ++i)
        delete [] data[i];

    delete [] data;
}
/*
MatrixGen::MatrixGen():col(3),row(3){
	data=new float *[row];
	for(int i=0;i<row;i++)
	{
	  data[i]=new float[col];
	}
}*/

MatrixGen::MatrixGen(int dim) {
	assert(dim != 0);
	col = dim;
	row = dim;
	data=new float *[row];
	for(int i=0;i<row;i++)
	{
	  data[i]=new float[col];
	}
}
/*
MatrixGen& MatrixGen::operator=(const MatrixGen& obj){
	this->col = obj.col;
	this->row = obj.row;
    for(int i = 0; i < row; ++i)
        delete [] data[i];

    delete [] data;
	this-> data = obj.data;
	return *this;
}*/


float MatrixGen::randBtwFloat(void)
{
	float range = (ub - lb);
	float div = RAND_MAX / range;
	return lb + (rand()/div);
}

void MatrixGen::fill(void)
{
	int i,j;
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			data[i][j] = this->randBtwFloat();
		}
	}
}

void MatrixGen::print(void)
{
	int i,j;
	cout << "Matrix is the following:"<< endl;
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			cout << data[i][j] << " ";
		}
		cout<<' '<<endl;
	}
}





