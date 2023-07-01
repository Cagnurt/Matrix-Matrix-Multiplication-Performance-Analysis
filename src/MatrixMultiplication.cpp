//============================================================================
// Name        : MatrixMultiplication.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "SquareMatMul.h"
#include "MatrixGen.h"
using namespace std;

int main() {
	int dim = 2;
	MatrixGen A(dim);
	A.fill();
	cout<<"First matrix is generated!"<<endl;
	A.print();
	MatrixGen B(dim);
	B.fill();
	B.print();
	cout<<"Second matrix is generated!"<<endl;
	SquareMatMul mult(A, B);
/*
	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim; j++)
		{
			cout<< mult.firstMat[i][j]<<" ";
		}
		cout <<" "<<endl;
	}

	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim; j++)
		{
			cout<< mult.secondMat[i][j]<<" ";
		}
		cout << " "<<endl;
	}
*/
#ifdef Experiment-1
	int increment = 1000;
	for (int dim=8000; dim > 0; dim-=increment){
		cout<<"Dimension="<<dim<<":"<< endl;
		MatrixGen A(dim);
		A.fill();
		cout<<"First matrix is generated!"<<endl;
		MatrixGen B(dim);
		B.fill();
		cout<<"Second matrix is generated!"<<endl;
		SquareMatMul mult(A, B);
		mult.possibleIndexingTimeAnalysis(A.data, B.data);
	}
#endif

#ifdef Experiment-2
	int increment = 1000;
	for (int dim=1000; dim < 6000; dim+=increment){
		cout<<"Dimension="<<dim<<":"<< endl;
		MatrixGen A(dim);
		A.fill();
		cout<<"First matrix is generated!"<<endl;
		MatrixGen B(dim);
		B.fill();
		cout<<"Second matrix is generated!"<<endl;
		SquareMatMul mult(A, B);
		mult.onlyUnrollingTimeAnalysis(A.data, B.data);
	}
#endif

#ifdef Experiment-3
	int increment = 1000;
	for (int dim=1000; dim < 6000; dim+=increment){
		cout<<"Dimension="<<dim<<":"<< endl;
		MatrixGen A(dim);
		A.fill();
		cout<<"First matrix is generated!"<<endl;
		MatrixGen B(dim);
		B.fill();
		cout<<"Second matrix is generated!"<<endl;
		SquareMatMul mult(A, B);
		mult.unrollingWithFissionTimeAnalysis(A.data, B.data);
	}
#endif

#ifdef Experiment-4
	int start = 4;
	int stop = 257;
	int stepSize = 2;
	int dim = 2000;
	cout<<"Dimension="<<dim<<":"<< endl;
	MatrixGen A(dim);
	A.fill();
	cout<<"First matrix is generated!"<<endl;
	MatrixGen B(dim);
	B.fill();
	cout<<"Second matrix is generated!"<<endl;
	MatMul mult(dim);
	mult.blockSizeTimeAnalysis(A.data,B.data,start, stop,stepSize);
	cout<< " "<< endl;
	dim = 2048;
	cout<<"Dimension="<<dim<<":"<< endl;
	MatrixGen A(dim);
	A.fill();
	cout<<"First matrix is generated!"<<endl;
	MatrixGen B(dim);
	B.fill();
	cout<<"Second matrix is generated!"<<endl;
	SquareMatMul mult(A, B);
	mult.blockSizeTimeAnalysis(A.data,B.data,start, stop,stepSize);
#endif

#ifdef Experiment-5
	int dim = 4;
	cout<<"Dimension="<<dim<<":"<< endl;
	MatrixGen A(dim);
	A.fill();
	cout<<"First matrix is generated!"<<endl;
	MatrixGen B(dim);
	B.fill();
	cout<<"Second matrix is generated!"<<endl;
	SquareMatMul mult(A, B);
	mult.changeBlockSize(2);
	mult.bijk();
	mult.printMat();
	mult.reset();
	cout<<" "<< endl;
	mult.bijk_updated();
	mult.printMat();
	mult.reset();
	cout<<" "<< endl;
	mult.bikj();
	mult.printMat();
	mult.reset();
	cout<<"Time analysis part:"<< endl;
	int increment = 1000;
	for (dim=1000; dim < 5000; dim+=increment){
		cout<<"Dimension="<<dim<<":"<< endl;
		MatrixGen A(dim);
		A.fill();
		cout<<"First matrix is generated!"<<endl;
		MatrixGen B(dim);
		B.fill();
		cout<<"Second matrix is generated!"<<endl;
		SquareMatMul mult(A, B);
		mult.changeBlockSize(64);
		mult.loopBlockingTimeAnalysis();
	}
#endif

#ifdef Experiment-6
	int increment = 1000;
	for (int dim=4000; dim < 6000; dim+=increment){
		cout<<"Dimension="<<dim<<":"<< endl;
		MatrixGen A(dim);
		A.fill();
		cout<<"First matrix is generated!"<<endl;
		MatrixGen B(dim);
		B.fill();
		cout<<"Second matrix is generated!"<<endl;
		SquareMatMul mult(A, B);
		mult.changeBlockSize(64);
		mult.bikjUnrolledTimeAnalysis();
	}

#endif
	return 0;
}




