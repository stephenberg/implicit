#pragma once

#ifndef __GRID__
#define __GRID__


#include <RcppEigen.h>

//container class for quantities relating to a 2d spatial grid
class Grid {

public:
	VectorXi internalPoints; //indicator for whether each point is internal or on the boundary
	int rows, cols, rows_internal, cols_internal; //number of internal rows, internal columns
	int nInternal, nFull;  //number of locations in internal, full grid 
	double lengthX, lengthY; //horizontal and vertical lengths of the grid, in units (m, km, miles, ...)
	bool regularGrid;

	Grid(int rows_, int cols_, double lengthX_ = 1, double lengthY_ = 1) {

		//set grid dimensions
		rows = rows_;
		cols = cols_;
		rows_internal = rows - 2;
		cols_internal = cols - 2;
		lengthX = lengthX_;
		lengthY = lengthY_;
		nFull = rows * cols;
		nInternal = rows_internal * cols_internal;

		//construct internal points: all points internal, except for
		//top and bottom rows, left and right columns
		internalPoints.setConstant(nFull, 1);

		Map<MatrixXi> internalMatrix(internalPoints.data(), rows, cols);
		internalMatrix.col(0).setZero();
		internalMatrix.col(cols - 1).setZero();
		internalMatrix.row(0).setZero();
		internalMatrix.row(rows - 1).setZero();
		regularGrid = true;//no internalPoints vector being passed in

	}

	Grid(int rows_, int cols_, VectorXi internalPoints_, double lengthX_ = 1, double lengthY_ = 1) {
		regularGrid = false; //assume the internalPoints vector being passed in is for a non-rectangular grid

		//set grid dimensions
		rows = rows_;
		cols = cols_;
		rows_internal = rows - 2;
		cols_internal = cols - 2;
		lengthX = lengthX_;
		lengthY = lengthY_;
		nFull = rows * cols;
		nInternal = rows_internal * cols_internal;

		internalPoints = internalPoints_;
	}
};

#endif