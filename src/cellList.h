#ifndef CELLLIST_H_
#define CELLLIST_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "atom.h"

/*!
 * cellList class. 
 */
class cellList {
private:
	void initNeighbours();
    
public:
	cellList (std::vector < double >, double, std::vector < atom* >);
    
	void sortIntoCells (std::vector < atom* >);
	void sortIntoCells (std::vector < atom >*);
	void insertParticle (atom*);
	void swapAndDeleteParticle (atom*, atom*);
	void translateParticle (atom*, std::vector < double >);
	
	int calcIndex (int, int, int);
	int calcIndexS (int, int, int);
	int calcIndex (double, double, double);
	int calcIndexS (double, double, double);
    
    std::vector < double > cellSize;
    std::vector < double > box;
    
    int cellCountX, cellCountY, cellCountZ, cellCountXY, cellCount;
    std::vector < std::vector < atom* > > cells;
    std::vector < std::vector < unsigned int > > neighbours;
};

inline int cellList::calcIndex (int k, int l, int m) {
	return (k + l*cellCountX + m*cellCountXY);
}

// this is the safe implementation which does not allow out of bound indices
inline int cellList::calcIndexS (int _k, int _l, int _m) {
	int k = _k;
	int l = _l;
	int m = _m;

	if (k >= cellCountX)
		k -= cellCountX;
	else if (k < 0)
		k += cellCountX;		
		
	if (l >= cellCountY)
		l -= cellCountY;
	else if (l < 0)
		l += cellCountY;		

	if (m >= cellCountZ)
		m -= cellCountZ;
	else if (m < 0)
		m += cellCountZ;
		
	return (k + l*cellCountX + m*cellCountXY);
}

inline int cellList::calcIndex (double _x, double _y, double _z) {
	return (floor(_x/cellSize[0]) + floor(_y/cellSize[1])*cellCountX +floor(_z/cellSize[2])*cellCountXY);
}

// this is the safe implementation which does not allow out of bound positions
inline int cellList::calcIndexS (double _x, double _y, double _z) {
	double x = _x;
	double y = _y;
	double z = _z;

	if (x >= box[0])
		x -= box[0];
	else if (x < 0.0)
		x += box[0];		
		
	if (y >= box[1])
		y -= box[1];
	else if (y < 0.0)
		y += box[1];			

	if (z >= box[2])
		z -= box[2];
	else if (z < 0.0)
		z += box[2];	
	
	return (floor(x/cellSize[0]) + floor(y/cellSize[1])*cellCountX +floor(z/cellSize[2])*cellCountXY);
}

#endif
