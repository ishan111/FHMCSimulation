#include "cellList.h"

cellList::cellList (std::vector < double > _box, double _cellSize, std::vector< atom* > _atoms) {
	box = _box;
		
	cellCountX = floor(_box[0]/_cellSize);
	cellCountY = floor(_box[1]/_cellSize);
	cellCountZ = floor(_box[2]/_cellSize);
	cellCountXY = cellCountX*cellCountY;
	cellCount = cellCountXY*cellCountZ;
	
	cellSize.resize(3);
	cellSize[0] = _box[0]/cellCountX;
	cellSize[1] = _box[1]/cellCountY;
	cellSize[2] = _box[2]/cellCountZ;

	// init neighbour list
	neighbours.assign(cellCount, std::vector < unsigned int > (26));
	initNeighbours();
	
	// init cell list
	cells.assign(cellCount, std::vector < atom* > (0));
	
	// clear all cells
	const unsigned int reserveCount = ceil(cellSize[0]*cellSize[1]*cellSize[2]);
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells[i].reserve(reserveCount);
		cells[i].clear();
	}
	
	if (_atoms.size() > 0) {
		//std::cout << "Sorting particles into cells." << std::endl;
		sortIntoCells(_atoms);
	}
}

void cellList::initNeighbours() {
	for (int k=0; k<cellCountX; k++) {
		for (int l=0; l<cellCountY; l++) {
			for (int m=0; m<cellCountZ; m++) {
				const unsigned int index = calcIndexS(k, l, m);
				neighbours[index].at(0) = calcIndexS(k+1, l+1, m+1);
				neighbours[index].at(1) = calcIndexS(k, l+1, m+1);
				neighbours[index].at(2) = calcIndexS(k-1, l+1, m+1);
				neighbours[index].at(3) = calcIndexS(k+1, l, m+1);
				neighbours[index].at(4) = calcIndexS(k, l, m+1);
				neighbours[index].at(5) = calcIndexS(k-1, l, m+1);
				neighbours[index].at(6) = calcIndexS(k+1, l-1, m+1);
				neighbours[index].at(7) = calcIndexS(k, l-1, m+1);
				neighbours[index].at(8) = calcIndexS(k-1, l-1, m+1);
				neighbours[index].at(9) = calcIndexS(k+1, l+1, m);
				neighbours[index].at(10) = calcIndexS(k, l+1, m);
				neighbours[index].at(11) = calcIndexS(k-1, l+1, m);
				neighbours[index].at(12) = calcIndexS(k+1, l, m);
				neighbours[index].at(13) = calcIndexS(k-1, l, m);
				neighbours[index].at(14) = calcIndexS(k+1, l-1, m);
				neighbours[index].at(15) = calcIndexS(k, l-1, m);
				neighbours[index].at(16) = calcIndexS(k-1, l-1, m);
				neighbours[index].at(17) = calcIndexS(k+1, l+1, m-1);
				neighbours[index].at(18) = calcIndexS(k, l+1, m-1);
				neighbours[index].at(19) = calcIndexS(k-1, l+1, m-1);
				neighbours[index].at(20) = calcIndexS(k+1, l, m-1);
				neighbours[index].at(21) = calcIndexS(k, l, m-1);
				neighbours[index].at(22) = calcIndexS(k-1, l, m-1);
				neighbours[index].at(23) = calcIndexS(k+1, l-1, m-1);
				neighbours[index].at(24) = calcIndexS(k, l-1, m-1);
				neighbours[index].at(25) = calcIndexS(k-1, l-1, m-1);
			}
		}
	}
}

void cellList::sortIntoCells (std::vector < atom* > _atoms) {
	// clear all cells
	for (unsigned int i = 0; i < cells.size(); i++)
		cells[i].clear();

	for (unsigned int i=0; i<_atoms.size(); i++) {
		const unsigned index = calcIndex(_atoms[i]->pos[0], _atoms[i]->pos[1], _atoms[i]->pos[2]);
		cells[index].push_back(_atoms[i]);
	}
}

void cellList::sortIntoCells (std::vector < atom >* _atoms) {
	// clear all cells
	for (unsigned int i = 0; i < cells.size(); i++)
		cells[i].clear();

	for (unsigned int i = 0; i < _atoms->size(); i++) {
		const unsigned index = calcIndex(_atoms->at(i).pos[0], _atoms->at(i).pos[1], _atoms->at(i).pos[2]);
		cells[index].push_back(&_atoms->at(i));
	}
}

void cellList::insertParticle (atom* _a) {
	const unsigned index = calcIndex(_a->pos[0], _a->pos[1], _a->pos[2]);
	cells[index].push_back(_a);
}


// swap values at addresses a and b, and then drop _a from cell list
void cellList::swapAndDeleteParticle (atom* _a, atom* _b) {
	const unsigned indexA = calcIndex(_a->pos[0], _a->pos[1], _a->pos[2]);
	const unsigned indexB = calcIndex(_b->pos[0], _b->pos[1], _b->pos[2]);
	
	unsigned int cellIndexA = 0, cellIndexB = 0;
	bool foundCellIndexA = false, foundCellIndexB = false;
	
	// locate position of atom _a in its cell
	for (unsigned int i = 0; i < cells[indexA].size(); i++) { // error?
		if (cells[indexA][i] == _a) {
			cellIndexA = i;
			foundCellIndexA = true;
			break;
		}
	}
	
	// locate position of atom _b in its cell
	for (unsigned int i = 0; i < cells[indexB].size(); i++) { // error ?
		if (cells[indexB][i] == _b) {
			cellIndexB = i;
			foundCellIndexB = true;
			break;
		}
	}
	
	if (!foundCellIndexA || !foundCellIndexB) {
		throw customException ("Failed to locate index in cell list properly");
	}
	
	// swap addresses
	cells[indexB][cellIndexB] = cells[indexA][cellIndexA];
	
	// remove _a from its cell
	cells[indexA].erase(cells[indexA].begin()+cellIndexA);
}

// translate a particle
void cellList::translateParticle (atom* _a, std::vector < double > _oldPos) {
	const unsigned indexOld = calcIndex(_oldPos[0], _oldPos[1], _oldPos[2]);
	const unsigned indexNew = calcIndex(_a->pos[0], _a->pos[1], _a->pos[2]);
	
	if (indexOld != indexNew) {
		unsigned int cellIndexOld = 0;
		bool foundCellIndexOld = false;
	
		// locate position of atom _a in its cell
		for (unsigned int i = 0; i < cells[indexOld].size(); i++) { //error?
			if (cells[indexOld][i] == _a) {
				cellIndexOld = i;
				foundCellIndexOld = true;
				break;
			}
		}
		
		if (!foundCellIndexOld) {
			throw customException ("Failed to locate cell index properly");
		}
		
		// remove _a from its cell
		cells[indexOld].erase(cells[indexOld].begin()+cellIndexOld);
		
		// insert _a into new cell
		cells[indexNew].push_back(_a);
	}
}