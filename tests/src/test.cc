#include "gtest-1.7.0/include/gtest/gtest.h"
#include <vector>
#include <iostream>
#include <limits>
#include "../../src/atom.h"
#include "../../src/bias.h"
#include "../../src/cellList.h"
#include "../../src/delete.h"
#include "../../src/global.h"
#include "../../src/insert.h"
#include "../../src/moves.h"
#include "../../src/potentials.h"
#include "../../src/swap.h"
#include "../../src/system.h"
#include "../../src/translate.h"
#include "../../src/utilities.h"

/* Test Atom Class to be 3D */
TEST (Initialize, Atom) {
	atom newAtom;
	EXPECT_EQ(newAtom.pos.size(), 3);
}

/* Test initializing moves */
TEST (Initialize, Moves) {
	moves test;
	int caught = 0;
	try {
		test.reportMoveStatistics(); // should fail since no moves added
	} catch (customException &ce) {
		caught = 1;
	}
	EXPECT_EQ (caught, 1);
}

/* Test initializing moves */
TEST (Initialize, MoveIndex) {
	moves test;
	const int ranIndex = 2;
	insertParticle insTest (ranIndex, "insert");
	EXPECT_EQ (insTest.whatType(), ranIndex);
	
	deleteParticle delTest (ranIndex, "insert");
	EXPECT_EQ (delTest.whatType(), ranIndex);
}

/* Test initializing moves */
TEST (Adding, Moves) {
	moves test;
	double prob1 = 0.4, prob2 = 0.6, tol = 1.0e-9;
	insertParticle insTest (0, "insert");
	test.addMove(&insTest, prob1);
	std::vector < double > pr = test.reportProbabilities();
	EXPECT_TRUE (fabs(pr[0] - 1.0) < tol);
	test.addMove(&insTest, prob2);
	std::vector < double > pr2 = test.reportProbabilities();
	EXPECT_TRUE (fabs(pr2[0] - prob1) < tol);
	EXPECT_TRUE (fabs(pr2[1] - (prob1+prob2)) < tol);
}

/* periodic boundary distances */
TEST (PBC, TwoVectors) {
	double tol = 1.0e-9, L = 10;
	std::vector < double > p1(3, 0), p2 (3, 1), p3 (3, L/2.0), p4 (3, 0.75*L), box(3, L);
	EXPECT_TRUE (fabs(pbc_dist2(p1, p2, box) - 3.0) < tol );
	EXPECT_TRUE (fabs(pbc_dist2(p1, p3, box) - (3.0*L*L/4.0)) < tol );
	EXPECT_TRUE (fabs(pbc_dist2(p1, p4, box) - (0.25*L*0.25*L*3.0)) < tol );
}

/* moving across periodic boundaries */
TEST (PBC, ReplaceInBox) {
	double tol = 1.0e-9, L = 10;
	std::vector < double > p1(3, 1.234*L), box(3, L);
	pbc(p1, box);
	for (unsigned int i = 0; i < p1.size(); ++i) {
		EXPECT_TRUE ( fabs(p1[i] - 0.234*L) < tol );
	}
}

/* Check potentials */
class lennardJonesTest : public ::testing::Test {
protected:
	pairPotential* lj;
	
	double eps, sigma, rcut, ushift, tol;
	std::vector < double > params;
	
	virtual void SetUp() {
		lj = new lennardJones;
		params.resize(5, 0);
		params[4] = 1;
		tol = 1.0e-9;
	}
};

TEST_F (lennardJonesTest, badParams0) {
	eps = -1.0;
	sigma = 1.0;
	rcut = 2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	bool caught = false;
	try {
		lj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete lj;
}

TEST_F (lennardJonesTest, badParams1) {
	eps = 1.0;
	sigma = -1.0;
	rcut = 2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	bool caught = false;
	try {
		lj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete lj;
}

TEST_F (lennardJonesTest, badParams2) {
	eps = 1.0;
	sigma = 1.0;
	rcut = -2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	bool caught = false;
	try {
		lj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete lj;
}

TEST_F (lennardJonesTest, StandardParams) {
	eps = 1.0;
	sigma = 1.0;
	rcut = 2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	lj->setParameters(params);
	atom a1, a2;
	std::vector < double > box(3, 2.1*rcut);
	
	a2.pos[2] = pow(2.0, 1./6.);	
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - -eps) < tol );
	
	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 0.9999*rcut;
	EXPECT_TRUE ( lj->energy(&a1, &a2, box) < 0.0 );
	
	a2.pos[2] = sigma;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 1.234;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - -0.812008901) < tol );
	
	delete lj;
}

TEST_F (lennardJonesTest, WCAparams) {
	eps = 1.0;
	sigma = 1.0;
	rcut = pow(2.0, 1./6.);
	ushift = 1.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	lj->setParameters(params);
	atom a1, a2;
	std::vector < double > box(3, 2.1*rcut);
	
	a2.pos[2] = pow(2.0, 1./6.);
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0.0) < tol );
	
	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 0.9999*rcut;
	EXPECT_TRUE ( lj->energy(&a1, &a2, box) > 0.0 );
	
	a2.pos[2] = 0.987;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 1.353386603) < tol );
	
	delete lj;
}

class squareWellTest : public ::testing::Test {
protected:
	squareWell *sw;
	double eps, sigma, width, tol;
	std::vector < double > params;
	
	virtual void SetUp() {
		sw = new squareWell;
		eps = 1.234;
		sigma = 2.345;
		width = 0.12;
		params.resize(4, 0);
		params[0] = sigma;
		params[1] = width;
		params[2] = eps;
		params[3] = 1; // Mtot
		tol = 1.0e-9;
	}
};

TEST_F (squareWellTest, badParams0) {
	params[0] = -1.0;
	bool caught = false;
	try {
		sw->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete sw;
}

TEST_F (squareWellTest, badParams1) {
	params[1] = -1.0;
	bool caught = false;
	try {
		sw->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete sw;
}

TEST_F (squareWellTest, badParams2) {
	params[2] = -1.0;
	bool caught = false;
	try {
		sw->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete sw;
}

TEST_F (squareWellTest, testInRangeLower) {	
	sw->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 1.0001*sigma;

	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - -eps) < tol);
	delete sw;
}

TEST_F (squareWellTest, testInRangeUpper) {
	sw->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 0.9999*(sigma+width);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - -eps) < tol);
	delete sw;
}

TEST_F (squareWellTest, testOutOfRange) {
	sw->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 1.0001*(sigma+width);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - 0) < tol);
	delete sw;
}

TEST_F (squareWellTest, testBadRange) {
	sw->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 0.9999*(sigma);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - NUM_INFINITY) < tol);
	delete sw;
}

TEST_F (squareWellTest, testRcut) {
	sw->setParameters(params);
	EXPECT_TRUE (fabs(sw->rcut() - (sigma+width)) < tol);
	delete sw;
}

TEST_F (squareWellTest, testTailCorrection) {
	sw->setParameters(params);
	EXPECT_TRUE (fabs(sw->tailCorrection(0.1234) - 0) < tol);
	delete sw;
}

class hardCoreTest : public ::testing::Test {
protected:
	hardCore *hc;
	double sigma, tol;
	std::vector < double > params;
	
	virtual void SetUp() {
		hc = new hardCore;
		sigma = 2.345;
		params.resize(2, 0);
		params[0] = sigma;
		params[1] = 1; // Mtot
		tol = 1.0e-9;
	}
};

TEST_F (hardCoreTest, badParams0) {
	params[0] = -1.0;
	bool caught = false;
	try {
		hc->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete hc;
}

TEST_F (hardCoreTest, testRange) {
	hc->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 1.0001*sigma;
	EXPECT_TRUE (fabs(hc->energy(&a1, &a2, box) - 0) < tol);
	delete hc;
}

TEST_F (hardCoreTest, testBadRange) {
	hc->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 0.9999*sigma;
	EXPECT_TRUE (fabs(hc->energy(&a1, &a2, box) - NUM_INFINITY) < tol);
	delete hc;
}

TEST_F (hardCoreTest, testRcut) {
	hc->setParameters(params);
	EXPECT_TRUE (fabs(hc->rcut() - sigma) < tol);
	delete hc;
}

TEST_F (hardCoreTest, testTailCorrection) {
	hc->setParameters(params);
	EXPECT_TRUE (fabs(hc->tailCorrection(0.1234) - 0) < tol);
	delete hc;
}

/* Check the System */
class InitializeSystem : public ::testing::Test {
protected:
	unsigned int nSpecies;
	double beta, tol;
	std::vector < double > box, mu;
	std::vector < int > maxSpecies, minSpecies;
	int tmmcSweepSize;

	virtual void SetUp() {
		nSpecies = 3;
		maxSpecies.resize(nSpecies), mu.resize(nSpecies), box.resize(3), minSpecies.resize(nSpecies);
		tol = 1.0e-9;
		beta = 1.234;
		tmmcSweepSize = 1;
		for (unsigned int i = 0; i < nSpecies; ++i) {
			box[i] = 2*i+1.0;
			mu[i] = 2*i+2.0;
			maxSpecies[i] = 10*i+10; // 10, 20, 30
			minSpecies[i] = 0;
		}
	}
};

TEST_F (InitializeSystem, nSpecies) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	EXPECT_EQ (mysys.nSpecies(), nSpecies);
	EXPECT_EQ (mysys.totNMax(), 60);
	EXPECT_EQ (mysys.totNMin(), 0);
}

TEST_F (InitializeSystem, beta) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	EXPECT_EQ (mysys.beta(), beta);
}

TEST_F (InitializeSystem, expandedStates) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 2);
	EXPECT_EQ (mysys.getTotalM(), 2);
	EXPECT_EQ (mysys.getCurrentM(), 0);
}

TEST_F (InitializeSystem, box) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	EXPECT_EQ (mysys.box(), box);
}

TEST_F (InitializeSystem, mu) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	for (unsigned int i = 0; i < mysys.nSpecies(); ++i) {
		EXPECT_EQ (mysys.mu(i), mu[i]);
	}
}

TEST_F (InitializeSystem, maxSpecies) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	for (unsigned int i = 0; i < mysys.nSpecies(); ++i) {
		EXPECT_EQ (mysys.maxSpecies(i), maxSpecies[i]);
	}
}

TEST_F (InitializeSystem, minSpecies) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	for (unsigned int i = 0; i < mysys.nSpecies(); ++i) {
		EXPECT_EQ (mysys.minSpecies(i), minSpecies[i]);
	}
}

TEST_F (InitializeSystem, setUpperWindowHigh) {
	std::vector <int> totNbounds (2, 0);
	totNbounds[0] = 0;
	totNbounds[1] = 20;
	
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	mysys.setTotNBounds (totNbounds);
	
	// should have resized species 3 from 30 to 20
	EXPECT_EQ (mysys.maxSpecies(2), totNbounds[1]);
	EXPECT_EQ (mysys.atoms[2].size(), totNbounds[1]);
	EXPECT_EQ (mysys.totNMax(), totNbounds[1]);
	EXPECT_EQ (mysys.totNMin(), totNbounds[0]);
}

TEST_F (InitializeSystem, setUpperWindowLower) {
	std::vector <int> totNbounds (2, 0);
	totNbounds[0] = 0;
	totNbounds[1] = 19;
	
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	mysys.setTotNBounds (totNbounds);
	
	// should have resized species 2,3 from 30,20 to 19,19
	EXPECT_EQ (mysys.maxSpecies(0), 10);
	EXPECT_EQ (mysys.maxSpecies(1), totNbounds[1]);
	EXPECT_EQ (mysys.maxSpecies(2), totNbounds[1]);
	EXPECT_EQ (mysys.atoms[1].size(), totNbounds[1]);
	EXPECT_EQ (mysys.atoms[2].size(), totNbounds[1]);
	EXPECT_EQ (mysys.totNMax(), totNbounds[1]);
	EXPECT_EQ (mysys.totNMin(), totNbounds[0]);
}

TEST_F (InitializeSystem, setInterfereWindowLower) {
	std::vector <int> totNbounds (2, 0);
	totNbounds[0] = 10;
	totNbounds[1] = 19;
	
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	mysys.setTotNBounds (totNbounds);
	
	// should have resized species 2,3 from 30,20 to 19,19
	EXPECT_EQ (mysys.maxSpecies(0), 10);
	EXPECT_EQ (mysys.maxSpecies(1), totNbounds[1]);
	EXPECT_EQ (mysys.maxSpecies(2), totNbounds[1]);
	EXPECT_EQ (mysys.atoms[1].size(), totNbounds[1]);
	EXPECT_EQ (mysys.atoms[2].size(), totNbounds[1]);
	EXPECT_EQ (mysys.totNMax(), totNbounds[1]);
	EXPECT_EQ (mysys.totNMin(), totNbounds[0]);
}

TEST_F (InitializeSystem, incrementEnergy) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	double dU = -0.1;
	mysys.incrementEnergy(dU);
	EXPECT_TRUE (fabs(mysys.energy() - dU) < tol);
	mysys.incrementEnergy(-dU);
	EXPECT_TRUE (fabs(mysys.energy() - 0) < tol);
}

TEST_F (InitializeSystem, addPotential) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	lennardJones ljtest;
	mysys.addPotential(0, 0, &ljtest);
	EXPECT_EQ (mysys.potentialIsSet (0, 0), true);
	EXPECT_EQ (mysys.potentialIsSet (0, 1), false);
	EXPECT_EQ (mysys.potentialIsSet (1, 1), false);
	EXPECT_EQ (mysys.potentialIsSet (1, 1), false);
}

TEST_F (InitializeSystem, insertAtom) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	atom typeA1, typeA2, typeB1, typeB2;
	std::vector < double > pos (3, 1.234), pos2 (3, 2.345);
	typeA1.pos = pos;
	typeA2.pos = pos2;
	
	const int nStart = mysys.getTotN();
	
	/* type A */
	EXPECT_TRUE(mysys.atoms[0].end() == (mysys.atoms[0].begin() + mysys.maxSpecies(0)));
	
	mysys.insertAtom (0, &typeA1);
	EXPECT_EQ(mysys.numSpecies[0], 1);
	
	mysys.insertAtom (0, &typeA2);
	EXPECT_EQ(mysys.numSpecies[0], 2);
	
	for (unsigned int i = 0; i < pos.size(); ++i) {
		EXPECT_TRUE(fabs(mysys.atoms[0][0].pos[i] - 1.234) < tol);
	}
	
	for (unsigned int i = 0; i < pos2.size(); ++i) {
		EXPECT_TRUE(fabs(mysys.atoms[0][1].pos[i] - 2.345) < tol);
	}
	
	/* type B */
	typeB1.pos = pos;
	typeB2.pos = pos2;
	
	EXPECT_TRUE(mysys.atoms[1].end() == (mysys.atoms[1].begin() + mysys.maxSpecies(1)));
	
	mysys.insertAtom (1, &typeB1);
	EXPECT_EQ(mysys.numSpecies[1], 1);
	
	mysys.insertAtom (1, &typeB2);
	EXPECT_EQ(mysys.numSpecies[1], 2);
	
	for (unsigned int i = 0; i < pos.size(); ++i) {
		EXPECT_TRUE(fabs(mysys.atoms[1][0].pos[i] - 1.234) < tol);
	}
	
	for (unsigned int i = 0; i < pos2.size(); ++i) {
		EXPECT_TRUE(fabs(mysys.atoms[1][1].pos[i] - 2.345) < tol);
	}
	
	EXPECT_EQ (mysys.getTotN(), nStart+4); // two of each type
}

TEST_F (InitializeSystem, deleteAtom) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	atom typeA1, typeA2, typeB1, typeB2;
	std::vector < double > pos (3, 1.234), pos2 (3, 2.345);
	typeA1.pos = pos;
	typeA2.pos = pos2;
	const int nStart = mysys.getTotN();
	
	mysys.insertAtom (0, &typeA1);
	mysys.insertAtom (0, &typeA2);
	mysys.deleteAtom (0, 0);	// should "replace" 0 with 1 information
	EXPECT_EQ (mysys.numSpecies[0], 1);
	for (unsigned int i = 0; i < pos2.size(); ++i) {
		EXPECT_TRUE(fabs(mysys.atoms[0][0].pos[i] - typeA2.pos[i]) < tol);
	}	
	mysys.deleteAtom (0, 0);	// now information technically should remain, but system should not record having this particle
	EXPECT_EQ (mysys.numSpecies[0], 0);
	for (unsigned int i = 0; i < pos2.size(); ++i) {
		EXPECT_TRUE(fabs(mysys.atoms[0][0].pos[i] - typeA2.pos[i]) < tol);
	}
	mysys.insertAtom (0, &typeA1);
	EXPECT_EQ (mysys.numSpecies[0], 1);
	for (unsigned int i = 0; i < pos2.size(); ++i) {
		EXPECT_TRUE(fabs(mysys.atoms[0][0].pos[i] - typeA1.pos[i]) < tol);
	}
	
	EXPECT_EQ (mysys.getTotN(), nStart+1); // total change of +1
}

TEST_F (InitializeSystem, scratchEnergy) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	lennardJones ljtest;
	std::vector < double > params (5), p1(3, 0.0), p2(3, 0.0);
	
	// wca
	params[0] = 1.0;
	params[1] = 1.0;
	params[2] = pow(2.0, 1./6.);
	params[3] = 1.0;
	params[4] = 1; // Mtot
	ljtest.setParameters(params);
	mysys.addPotential(0, 0, &ljtest);
	
	atom a, b;
	p2[2] = 1.01*pow(2.0, 1./6.);
	b.pos = p2;
	a.pos = p1;
	mysys.insertAtom(0, &a);
	mysys.insertAtom(0, &b);
	EXPECT_EQ (mysys.scratchEnergy(), 0.0);
		
	mysys.atoms[0][1].pos[2] = 0.99*pow(2.0, 1./6.);
	EXPECT_TRUE (mysys.scratchEnergy() > 0.0);

	// lj
	params[0] = 1.0;
	params[1] = 1.0;
	params[2] = 2.5;
	params[3] = 0.0;
	ljtest.setParameters(params);
	
	mysys.atoms[0][1].pos[2] = pow(2.0, 1./6.);
	EXPECT_TRUE (fabs(mysys.scratchEnergy() - -1.0) < tol);

	mysys.atoms[0][1].pos[2] = 1.0;
	EXPECT_TRUE (fabs(mysys.scratchEnergy() - 0.0) < tol);
}


TEST (testTMMC, tmmcGoodInit) {
	const int Nmax = 100, Nmin = 10, sweepSize = 100;
	std::vector < double > dummyBox (3, 0);
	bool pass = true;
	try {
		tmmc test (Nmax, Nmin, 1, sweepSize, dummyBox);
	} catch (customException &ce) {
		pass = false;
	}
	EXPECT_TRUE(pass);
}

TEST (testTMMC, tmmcBadInitNegative) {
	const int Nmax = 100, Nmin = -10, sweepSize = 100;
	std::vector < double > dummyBox (3, 0);
	bool pass = true;
	try {
		tmmc test (Nmax, Nmin, 1, sweepSize, dummyBox);
	} catch (customException &ce) {
		pass = false;
	}
	EXPECT_TRUE(!pass);
}

TEST (testTMMC, tmmcBadInitOrder) {
	const int Nmax = 10, Nmin = 100, sweepSize = 100;
	std::vector < double > dummyBox (3, 0);
	bool pass = true;
	try {
		tmmc test (Nmax, Nmin, 1, sweepSize, dummyBox);
	} catch (customException &ce) {
		pass = false;
	}
	EXPECT_TRUE(!pass);
}

class tmmBiasC : public ::testing::Test {
protected:
	tmmc* tmmcBias;
	int Nmin, Nmax, sweepSize;
	std::vector < double > dummyBox;
	
	virtual void SetUp() {
		Nmin = 3;
		Nmax = 10;
		dummyBox.resize(3, 0);
		sweepSize = 100;
		tmmcBias = new tmmc (Nmax, Nmin, 1, sweepSize, dummyBox);
	}
};

TEST_F (tmmBiasC, iterateForward) {
	int sweep1 = tmmcBias->numSweeps();
	tmmcBias->iterateForward();
	int sweep2 = tmmcBias->numSweeps();
	EXPECT_EQ (sweep1, 0);
	EXPECT_EQ (sweep2, 1);
}

TEST_F (tmmBiasC, cSize) {
	std::vector <double> collMat = tmmcBias->getC();
	EXPECT_TRUE (collMat.size() == 3*(Nmax-Nmin+1));
}

TEST_F (tmmBiasC, updateCMin) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmin, Nmin, 0, 0, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 1; i < collMat.size(); ++i) {
		EXPECT_EQ(collMat[i], 0);
	}
	EXPECT_TRUE(fabs(collMat[0] - 1.0) < 1.0e-9);
	
	delete tmmcBias;
}

TEST_F (tmmBiasC, updateCMax) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmax, Nmax, 0, 0, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 0; i < collMat.size(); ++i) {
		if (i != collMat.size()-3) {
			EXPECT_EQ(collMat[i], 0);
		}
	}
	EXPECT_TRUE(fabs(collMat[collMat.size()-3] - 1.0) < 1.0e-9);
	
	delete tmmcBias;
}

TEST_F (tmmBiasC, updateCForward) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmin, Nmin+1, 0, 0, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 2; i < collMat.size(); ++i) {
		EXPECT_EQ(collMat[i], 0);
	}
	EXPECT_TRUE(fabs(collMat[1] - pu) < 1.0e-9);
	EXPECT_TRUE(fabs(collMat[0] - (1-pu)) < 1.0e-9);
	
	delete tmmcBias;
}

TEST_F (tmmBiasC, updateCBackward) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmin+1, Nmin, 0, 0, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 0; i < collMat.size(); ++i) {
		if (i != 5 && i != 3) {
			EXPECT_EQ(collMat[i], 0);
		}
	}
	EXPECT_TRUE(fabs(collMat[5] - pu) < 1.0e-9);
	EXPECT_TRUE(fabs(collMat[3] - (1-pu)) < 1.0e-9);
	delete tmmcBias;
}

class tmmBiaslnPI : public ::testing::Test {
protected:
	tmmc* tmmcBias;
	int Nmin, Nmax, sweepSize;
	double pu;
	std::vector < double > dummyBox;
	
	virtual void SetUp() {
		pu = 0.123;
		Nmin = 3;
		Nmax = 5;
		dummyBox.resize(3, 0);
		sweepSize = 1;
		tmmcBias = new tmmc (Nmax, Nmin, 1, sweepSize, dummyBox);
		tmmcBias->updateC(Nmin, Nmin+1, 0, 0, pu);
		tmmcBias->updateC(Nmin+1, Nmin+2, 0, 0, pu);
		tmmcBias->updateC(Nmin+2, Nmin+1, 0, 0, pu);
	}
};

TEST_F (tmmBiaslnPI, incompleteC) {
	bool caught = false;
	try {
		tmmcBias->calculatePI();
	} catch (customException& ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	delete tmmcBias;
}

TEST_F (tmmBiaslnPI, completeC) {
	bool caught = false;
	tmmcBias->updateC(Nmin+1, Nmin, 0, 0, pu);
	try {
		tmmcBias->calculatePI();
	} catch (customException& ce) {
		caught = true;
	}
	EXPECT_TRUE(!caught);
	delete tmmcBias;
}

TEST_F (tmmBiaslnPI, checkVisitedNo) {
	bool pass = tmmcBias->checkFullyVisited();
	EXPECT_TRUE(!pass);
	delete tmmcBias;
}

TEST_F (tmmBiaslnPI, checkVisitedYes) {
	tmmcBias->updateC(Nmin+1, Nmin, 0, 0, pu);
	tmmcBias->updateC(Nmin, Nmin, 0, 0, pu); // also have to include visiting self
	tmmcBias->updateC(Nmin+1, Nmin+1, 0, 0, pu);
	tmmcBias->updateC(Nmin+2, Nmin+2, 0, 0, pu);	
	bool pass = tmmcBias->checkFullyVisited();
	EXPECT_TRUE(pass);
	delete tmmcBias;
}

TEST_F (tmmBiaslnPI, setBias) {
	std::vector < double > lnPIguess (3, 1.234);
	tmmcBias->setlnPI(lnPIguess);
	for (unsigned int i = 0; i < lnPIguess.size(); ++i) {
		EXPECT_TRUE (fabs(lnPIguess[i] - -tmmcBias->getBias(i)) < 1.0e-9);
	}
}

TEST_F (tmmBiaslnPI, checkPrintAndRead) {
	tmmcBias->updateC(Nmin+1, Nmin, 0, 0, pu);
	tmmcBias->calculatePI();
	tmmcBias->print("tmmBiaslnPI_checkPrint", true);
	std::vector < double > C1 = tmmcBias->getC();
#ifdef NETCDF_CAPABLE
	tmmcBias->readC("tmmBiaslnPI_checkPrint_C.nc");
#else
	tmmcBias->readC("tmmBiaslnPI_checkPrint_C.dat");
#endif
	std::vector < double > C2 = tmmcBias->getC();
	EXPECT_EQ (C2.size(), C1.size());
	for (unsigned int i = 0; i < C1.size(); ++i) {
		EXPECT_TRUE (fabs(C1[i] - C2[i]) < 1.0e-9);
	}
	
	delete tmmcBias;
}

TEST_F (tmmBiaslnPI, checkCAddresses) {
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin, 0, 0), 0);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin+1, 0, 0), 1);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin-1, 0, 0), 2);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+1, 0, 0), 3);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+2, 0, 0), 4);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin, 0, 0), 5);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+2, 0, 0), 6);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+3, 0, 0), 7);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+1, 0, 0), 8);
	
	delete tmmcBias;
}

TEST_F (tmmBiaslnPI, checkLnPIAddresses) {
	EXPECT_EQ (tmmcBias->getAddress(Nmin, 0), 0);
	EXPECT_EQ (tmmcBias->getAddress(Nmin+1, 0), 1);
	EXPECT_EQ (tmmcBias->getAddress(Nmin+2, 0), 2);
	
	delete tmmcBias;
}

TEST_F (tmmBiaslnPI, checkLnPI) {
	tmmcBias->updateC(Nmin+1, Nmin, 0, 0, pu);
	tmmcBias->calculatePI();
	
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin, 0)) - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin+1, 0)) - log(2.0)) < 1.0e-9);
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin+2, 0)) - (log(2.0)+log(0.5))) < 1.0e-9);
	
	delete tmmcBias;
}

TEST (testWALA, walaInitBadlnF) {
	const double s = 0.8, g = 0.5, lnF = -1.0;
	const int Nmax = 5, Nmin = 3;
	bool caught = false;
	std::vector < double > dummyBox (3, 0);
	
	try {
		wala (lnF, g, s, Nmax, Nmin, 1, dummyBox);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST (testWALA, walaInitBadGNeg) {
	const double s = 0.8, g = -0.5, lnF = 1.0;
	const int Nmax = 5, Nmin = 3;
	bool caught = false;
	std::vector < double > dummyBox (3, 0);
	try {
		wala (lnF, g, s, Nmax, Nmin, 1, dummyBox);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST (testWALA, walaInitBadGPos) {
	const double s = 0.8, g = 1.5, lnF = 1.0;
	const int Nmax = 5, Nmin = 3;
	bool caught = false;
	std::vector < double > dummyBox (3, 0);
	try {
		wala (lnF, g, s, Nmax, Nmin, 1, dummyBox);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST (testWALA, walaInitBadSNeg) {
	const double s = -0.8, g = 0.5, lnF = 1.0;
	const int Nmax = 5, Nmin = 3;
	bool caught = false;
	std::vector < double > dummyBox (3, 0);
	try {
		wala (lnF, g, s, Nmax, Nmin, 1, dummyBox);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST (testWALA, walaInitBadSPos) {
	const double s = 1.8, g = 0.5, lnF = 1.0;
	const int Nmax = 5, Nmin = 3;
	bool caught = false;
	std::vector < double > dummyBox (3, 0);
	try {
		wala (lnF, g, s, Nmax, Nmin, 1, dummyBox);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST (testWALA, walaInitSwitchedBounds) {
	const double s = 0.8, g = -0.5, lnF = 1.0;
	const int Nmax = 5, Nmin = 3;
	bool caught = false;
	std::vector < double > dummyBox (3, 0);
	try {
		wala (lnF, g, s, Nmin, Nmax, 1, dummyBox);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST (testWALA, walaInitNegLowerBound) {
	const double s = 0.8, g = -0.5, lnF = 1.0;
	const int Nmax = 5, Nmin = -3;
	bool caught = false;
	std::vector < double > dummyBox (3, 0);
	try {
		wala (lnF, g, s, Nmax, Nmin, 1, dummyBox);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

class testWalaBias : public ::testing::Test {
protected:
	wala* walaBias;
	int Nmin, Nmax;
	double s, g, lnF, pu;
	std::vector < double > dummyBox;
	
	virtual void SetUp() {
		pu = 0.123;
		Nmin = 3;
		Nmax = 5;
		s = 0.8;
		g = 0.5;
		lnF = 1.0;
		dummyBox.resize(3, 0);
		walaBias = new wala (lnF, g, s, Nmax, Nmin, 1, dummyBox);
	}
};

TEST_F (testWalaBias, testMatrixSizes) {
	std::vector < double > H = walaBias->getH(), lnPI = walaBias->getlnPI();
	EXPECT_EQ (H.size(), 3);
	EXPECT_EQ (lnPI.size(), 3);
}

TEST_F (testWalaBias, testUpdateSame) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin, 0);
	std::vector < double > H = walaBias->getH(), lnPI = walaBias->getlnPI();
	EXPECT_TRUE (fabs(H[0] - 3.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[0] - 3.0*lnF) < 1.0e-9);
}

TEST_F (testWalaBias, testUpdateDiff) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	std::vector < double > H = walaBias->getH(), lnPI = walaBias->getlnPI();
	EXPECT_TRUE (fabs(H[0] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(H[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(H[2] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[0] - lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[1] - lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[2] - lnF) < 1.0e-9);
}

TEST_F (testWalaBias, GetGoodAddress) {
	EXPECT_EQ (walaBias->getAddress(Nmin, 0), 0);
	EXPECT_EQ (walaBias->getAddress(Nmin+1, 0), 1);
	EXPECT_EQ (walaBias->getAddress(Nmin+2, 0), 2);
}

TEST_F (testWalaBias, GetBadAddress) {
	bool caught = false;
	try {
		EXPECT_EQ (walaBias->getAddress(Nmin+3, 0), 0);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST_F (testWalaBias, getBias) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin, 0)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin+1, 0)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin+2, 0)) - -lnF) < 1.0e-9);
}

TEST_F (testWalaBias, checkIterateForward) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	double lnF1 = 0, lnF2 = 0;
	std::vector < double > H1 (3, 0), H2 (3, 0);
	lnF1 = walaBias->lnF();
	H1 = walaBias->getH(); // ensure NOT originally all zero
	for (unsigned int i = 0; i < H1.size(); ++i) {
		EXPECT_TRUE (H1[i] > 0);
	}
	
	// this clears the H_ matrix and resets lnF
	walaBias->iterateForward();
	H2 = walaBias->getH();
	lnF2 = walaBias->lnF();
	for (unsigned int i = 0; i < H2.size(); ++i) {
		EXPECT_TRUE (fabs(H2[i] - 0) < 1.0e-9);
	}
	EXPECT_TRUE (fabs(lnF1*g - lnF2) < 1.0e-9);
}

TEST_F (testWalaBias, checkPrintReadlnPI) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	
	std::vector < double > lnPI_ref = walaBias->getlnPI();
	walaBias->print("walaBiaslnPI_checkPrint", false); // only print lnPI matrix
	
	// set lnPI to something random
	std::vector < double > lnPI_random (lnPI_ref.size(), 0.123456), lnPI_check (lnPI_ref.size(), 0);
	walaBias->setlnPI (lnPI_random);
	lnPI_check = walaBias->getlnPI();
	for (unsigned int i = 0; i < lnPI_check.size(); ++i) {
		EXPECT_TRUE (fabs(lnPI_check[i] - 0.123456) < 1.0e-9);
	}
	
	// read in and check again
#ifdef NETCDF_CAPABLE
	walaBias->readlnPI("walaBiaslnPI_checkPrint_lnPI.nc");
#else
	walaBias->readlnPI("walaBiaslnPI_checkPrint_lnPI.dat");
#endif
	
	std::vector < double > lnPI_new = walaBias->getlnPI();
	for (unsigned int i = 0; i < lnPI_new.size(); ++i) {
		EXPECT_TRUE (fabs(lnPI_new[i] - lnPI_ref[i]) < 1.0e-6); // read loses precision
	}
}

TEST_F (testWalaBias, checkEvaluateFlatnessNo) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	// average of [1, 2, 1] = 1.333, average * s (=0.8) = 1.06 > min (=1)
	bool flat = walaBias->evaluateFlatness();
	EXPECT_TRUE (!flat);
}

TEST_F (testWalaBias, checkEvaluateFlatnessYes) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	// average of [1, 1, 1] = 1, average * s (=0.8) = 0.8 < min (=1)
	bool flat = walaBias->evaluateFlatness();
	EXPECT_TRUE (flat);
}

TEST_F (InitializeSystem, setWALAbias) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	mysys.startWALA (1.0, 0.5, 0.8, 1);
	EXPECT_TRUE (mysys.useWALA);
}

TEST_F (InitializeSystem, unsetWALAbias) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	mysys.startWALA (1.0, 0.5, 0.8, 1);
	EXPECT_TRUE (mysys.useWALA);
	mysys.stopWALA();
	EXPECT_TRUE (!mysys.useWALA);
}

TEST_F (InitializeSystem, setTMMCbias) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	mysys.startTMMC (tmmcSweepSize, 1);
	EXPECT_TRUE (mysys.useTMMC);
}

TEST_F (InitializeSystem, unsetTMMCbias) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies, 1);
	mysys.startTMMC (tmmcSweepSize, 1);
	EXPECT_TRUE (mysys.useTMMC);
	mysys.stopTMMC();
	EXPECT_TRUE (!mysys.useTMMC);
}

class testComputeBias : public ::testing::Test {
protected:
	double s, g, lnF, pu;
	std::vector < int > specNmax, specNmin, bounds;
	int tmmcSweepSize;

	virtual void SetUp() {
		pu = 0.123;
		s = 0.8;
		g = 0.5;
		lnF = 1.0;
		tmmcSweepSize = 1;
		specNmin.resize(2, 0);
		specNmax.resize(2, 3);
		bounds.resize(2, 0);
		bounds[1] = 2;
	}
};

TEST_F (testComputeBias, setCalculateWALABias) {
	std::vector < double > ib (3, 10), mu (2, 1.0);
	simSystem mysys (2, 1.0, ib, mu, specNmax, specNmin, 1);
	mysys.setTotNBounds (bounds);
	mysys.startWALA (1.0, 0.5, 0.8, 1);
	EXPECT_TRUE (mysys.useWALA);
	
	// some simple updates
	mysys.getWALABias()->update(bounds[0], 0);
	mysys.getWALABias()->update(bounds[0]+1, 0);
	mysys.getWALABias()->update(bounds[0]+1, 0);
	mysys.getWALABias()->update(bounds[0]+2, 0);
	
	// WALA bias
	atom a1;
	
	// 0 atoms --> 1 atom
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 0) - exp(-1)) < 1.0e-9);
	mysys.insertAtom(0, &a1);
			
	// 1 atom --> 2 atoms (same type)
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+2, 0) - exp(1)) < 1.0e-9);
	
	// delete, then move to 1 atom total --> +1 atom of another type
	mysys.deleteAtom (0, 0);
	mysys.insertAtom(1, &a1);
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+2, 0) - exp(1)) < 1.0e-9); // should be same result
}

TEST_F (testComputeBias, setCalculateTMMCBias) {
	std::vector < double > ib (3, 10), mu (2, 1.0);
	simSystem mysys (2, 1.0, ib, mu, specNmax, specNmin, 1);
	mysys.setTotNBounds (bounds);
	mysys.startTMMC (tmmcSweepSize, 1);
	EXPECT_TRUE (mysys.useTMMC);
	
	// some simple updates
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0], 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0]+1, 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+2, 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0], 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+2, bounds[0]+2, 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+2, bounds[0]+1, 0, 0, std::min(1.0, pu));
	EXPECT_TRUE (mysys.getTMMCBias()->checkFullyVisited());
	mysys.getTMMCBias()->calculatePI();
	
	// TMMC bias
	atom a1;
	
	// 0 atoms --> 1 atom
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 0) - 2.0/3.0) < 1.0e-9); // 0 assigned as a references
	mysys.insertAtom(0, &a1);
			
	// 1 atom --> 2 atoms (same type)
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+2, 0) - 3.0/2.0) < 1.0e-9);
	
	// delete, then move to 1 atom total --> +1 atom of another type
	mysys.deleteAtom (0, 0);
	mysys.insertAtom(1, &a1);
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+2, 0) - 3.0/2.0) < 1.0e-9); // should be same result	
}

TEST_F (testComputeBias, testInSituWALASingleComponent) {
	std::vector < double > mu (1, std::numeric_limits<double>::max()); // force an insertion to an empty system
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (1, 3), nmin (1, 0);
	lnF = 3.1; // something random this time
	simSystem mysys (1, 1.0, ib, mu, nmax, nmin, 1);
	mysys.startWALA (lnF, 0.5, 0.8, 1);
	EXPECT_TRUE (mysys.useWALA);
		
	hardCore hc;
	std::vector < double > params (2, 1.0);
	hc.setParameters (params);
	mysys.addPotential (0, 0, &hc, false);
	
	moves usedMoves;
	insertParticle newIns (0, "insert");
	usedMoves.addMove(&newIns, 1.0);
	
	// will insert
	usedMoves.makeMove(mysys);
	
	// check WALA properties - should have incremented where the system ENDED (at N = 1)
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(1) - -lnF) < 1.0e-9);
	
	// all the rest should be 0
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(0) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(2) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(3) - 0) < 1.0e-9);
}

TEST_F (testComputeBias, testInSituWALAMultiComponent) {
	std::vector < double > mu (2, std::numeric_limits<double>::max()); // force an insertion to an empty system
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (2, 3), nmin (2, 0);
	lnF = 3.1; // something random this time
	simSystem mysys (2, 1.0, ib, mu, nmax, nmin, 1);
	mysys.startWALA (lnF, 0.5, 0.8, 1);
	EXPECT_TRUE (mysys.useWALA);
	
	hardCore hc11, hc12, hc22;
	std::vector < double > params (2, 1.0);
	hc11.setParameters (params);
	params[0] = 0.0; // ergo 1 and 2 can sit on top of each other
	hc12.setParameters (params);
	params[1] = 2.0;
	hc22.setParameters (params);
	mysys.addPotential (0, 0, &hc11, false);
	mysys.addPotential (0, 1, &hc12, false);
	mysys.addPotential (1, 1, &hc22, false);
		
	moves usedMoves;
	insertParticle newIns (0, "insert");
	usedMoves.addMove(&newIns, 1.0);
	
	// will insert first species
	usedMoves.makeMove(mysys);
	
	// check WALA properties - should have incremented where the system ENDED (at N = 1)
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(1) - -lnF) < 1.0e-9);
	
	// all the rest should be 0
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(0) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(2) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(3) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(4) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(5) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(6) - 0) < 1.0e-9);
	
	// will insert the second species
	moves usedMoves2;
	insertParticle newIns2 (1, "insert");
	usedMoves2.addMove(&newIns2, 1.0);
	
	usedMoves2.makeMove(mysys);
	
	// check WALA properties - should have incremented where the system ENDED (at N_tot = 2)
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(2) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(1) - -lnF) < 1.0e-9); // left over from first insertion
	
	// all the rest should be 0
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(0) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(3) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(4) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(5) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(6) - 0) < 1.0e-9);
}

TEST_F (testComputeBias, testInSituTMMCSingleComponent) {
	std::vector < double > mu (1, std::numeric_limits<double>::max()); // force an insertion to an empty system
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (1, 3), nmin (1, 0);
	simSystem mysys (1, 1.0, ib, mu, nmax, nmin, 1);
	mysys.startTMMC (tmmcSweepSize, 1);
	EXPECT_TRUE (mysys.useTMMC);
		
	hardCore hc;
	std::vector < double > params (2, 1.0);
	hc.setParameters (params);
	mysys.addPotential (0, 0, &hc, false);
	
	moves usedMoves;
	insertParticle newIns (0, "insert");
	usedMoves.addMove(&newIns, 1.0);
	
	// will insert
	usedMoves.makeMove(mysys);
	
	// check TMMC properties - should have incremented where the system started from (N = 0)
	std::vector < double > C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 2; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// insert again
	usedMoves.makeMove(mysys);
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 6; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
}

TEST_F (testComputeBias, testInSituTMMCMultiComponent) {
	std::vector < double > mu (2, std::numeric_limits<double>::max()); // force an insertion to an empty system
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (2, 3), nmin (2, 0);
	simSystem mysys (2, 1.0, ib, mu, nmax, nmin, 1);
	mysys.startTMMC (tmmcSweepSize, 1);
	EXPECT_TRUE (mysys.useTMMC);
		
	hardCore hc11, hc12, hc22;
	std::vector < double > params (2, 1.0);
	hc11.setParameters (params);
	params[0] = 0.0;
	hc12.setParameters (params);
	params[0] = 0.0;
	hc22.setParameters (params);
	mysys.addPotential (0, 0, &hc11, false);
	mysys.addPotential (0, 1, &hc12, false);
	mysys.addPotential (1, 1, &hc22, false);
	
	moves usedMoves;
	insertParticle newIns (0, "insert");
	usedMoves.addMove(&newIns, 1.0);
	
	// will insert first species
	usedMoves.makeMove(mysys);
	
	// check TMMC properties - should have incremented where the system started from (N = 0)
	std::vector < double > C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 2; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// insert same species again
	usedMoves.makeMove(mysys);
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 6; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// do insertions with the other species
	moves usedMoves2;
	insertParticle newIns2 (1, "insert");
	usedMoves2.addMove(&newIns2, 1.0);
	
	// insert 1 atom from second species
	usedMoves2.makeMove(mysys);
	
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[6] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[7] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[8] - 0.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 9; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
}

TEST (testSwapMove, twoComponents) {
	std::vector < double > mu (3, 0); 
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (3, 3), nmin (3, 0);
	simSystem mysys (3, 1.0, ib, mu, nmax, nmin, 1);

	squareWell sw12, sw13, sw23; // only bother to set cross interactions here
	std::vector < double > params (4, 0.0);
	params[0] = 1.0; params[1] = 0.5; params[2] = 1.0, params[3] = 1;
	sw12.setParameters (params);
	params[0] = 1.0; params[1] = 0.5; params[2] = 0.3, params[3] = 1;
	sw13.setParameters (params);
	params[0] = 1.0; params[1] = 0.5; params[2] = 2.0, params[3] = 1;
	sw23.setParameters (params);
	mysys.addPotential (0, 1, &sw12, true);
	mysys.addPotential (0, 2, &sw13, true);
	mysys.addPotential (1, 2, &sw23, true);
	
	atom a1, a2, a3;
	a1.pos[0] = 2.5;
	a1.pos[1] = 5;
	a1.pos[2] = 5;
	
	a2.pos[0] = 3.75;
	a2.pos[1] = 5;
	a2.pos[2] = 5;
	
	a3.pos[0] = 7.5;
	a3.pos[1] = 5;
	a3.pos[2] = 7.5;
	
	mysys.insertAtom(0, &a1);
	mysys.insertAtom(1, &a2);
	mysys.insertAtom(2, &a3);
		
	moves usedMoves;
	swapParticles newSwap (0, 1, "swap"); // only can swap a1 and a2
	usedMoves.addMove(&newSwap, 1.0);
	
	// swap a1 & a2 should result in no change since they are an isolated pair
	bool noMove = true;
	double lastAns = 0, U_save = mysys.scratchEnergy();
	while (noMove) {
		usedMoves.makeMove (mysys);
		std::vector < double > ans = usedMoves.reportMoveStatistics();
		if (ans[0] > lastAns) {
			noMove = false;
		} 
		lastAns = ans[0];
	}
	double dU1 = mysys.scratchEnergy();
	EXPECT_TRUE (fabs(U_save - dU1) < 1.0e-9);

	// now try to swap a1 and a3 after swapping a1 and a2
	moves usedMoves2;
	swapParticles newSwap2 (1, 2, "swap"); // swap a2 and a3
	usedMoves2.addMove(&newSwap2, 1.0);
	noMove = true;
	lastAns = 0;
	U_save = mysys.scratchEnergy();
	while (noMove) {
		usedMoves2.makeMove (mysys);
		std::vector < double > ans = usedMoves2.reportMoveStatistics();
		if (ans[0] > lastAns) {
			noMove = false;
		} 
		lastAns = ans[0];
	}	
	EXPECT_TRUE (fabs(mysys.scratchEnergy() - -0.3) < 1.0e-9); // only 1-3 interaction remains
		
	// swap a1 and a3 which are now an isolated pair
	moves usedMoves3;
	swapParticles newSwap3 (0, 2, "swap"); // only can swap a1 and a3
	usedMoves3.addMove(&newSwap3, 1.0);
	noMove = true;
	lastAns = 0;
	U_save = mysys.scratchEnergy();
	while (noMove) {
		usedMoves3.makeMove (mysys);
		std::vector < double > ans = usedMoves3.reportMoveStatistics();
		if (ans[0] > lastAns) {
			noMove = false;
		} 
		lastAns = ans[0];
	}	
	EXPECT_TRUE (fabs(U_save - mysys.scratchEnergy() - 0) < 1.0e-9); // no change in energy
}

/* Expanded Ensemble Tests */
TEST (testExpandedTMMC, tmmcBadInitOrder) {
	const int Nmax = 10, Nmin = 100, sweepSize = 100, Mtot = 3;
	std::vector < double > dummyBox (3, 0);
	bool pass = true;
	try {
		tmmc test (Nmax, Nmin, Mtot, sweepSize, dummyBox);
	} catch (customException &ce) {
		pass = false;
	}
	EXPECT_TRUE(!pass);
}

class tmmcExpandedBiasC : public ::testing::Test {
protected:
	tmmc* tmmcBias;
	int Nmin, Nmax, sweepSize, Mtot;
	std::vector < double > dummyBox;
	
	virtual void SetUp() {
		Nmin = 3;
		Nmax = 10;
		dummyBox.resize(3, 0);
		sweepSize = 100;
		Mtot = 3;
		tmmcBias = new tmmc (Nmax, Nmin, Mtot, sweepSize, dummyBox);
	}
};

TEST_F (tmmcExpandedBiasC, iterateForward) {
	int sweep1 = tmmcBias->numSweeps();
	tmmcBias->iterateForward();
	int sweep2 = tmmcBias->numSweeps();
	EXPECT_EQ (sweep1, 0);
	EXPECT_EQ (sweep2, 1);
}

TEST_F (tmmcExpandedBiasC, cSize) {
	std::vector <double> collMat = tmmcBias->getC();
	EXPECT_TRUE (collMat.size() == 3*(Nmax-Nmin+1)*Mtot);
}

TEST_F (tmmcExpandedBiasC, updateCMinNoMove) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmin, Nmin, 1, 1, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 0; i < collMat.size(); ++i) {
		if (i == 3) {
			continue;
		} else {
			EXPECT_EQ(collMat[i], 0);
		}
	}
	EXPECT_TRUE(fabs(collMat[3] - 1.0) < 1.0e-9);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, updateCMinWithMove) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmin, Nmin, 1, 2, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 0; i < collMat.size(); ++i) {
		if (i == 3 || i == 4) {
			continue;
		} else {
			EXPECT_EQ(collMat[i], 0);
		}
	}
	EXPECT_TRUE(fabs(collMat[3] - (1-pu)) < 1.0e-9);
	EXPECT_TRUE(fabs(collMat[4] - pu) < 1.0e-9);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, updateCMaxNoMove) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmax, Nmax, 1, 1, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 0; i < collMat.size(); ++i) {
		if (i != 3*((10-3)*3 + 1)) {
			EXPECT_EQ(collMat[i], 0);
		}
	}
	EXPECT_TRUE(fabs(collMat[3*((10-3)*3 + 1)] - 1.0) < 1.0e-9);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, updateCMaxWithMove) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmax, Nmax, 1, 2, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 0; i < collMat.size(); ++i) {
		if (i != 3*((10-3)*3 + 1) && i != 3*((10-3)*3 + 1)+1) {
			EXPECT_EQ(collMat[i], 0);
		}
	}
	EXPECT_TRUE(fabs(collMat[3*((10-3)*3 + 1)] - (1-pu)) < 1.0e-9);
	EXPECT_TRUE(fabs(collMat[3*((10-3)*3 + 1)+1] - pu) < 1.0e-9);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, updateCForward) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmin, Nmin+1, 2, 0, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 0; i < collMat.size(); ++i) {
		if (i != 7 && i != 6) {
			EXPECT_EQ(collMat[i], 0);
		}
	}
	EXPECT_TRUE(fabs(collMat[7] - pu) < 1.0e-9);
	EXPECT_TRUE(fabs(collMat[6] - (1-pu)) < 1.0e-9);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, updateCBackward) {
	const double pu = 0.123;
	tmmcBias->updateC(Nmin+1, Nmin, 0, 2, pu);
	std::vector <double> collMat = tmmcBias->getC();
	
	for (unsigned int i = 0; i < collMat.size(); ++i) {
		if (i != 9 && i != 11) {
			EXPECT_EQ(collMat[i], 0);
		}
	}
	EXPECT_TRUE(fabs(collMat[11] - pu) < 1.0e-9);
	EXPECT_TRUE(fabs(collMat[9] - (1-pu)) < 1.0e-9);
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchHighM) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin, Nmin+1, 2, 3, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchLowM) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin+1, Nmin, 0, -1, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchBadM1) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin+1, Nmin, 0, 0, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchBadM2) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin, Nmin+1, 0, 2, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchBadM3) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin, Nmin+1, 1, 2, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchBadM4) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin, Nmin+1, 2, 1, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchBadM5) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin+1, Nmin, 1, 1, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchBadM6) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin+1, Nmin, 1, 0, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchBadM7) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin+1, Nmin, 2, 0, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, catchBadM8) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin+1, Nmin, 2, 1, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, allowLowBoundaryN) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmin, Nmin-1, 0, 2, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(!caught);
	
	delete tmmcBias;
}

TEST_F (tmmcExpandedBiasC, allowHighBoundaryN) {
	const double pu = 0.123;
	bool caught = false;
	try {
		tmmcBias->updateC(Nmax, Nmax+1, 2, 0, pu);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE(!caught);
	
	delete tmmcBias;
}

class tmmBiasExpandedlnPI : public ::testing::Test {
protected:
	tmmc* tmmcBias;
	int Nmin, Nmax, sweepSize, Mtot;
	double pu;
	std::vector < double > dummyBox;
	
	virtual void SetUp() {
		pu = 0.123;
		Nmin = 3;
		Nmax = 5;
		dummyBox.resize(3, 0);
		sweepSize = 1;
		Mtot = 3;
		tmmcBias = new tmmc (Nmax, Nmin, Mtot, sweepSize, dummyBox);
		
		tmmcBias->updateC(Nmin, Nmin, 0, 0, pu);
		tmmcBias->updateC(Nmin, Nmin, 0, 1, pu);

		tmmcBias->updateC(Nmin, Nmin, 1, 1, pu);
		tmmcBias->updateC(Nmin, Nmin, 1, 2, pu);
		tmmcBias->updateC(Nmin, Nmin, 1, 0, pu);
		
		tmmcBias->updateC(Nmin, Nmin, 2, 2, pu);
		tmmcBias->updateC(Nmin, Nmin+1, 2, 0, pu);
		tmmcBias->updateC(Nmin, Nmin, 2, 1, pu);

		tmmcBias->updateC(Nmin+1, Nmin+1, 0, 0, pu);
		tmmcBias->updateC(Nmin+1, Nmin+1, 0, 1, pu);
		tmmcBias->updateC(Nmin+1, Nmin, 0, 2, pu);
		
		tmmcBias->updateC(Nmin+1, Nmin+1, 1, 1, pu);
		tmmcBias->updateC(Nmin+1, Nmin+1, 1, 2, pu);
		tmmcBias->updateC(Nmin+1, Nmin+1, 1, 0, pu);
		
		tmmcBias->updateC(Nmin+1, Nmin+1, 2, 2, pu);
		tmmcBias->updateC(Nmin+1, Nmin+2, 2, 0, pu);
		tmmcBias->updateC(Nmin+1, Nmin+1, 2, 1, pu);
		
		// "missing states"
		/*tmmcBias->updateC(Nmin+2, Nmin+2, 0, 0, pu);
		tmmcBias->updateC(Nmin+2, Nmin+1, 0, 2, pu);*/
	}
};

TEST_F (tmmBiasExpandedlnPI, checkVisitedNo) {
	bool pass = tmmcBias->checkFullyVisited();
	EXPECT_TRUE(!pass);
	delete tmmcBias;
}

TEST_F (tmmBiasExpandedlnPI, checkVisitedYes) {	
	tmmcBias->updateC(Nmin+2, Nmin+2, 0, 0, pu);
	tmmcBias->updateC(Nmin+2, Nmin+1, 0, 2, pu);
	bool pass = tmmcBias->checkFullyVisited();
	EXPECT_TRUE(pass);
	delete tmmcBias;
}

TEST_F (tmmBiasExpandedlnPI, incompleteC) {
	bool caught = false;
	try {
		tmmcBias->calculatePI();
	} catch (customException& ce) {
		caught = true;
	}
	EXPECT_TRUE(caught);
	delete tmmcBias;
}

TEST_F (tmmBiasExpandedlnPI, completeC) {
	bool caught = false;
	tmmcBias->updateC(Nmin+2, Nmin+2, 0, 0, pu);
	tmmcBias->updateC(Nmin+2, Nmin+1, 0, 2, pu);
	try {
		tmmcBias->calculatePI();
	} catch (customException& ce) {
		caught = true;
	}
	EXPECT_TRUE(!caught);
	delete tmmcBias;
}

TEST_F (tmmBiasExpandedlnPI, setBias) {
	std::vector < double > lnPIguess (9, 1.234);
	tmmcBias->setlnPI(lnPIguess);
	for (unsigned int i = 0; i < lnPIguess.size(); ++i) {
		EXPECT_TRUE (fabs(lnPIguess[i] - -tmmcBias->getBias(i)) < 1.0e-9);
	}
}

TEST_F (tmmBiasExpandedlnPI, checkPrintAndRead) {
	tmmcBias->updateC(Nmin+2, Nmin+2, 0, 0, pu);	// add missing states
	tmmcBias->updateC(Nmin+2, Nmin+1, 0, 2, pu);
	tmmcBias->calculatePI();
	tmmcBias->print("tmmBiasExpandedlnPI_checkPrint", true);
	std::vector < double > C1 = tmmcBias->getC();
#ifdef NETCDF_CAPABLE
	tmmcBias->readC("tmmBiasExpandedlnPI_checkPrint_C.nc");
#else
	tmmcBias->readC("tmmBiasExpandedlnPI_checkPrint_C.dat");
#endif
	std::vector < double > C2 = tmmcBias->getC();
	EXPECT_EQ (C2.size(), C1.size());
	for (unsigned int i = 0; i < C1.size(); ++i) {
		EXPECT_TRUE (fabs(C1[i] - C2[i]) < 1.0e-9);
	}
	
	delete tmmcBias;
}

TEST_F (tmmBiasExpandedlnPI, checkCAddresses) {
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin, 0, 0), 0);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin, 0, 1), 1);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin-1, 0, 2), 2);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin, 1, 1), 3);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin, 1, 2), 4);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin, 1, 0), 5);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin, 2, 2), 6);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin+1, 2, 0), 7);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin, Nmin, 2, 1), 8);
	
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+1, 0, 0), 9);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+1, 0, 1), 10);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin, 0, 2), 11);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+1, 1, 1), 12);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+1, 1, 2), 13);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+1, 1, 0), 14);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+1, 2, 2), 15);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+2, 2, 0), 16);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+1, Nmin+1, 2, 1), 17);
	
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+2, 0, 0), 18);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+2, 0, 1), 19);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+1, 0, 2), 20);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+2, 1, 1), 21);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+2, 1, 2), 22);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+2, 1, 0), 23);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+2, 2, 2), 24);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+3, 2, 0), 25);
	EXPECT_EQ (tmmcBias->getTransitionAddress(Nmin+2, Nmin+2, 2, 1), 26);
	
	delete tmmcBias;
}

TEST_F (tmmBiasExpandedlnPI, checkLnPIAddresses) {
	EXPECT_EQ (tmmcBias->getAddress(Nmin, 0), 0);
	EXPECT_EQ (tmmcBias->getAddress(Nmin, 1), 1);
	EXPECT_EQ (tmmcBias->getAddress(Nmin, 2), 2);
	EXPECT_EQ (tmmcBias->getAddress(Nmin+1, 0), 3);
	EXPECT_EQ (tmmcBias->getAddress(Nmin+1, 1), 4);
	EXPECT_EQ (tmmcBias->getAddress(Nmin+1, 2), 5);
	EXPECT_EQ (tmmcBias->getAddress(Nmin+2, 0), 6);
	
	// these will never be used, but are still present in the matrix
	EXPECT_EQ (tmmcBias->getAddress(Nmin+2, 1), 7);
	EXPECT_EQ (tmmcBias->getAddress(Nmin+2, 2), 8);
	
	delete tmmcBias;
}


TEST_F (tmmBiasExpandedlnPI, checkLnPI) {
	tmmcBias->updateC(Nmin+2, Nmin+2, 0, 0, pu);
	tmmcBias->updateC(Nmin+2, Nmin+1, 0, 2, pu);
	tmmcBias->calculatePI();
	
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin, 0)) - 0.0) < 1.0e-9); // reference
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin, 1)) - (log(3/2.))) < 1.0e-9); 
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin, 2)) - (log(3/2.))) < 1.0e-9);
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin+1, 0)) - (log(3/2.))) < 1.0e-9);
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin+1, 1)) - (log(3/2.))) < 1.0e-9);
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin+1, 2)) - (log(3/2.))) < 1.0e-9);
	EXPECT_TRUE (fabs(-tmmcBias->getBias(tmmcBias->getAddress(Nmin+2, 0)) - 0.0) < 1.0e-9);

	delete tmmcBias;
}

class testExpandedWalaBias : public ::testing::Test {
protected:
	wala* walaBias;
	int Nmin, Nmax, Mtot;
	double s, g, lnF, pu;
	std::vector < double > dummyBox;
	
	virtual void SetUp() {
		pu = 0.123;
		Nmin = 3;
		Nmax = 5;
		s = 0.8;
		g = 0.5;
		lnF = 1.0;
		Mtot = 3;
		dummyBox.resize(3, 0);
		walaBias = new wala (lnF, g, s, Nmax, Nmin, Mtot, dummyBox);
	}
};

TEST_F (testExpandedWalaBias, testMatrixSizes) {
	std::vector < double > H = walaBias->getH(), lnPI = walaBias->getlnPI();
	EXPECT_EQ (H.size(), (Nmax-Nmin+1)*Mtot);
	EXPECT_EQ (lnPI.size(), (Nmax-Nmin+1)*Mtot);
}

TEST_F (testExpandedWalaBias, testUpdateSame) {
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin, 1);
	std::vector < double > H = walaBias->getH(), lnPI = walaBias->getlnPI();
	EXPECT_TRUE (fabs(H[1] - 3.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[1] - 3.0*lnF) < 1.0e-9);
}

TEST_F (testExpandedWalaBias, testUpdateDiff) {
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin+1, 2);
	walaBias->update(Nmin+2, 0);
	std::vector < double > H = walaBias->getH(), lnPI = walaBias->getlnPI();
	EXPECT_TRUE (fabs(H[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(H[5] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(H[6] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[0] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[1] - lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[2] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[4] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[5] - lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[6] - lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[7] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[8] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[9] - 0.0) < 1.0e-9);
}

TEST_F (testExpandedWalaBias, GetGoodAddress) {
	EXPECT_EQ (walaBias->getAddress(Nmin, 0), 0);
	EXPECT_EQ (walaBias->getAddress(Nmin+1, 0), 3);
	EXPECT_EQ (walaBias->getAddress(Nmin+2, 0), 6);
}

TEST_F (testExpandedWalaBias, GetBadAddress) {
	bool caught = false;
	try {
		EXPECT_EQ (walaBias->getAddress(Nmin+3, 0), 0);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST_F (testExpandedWalaBias, getBias) {
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin+1, 2);
	walaBias->update(Nmin+2, 0);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin, 1)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin+1, 2)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin+2, 0)) - -lnF) < 1.0e-9);
}

TEST_F (testExpandedWalaBias, checkIterateForward) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin, 2);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+1, 1);
	walaBias->update(Nmin+1, 2);
	walaBias->update(Nmin+2, 0);
	walaBias->update(Nmin+2, 1);
	walaBias->update(Nmin+2, 2);
	double lnF1 = 0, lnF2 = 0;
	std::vector < double > H1 (9, 0), H2 (9, 0);
	lnF1 = walaBias->lnF();
	H1 = walaBias->getH(); // ensure NOT originally all zero
	for (unsigned int i = 0; i < H1.size(); ++i) {
		EXPECT_TRUE (H1[i] > 0);
	}
	
	// this clears the H_ matrix and resets lnF
	walaBias->iterateForward();
	H2 = walaBias->getH();
	lnF2 = walaBias->lnF();
	for (unsigned int i = 0; i < H2.size(); ++i) {
		EXPECT_TRUE (fabs(H2[i] - 0) < 1.0e-9);
	}
	EXPECT_TRUE (fabs(lnF1*g - lnF2) < 1.0e-9);
}

TEST_F (testExpandedWalaBias, checkPrintReadlnPI) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin, 2);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+1, 1);
	walaBias->update(Nmin+1, 2);
	walaBias->update(Nmin+2, 0);
	walaBias->update(Nmin+2, 1);
	walaBias->update(Nmin+2, 2);
	
	std::vector < double > lnPI_ref = walaBias->getlnPI();
	walaBias->print("walaBiasExpandedlnPI_checkPrint", false); // only print lnPI matrix
	
	// set lnPI to something random
	std::vector < double > lnPI_random (lnPI_ref.size(), 0.123456), lnPI_check (lnPI_ref.size(), 0);
	walaBias->setlnPI (lnPI_random);
	lnPI_check = walaBias->getlnPI();
	for (unsigned int i = 0; i < lnPI_check.size(); ++i) {
		EXPECT_TRUE (fabs(lnPI_check[i] - 0.123456) < 1.0e-9);
	}
	
	// read in and check again
#ifdef NETCDF_CAPABLE
	walaBias->readlnPI("walaBiasExpandedlnPI_checkPrint_lnPI.nc");
#else
	walaBias->readlnPI("walaBiasExpandedlnPI_checkPrint_lnPI.dat");
#endif
	
	std::vector < double > lnPI_new = walaBias->getlnPI();
	for (unsigned int i = 0; i < lnPI_new.size(); ++i) {
		EXPECT_TRUE (fabs(lnPI_new[i] - lnPI_ref[i]) < 1.0e-6); // read loses precision
	}
}

TEST_F (testExpandedWalaBias, checkEvaluateFlatnessNo) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin, 2);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+1, 1);
	walaBias->update(Nmin+1, 1);
	walaBias->update(Nmin+1, 1);
	walaBias->update(Nmin+1, 1);
	walaBias->update(Nmin+1, 2);
	walaBias->update(Nmin+2, 0);
	walaBias->update(Nmin+2, 1);
	walaBias->update(Nmin+2, 2);
	// average = 8*1+1*4 / 9 = 1.333, average * s (=0.8) = 1.06 > min (=1)
	bool flat = walaBias->evaluateFlatness();
	EXPECT_TRUE (!flat);
}

TEST_F (testExpandedWalaBias, checkEvaluateFlatnessYes) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin, 2);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+1, 1);
	walaBias->update(Nmin+1, 2);
	walaBias->update(Nmin+2, 0);
	walaBias->update(Nmin+2, 1);
	walaBias->update(Nmin+2, 2);
	// average = 1, average * s (=0.8) = 0.8 < min (=1)
	bool flat = walaBias->evaluateFlatness();
	EXPECT_TRUE (flat);
}

/* Expanded Ensemble Single Component Insertion/Deletion */

class testBookkeepingExpanded : public ::testing::Test {
protected:
	//double s, g, lnF, pu;
	std::vector < int > specNmax, specNmin, bounds;
	int Mtot, Nspec;

	virtual void SetUp() {
		Mtot = 3;
		Nspec = 1;
		specNmin.resize(Nspec, 0);
		specNmax.resize(Nspec, 2);
		bounds.resize(2, 0);
		bounds[1] = 2; // window to [0, 2] same as total range [0, 2]
	}
};

TEST_F (testBookkeepingExpanded, partialInsert) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	atom a1;

	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (a1.mState, 0);
	
	// this should increase M state of atom and system
	EXPECT_EQ (mysys.atoms[0][0].mState, 0); // sys.atoms originally just empty
	mysys.insertAtom(0, &a1); // atoms makes a copy
		
	EXPECT_EQ (mysys.numSpecies[0], 0); // N stays the same since not fully inserted yet
	EXPECT_EQ (mysys.getCurrentM(), 1); // M should increase 
	EXPECT_EQ (mysys.atoms[0][0].mState, 1); // atom is copied and its M state increased	
}

TEST_F (testBookkeepingExpanded, sequentialInsert) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	atom a1;

	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (a1.mState, 0);
	
	// this should increase M state of atom and system
	EXPECT_EQ (mysys.atoms[0][0].mState, 0); // sys.atoms originally just empty
	mysys.insertAtom(0, &a1); // atoms makes a copy
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() != &a1);
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 0); // N stays the same since not fully inserted yet
	EXPECT_EQ (mysys.getCurrentM(), 2); // M should increase 
	EXPECT_EQ (mysys.atoms[0][0].mState, 2); // atom M state increased	
		
	// now "fully inserted
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N increased
	EXPECT_EQ (mysys.getCurrentM(), 0); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 0); // atom M state reset	
	
	// insert a second atom (which will hit the upper limit)
	mysys.insertAtom(0, &a1); // atoms makes a copy
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() != &a1);
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][1]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N same
	EXPECT_EQ (mysys.getCurrentM(), 1); // M increased
	EXPECT_EQ (mysys.atoms[0][1].mState, 1); // atom M state increased		
		
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][1]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N same
	EXPECT_EQ (mysys.getCurrentM(), 2); // M increased
	EXPECT_EQ (mysys.atoms[0][1].mState, 2); // atom M state increased	
	
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	
	EXPECT_EQ (mysys.numSpecies[0], 2); // N increased
	EXPECT_EQ (mysys.getCurrentM(), 0); // M reset
	EXPECT_EQ (mysys.atoms[0][1].mState, 0); // atom M state reset	
		
	// hitting upper bound now, can't have a third one
	bool caught = false;
	try {
		mysys.insertAtom(0, &a1);
	} catch (customException &ce) {
		caught = true;
	}
	
	EXPECT_TRUE (caught);
}

TEST_F (testBookkeepingExpanded, sequentialDelete) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	atom a1;

	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (a1.mState, 0);
	
	EXPECT_EQ (mysys.atoms[0][0].mState, 0); // sys.atoms originally just empty
	mysys.insertAtom(0, &a1); // atoms makes a copy
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	// now delete
	mysys.deleteAtom(0, 0);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 0); // N decreased
	EXPECT_EQ (mysys.getCurrentM(), 2); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 2); // atom M state reset	
	
	mysys.deleteAtom(0, 0);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 0); // N same
	EXPECT_EQ (mysys.getCurrentM(), 1); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 1); // atom M state reset	

	mysys.deleteAtom(0, 0);
	
	EXPECT_EQ (mysys.numSpecies[0], 0); // N same
	EXPECT_EQ (mysys.getCurrentM(), 0); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 0); // atom M state reset	
	
	// can't delete anymore
	bool caught = false;
	try {
		mysys.deleteAtom(0, 0);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST_F (testBookkeepingExpanded, deleteOutOfOrder) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	atom a1;

	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (a1.mState, 0);
	EXPECT_EQ (mysys.atoms[0][0].mState, 0);
	
	// insert first atom
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	// insert a second atom
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	
	// now delete the first atom
	mysys.deleteAtom(0, 0);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N decreased
	EXPECT_EQ (mysys.getCurrentM(), 2); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 2); // atom M state reset	
	
	mysys.deleteAtom(0, 0);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N same
	EXPECT_EQ (mysys.getCurrentM(), 1); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 1); // atom M state reset	

	mysys.deleteAtom(0, 0);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N same
	EXPECT_EQ (mysys.getCurrentM(), 0); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 0); // atom M state reset	
}

TEST_F (testBookkeepingExpanded, partialInsertDelete) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	atom a1;

	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (a1.mState, 0);
	EXPECT_EQ (mysys.atoms[0][0].mState, 0);
	
	// insert first atom
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	// insert a second atom
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	
	// now begin deleting the first atom
	mysys.deleteAtom(0, 0);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N decreased
	EXPECT_EQ (mysys.getCurrentM(), 2); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 2); // atom M state reset	
	
	mysys.deleteAtom(0, 0);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N same
	EXPECT_EQ (mysys.getCurrentM(), 1); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 1); // atom M state reset	

	// now insert again
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); 
	EXPECT_EQ (mysys.getCurrentM(), 2);
	EXPECT_EQ (mysys.atoms[0][0].mState, 2); 
	
	// and again
	mysys.insertAtom(0, &mysys.atoms[0][0]);
		
	EXPECT_EQ (mysys.numSpecies[0], 2); 
	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (mysys.atoms[0][0].mState, 0); 
}

TEST_F (testBookkeepingExpanded, badDelete) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	atom a1;

	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (a1.mState, 0);
	EXPECT_EQ (mysys.atoms[0][0].mState, 0);
	
	// insert first atom
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	// insert a second atom
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	
	// now begin deleting the first atom
	mysys.deleteAtom(0, 0);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N decreased
	EXPECT_EQ (mysys.getCurrentM(), 2); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 2); // atom M state reset	
	
	// accidentally try to delete the second atom no the first
	bool caught = false;
	try {
		mysys.deleteAtom(0, 1);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
}

TEST_F (testBookkeepingExpanded, badInsert) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	atom a1;

	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (a1.mState, 0);
	EXPECT_EQ (mysys.atoms[0][0].mState, 0);
	
	// insert first atom
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	// insert a second atom
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	
	// now begin deleting the first atom
	mysys.deleteAtom(0, 0);
	
	// check pointer is correct
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	
	EXPECT_EQ (mysys.numSpecies[0], 1); // N decreased
	EXPECT_EQ (mysys.getCurrentM(), 2); // M should reset
	EXPECT_EQ (mysys.atoms[0][0].mState, 2); // atom M state reset	
	
	// accidentally try to insert any new atoms
	bool caught = false;
	try {
		mysys.insertAtom(0, &a1);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	
	// correctly try to insert the fractional atom
	caught = false;
	try {
		mysys.insertAtom(0, &mysys.atoms[0][0]);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (!caught);
}

class testMulticomponentExpandedMCMove : public ::testing::Test {
protected:
	std::vector < int > specNmax, specNmin, bounds;
	int Mtot, Nspec;

	virtual void SetUp() {
		Mtot = 4;
		Nspec = 2;
		specNmin.resize(Nspec, 0);
		specNmax.resize(Nspec, 2);
		bounds.resize(2, 0);
		bounds[1] = 2; // window to [0, 2] same as total range [0, 2]
	}
};

TEST_F (testMulticomponentExpandedMCMove, selectSpec1) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	hardCore hc11, hc12, hc22;
	std::vector < double > params (2, Mtot);
	hc11.setParameters (params);
	params[0] = 1.0;
	hc12.setParameters (params);
	params[0] = 2.0;
	hc22.setParameters (params);
	mysys.addPotential (0, 0, &hc11, false);
	mysys.addPotential (0, 1, &hc12, false);
	mysys.addPotential (1, 1, &hc22, false);
		
	moves mover;
	
	insertParticle insOne (0, "in1"), insTwo (1, "in2");
	deleteParticle delOne (0, "del1"), delTwo (1, "del2");
	
	mover.addMove (&insOne, 1.0);
	mover.addMove (&insTwo, 1.0);
	mover.addMove (&delOne, 1.0);
	mover.addMove (&delTwo, 1.0);
	
	// manually insert two particles of each into the system
	atom a1;
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(1, &a1); 
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &a1); 
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	
	// now say we choose to set species 2, atom 0 as the partial by deleting it
	mysys.deleteAtom(1, 0);
	mysys.deleteAtom(1, 0); // set to m state = 2 (middle between inserted and deleted)
	
	// have to check this stochastically
	const int iters = 100;
	for (unsigned int i = 0; i < iters; ++i) {
		EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[1][0]);
		EXPECT_TRUE (mysys.getFractionalAtomType() == 1);
		
		mover.makeMove(mysys); // this MUST choose to operate on species 2 
		
		EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[1][0]);
		EXPECT_TRUE (mysys.getFractionalAtomType() == 1);
		
		if (mysys.atoms[1][0].mState == 3) {
			// was inserted, delete to reset
			mysys.deleteAtom(1, 0);
		}
		if (mysys.atoms[1][0].mState == 1) {
			// was deleted, insert a new one to reset
			mysys.insertAtom(1, &mysys.atoms[1][0]);
		}
	}
}


TEST_F (testMulticomponentExpandedMCMove, selectSpec1_moved) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	
	hardCore hc11, hc12, hc22;
	std::vector < double > params (2, Mtot);
	hc11.setParameters (params);
	params[0] = 1.0;
	hc12.setParameters (params);
	params[0] = 2.0;
	hc22.setParameters (params);
	mysys.addPotential (0, 0, &hc11, false);
	mysys.addPotential (0, 1, &hc12, false);
	mysys.addPotential (1, 1, &hc22, false);
		
	moves mover;
	
	insertParticle insOne (0, "in1"), insTwo (1, "in2");
	deleteParticle delOne (0, "del1"), delTwo (1, "del2");
	
	mover.addMove (&insOne, 1.0);
	mover.addMove (&insTwo, 1.0);
	mover.addMove (&delOne, 1.0);
	mover.addMove (&delTwo, 1.0);
	
	// manually insert two particles of each into the system
	atom a1;
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &a1); 
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(1, &a1); 
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &a1); 
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	
	// now say we choose to set species 2, atom 0 as the partial by deleting it
	mysys.deleteAtom(1, 0);
	mysys.deleteAtom(1, 0);
	mysys.deleteAtom(1, 0); // set to m state = 1 this time
	
	// have to check this stochastically
	const int iters = 100;
	bool moved = false;
	atom* mypointer = &mysys.atoms[1][0];
	for (unsigned int i = 0; i < iters; ++i) {
		EXPECT_TRUE (mysys.getFractionalAtom() == mypointer);
		EXPECT_TRUE (mysys.getFractionalAtomType() == 1);
		
		mover.makeMove(mysys); // this MUST choose to operate on species 2 
		
		EXPECT_TRUE (mysys.getFractionalAtom() == mypointer);
		EXPECT_TRUE (mysys.getFractionalAtomType() == 1);
		
		if (mysys.getFractionalAtom()->mState == 2) {
			EXPECT_TRUE (mysys.numSpecies[0] == 2);
			EXPECT_TRUE (mysys.numSpecies[1] == 1);
			EXPECT_TRUE (mysys.getCurrentM() == 2);
			EXPECT_TRUE (mysys.getTotN() == 3);
			
			// was inserted, delete to reset
			if (!moved) {
				mysys.deleteAtom(1, 0);
			} else {
				mysys.deleteAtom(1, 1);
			}
			
			EXPECT_TRUE (mysys.getCurrentM() == 1);
		}
		if (mysys.getFractionalAtom()->mState == 0) {
			EXPECT_TRUE (mysys.numSpecies[0] == 2);
			EXPECT_TRUE (mysys.numSpecies[1] == 1);
			EXPECT_TRUE (mysys.getCurrentM() == 0);
			EXPECT_TRUE (mysys.getTotN() == 3);
			
			// was deleted, insert a new one to reset - this is now at the END of this list with mState = 1
			mysys.insertAtom(1, &a1);

			EXPECT_TRUE (mysys.getCurrentM() == 1);
			
			// thus the partial particle is now (and will always be here)
			mypointer = &mysys.atoms[1][1];
			moved = true;
		}
	}
}

class testComputeBiasExpanded : public ::testing::Test {
protected:
	double s, g, lnF, pu;
	std::vector < int > specNmax, specNmin, bounds;
	int tmmcSweepSize, Mtot, Nspec;

	virtual void SetUp() {
		pu = 0.123;
		s = 0.8;
		g = 0.5;
		lnF = 1.0;
		Mtot = 3;
		tmmcSweepSize = 1;
		Nspec = 2;
		specNmin.resize(Nspec, 0);
		specNmax.resize(Nspec, 2);
		bounds.resize(2, 0);
		bounds[1] = 2; // window to [0, 2] same as total range [0, 2]
	}
};

TEST_F (testComputeBiasExpanded, setCalculateWALABias) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	mysys.setTotNBounds (bounds);
	mysys.startWALA (lnF, 0.5, 0.8, Mtot);
	EXPECT_TRUE (mysys.useWALA);
	
	mysys.getWALABias()->update(bounds[0], 0);
	mysys.getWALABias()->update(bounds[0], 1);
	mysys.getWALABias()->update(bounds[0], 2);
	mysys.getWALABias()->update(bounds[0], 2);
	mysys.getWALABias()->update(bounds[0], 2);
	
	mysys.getWALABias()->update(bounds[0]+1, 0);
	mysys.getWALABias()->update(bounds[0]+1, 0);
	mysys.getWALABias()->update(bounds[0]+1, 1);
	mysys.getWALABias()->update(bounds[0]+1, 2);
	
	mysys.getWALABias()->update(bounds[0]+2, 0);
	mysys.getWALABias()->update(bounds[0]+2, 1);
	mysys.getWALABias()->update(bounds[0]+2, 1);
	mysys.getWALABias()->update(bounds[0]+2, 2);
	mysys.getWALABias()->update(bounds[0]+2, 2);
	mysys.getWALABias()->update(bounds[0]+2, 2);
	
	atom a1;
	
	// N = 0, M = 0 state --> N = 0, M = 1
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0], 1) - exp(-lnF - -lnF)) < 1.0e-9);
	mysys.insertAtom(0, &a1);
	
	// N = 0, M = 1 state --> N = 0, M = 2
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0], 2) - exp(-3*lnF - -lnF)) < 1.0e-9);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
		
	// N = 0, M = 2 state --> N = 1, M = 0
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 0) - exp(-2*lnF - -3*lnF)) < 1.0e-9);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
		
	// now choose to insert atom into a different species

	// N = 1, M = 0 state --> N = 1, M = 1
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 1) - exp(-lnF - -2*lnF)) < 1.0e-9);
	mysys.insertAtom(1, &a1);
	
	// N = 1, M = 1 state --> N = 1, M = 2
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 2) - exp(-lnF - -lnF)) < 1.0e-9);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	
	// N = 1, M = 2 state --> N = 2, M = 0
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 2) - exp(-lnF - -lnF)) < 1.0e-9);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	
	// now delete (completely the first atom), so now in N = Nmin+1, M = 0 state
	mysys.deleteAtom(0, 0);
	mysys.deleteAtom(0, 0);
	mysys.deleteAtom(0, 0);
	
	EXPECT_EQ (mysys.getCurrentM(), 0);
	EXPECT_EQ (mysys.getTotN(), 1);
	
	// propose deletion of other particle --> N = Nmin (0), M = 2
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0], 2) - exp(-3*lnF - -2*lnF)) < 1.0e-9);
	mysys.deleteAtom(1, 0);
	
	// N = 0, M = 2 --> N = 0, M = 1
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0], 1) - exp(-1*lnF - -3*lnF)) < 1.0e-9);
	mysys.deleteAtom(1, 0);
	
	// N = 0, M = 1 --> N = 0, M = 0
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0], 0) - exp(-lnF - -lnF)) < 1.0e-9);
}
	
TEST_F (testComputeBiasExpanded, setCalculateTMMCBias) {
	std::vector < double > ib (3, 10), mu (Nspec, 1.0);
	simSystem mysys (Nspec, 1.0, ib, mu, specNmax, specNmin, Mtot);
	mysys.setTotNBounds (bounds);
	mysys.startTMMC (tmmcSweepSize, Mtot);
	EXPECT_TRUE (mysys.useTMMC);
	
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0], 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0], 0, 1, std::min(1.0, pu));
	
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0], 1, 1, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0], 1, 2, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0], 1, 0, std::min(1.0, pu));
	
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0], 2, 2, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0]+1, 2, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0], bounds[0], 2, 1, std::min(1.0, pu));
	
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 0, 1, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0], 0, 2, std::min(1.0, pu));
	
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 1, 1, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 1, 1, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 1, 1, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 1, 2, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 1, 0, std::min(1.0, pu));
	
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 2, 2, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 2, 1, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+1, 2, 1, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+1, bounds[0]+2, 2, 0, std::min(1.0, pu));
	
	mysys.getTMMCBias()->updateC(bounds[0]+2, bounds[0]+2, 0, 0, std::min(1.0, pu));
	mysys.getTMMCBias()->updateC(bounds[0]+2, bounds[0]+1, 0, 2, std::min(1.0, pu));
		
	mysys.getTMMCBias()->calculatePI();

	// some simple updates 
	
	atom a1;
	
	// Insert atom a1 through stages up to max
	
	// 0 --> 1/3
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0], 1) - exp(-log(3./2.) - 0)) < 1.0e-9); // by construction first lnPI = 0
	mysys.insertAtom(0, &a1);
	
	// 1/3 --> 2/3
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0], 2) - exp(-log(3./2.) - -log(3./2.))) < 1.0e-9);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	// 2/3 --> 1
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 0) - exp(-log(3./2.) - -log(3./2.))) < 1.0e-9);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	
	// 1 --> 4/3
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 1) - exp(-log(5./2.) - -log(3./2.))) < 1.0e-9);
	mysys.insertAtom(0, &a1);
	
	// 4/3 --> 5/3
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 2) - exp(-0 - -log(5./2.))) < 1.0e-9);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	
	// 5/3 --> 2
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+2, 0) - exp(-log(2./4.) - -0)) < 1.0e-9);
	// 5/3 --> 4/3 (backwards check)
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 1) - exp(-log(5./2.) - -0)) < 1.0e-9);

	// remove one of these from the system
	mysys.deleteAtom(0, 1);
	mysys.deleteAtom(0, 1);
	
	// add in a different species to check that bias indepenent of species, just depends on Ntot
	// 1 --> 4/3
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 1) - exp(-log(5./2.) - -log(3./2.))) < 1.0e-9);
	mysys.insertAtom(1, &a1);
	
	// 4/3 --> 5/3
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 2) - exp(-0 - -log(5./2.))) < 1.0e-9);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	
	// 5/3 --> 2
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+2, 0) - exp(-log(2./4.) - -0)) < 1.0e-9);
	// 5/3 --> 4/3 (backwards check)
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+1, 1) - exp(-log(5./2.) - -0)) < 1.0e-9);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	
	// 2 --> 2
	EXPECT_TRUE (fabs(calculateBias (mysys, bounds[0]+2, 0) - exp(0)) < 1.0e-9);
}

TEST_F (testComputeBiasExpanded, testInSituWALASingleComponent) {
	std::vector < double > mu (1, std::numeric_limits<double>::max()); // force an insertion to an empty system
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (1, 3), nmin (1, 0);
	lnF = 3.1; // something random this time
	simSystem mysys (1, 1.0, ib, mu, nmax, nmin, Mtot);
	mysys.startWALA (lnF, 0.5, 0.8, Mtot);
	EXPECT_TRUE (mysys.useWALA);
		
	hardCore hc;
	std::vector < double > params (2, 1.0);
	hc.setParameters (params);
	mysys.addPotential (0, 0, &hc, false);
	
	moves usedMoves;
	insertParticle newIns (0, "insert");
	usedMoves.addMove(&newIns, 1.0);
	
	// will insert
	usedMoves.makeMove(mysys);
	
	// check WALA properties - should have incremented where the system ENDED (at N = 0, M = 1)
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 1)) - -lnF) < 1.0e-9);
	
	// all the rest should be 0
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 0)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 2)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 0)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 1)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 2)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 0)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 1)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 2)) - 0) < 1.0e-9);
	
}

TEST_F (testComputeBiasExpanded, testInSituWALAMultiComponent) {
	std::vector < double > mu (2, std::numeric_limits<double>::max()); // force an insertion to an empty system
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (2, 3), nmin (2, 0);
	lnF = 3.1; // something random this time
	simSystem mysys (2, 1.0, ib, mu, nmax, nmin, Mtot);
	mysys.startWALA (lnF, 0.5, 0.8, Mtot);
	EXPECT_TRUE (mysys.useWALA);
	
	hardCore hc11, hc12, hc22;
	std::vector < double > params (2, Mtot);
	hc11.setParameters (params);
	params[0] = 0.0; // ergo 1 and 2 can sit on top of each other
	hc12.setParameters (params);
	params[1] = 2.0;
	hc22.setParameters (params);
	mysys.addPotential (0, 0, &hc11, false);
	mysys.addPotential (0, 1, &hc12, false);
	mysys.addPotential (1, 1, &hc22, false);
		
	moves usedMoves;
	insertParticle newIns (0, "insert");
	usedMoves.addMove(&newIns, 1.0);
	
	// will insert first species
	usedMoves.makeMove(mysys);
	
	// check WALA properties - should have incremented where the system ends and nowhere else
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 1)) - -lnF) < 1.0e-9);
	
	usedMoves.makeMove(mysys);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 1)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 2)) - -lnF) < 1.0e-9);
	
	usedMoves.makeMove(mysys);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 1)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 2)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 0)) - -lnF) < 1.0e-9);
	
	// all the rest should be 0
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 0)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 1)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 2)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 0)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 1)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 2)) - 0) < 1.0e-9);
	
	// will insert the second species
	moves usedMoves2;
	insertParticle newIns2 (1, "insert");
	usedMoves2.addMove(&newIns2, 1.0);
	
	usedMoves2.makeMove(mysys);
	
	// check WALA properties - should have incremented where the system ends (N = 1, M = 1)
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 1)) - -lnF) < 1.0e-9);
	
	usedMoves2.makeMove(mysys);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 2)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 1)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 1)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 2)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 0)) - -lnF) < 1.0e-9);
	
	usedMoves2.makeMove(mysys);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 0)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 2)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 1)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 1)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(0, 2)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(1, 0)) - -lnF) < 1.0e-9);
	
	// only "junk space" should remain 0
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 1)) - 0) < 1.0e-9);
	EXPECT_TRUE (fabs(mysys.getWALABias()->getBias(mysys.getWALABias()->getAddress(2, 2)) - 0) < 1.0e-9);
}

TEST_F (testComputeBiasExpanded, testInSituTMMCSingleComponent) {
	std::vector < double > mu (1, std::numeric_limits<double>::max()); // force an insertion to an empty system
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (1, 3), nmin (1, 0);
	simSystem mysys (1, 1.0, ib, mu, nmax, nmin, Mtot);
	mysys.startTMMC (tmmcSweepSize, Mtot);
	EXPECT_TRUE (mysys.useTMMC);
		
	hardCore hc;
	std::vector < double > params (2, Mtot);
	hc.setParameters (params);
	mysys.addPotential (0, 0, &hc, false);

	moves usedMoves;
	insertParticle newIns (0, "insert");
	usedMoves.addMove(&newIns, 1.0);
	
	// will insert
	usedMoves.makeMove(mysys);
	
	// check TMMC properties - should have incremented where the system started from (N = 0)
	std::vector < double > C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 2; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// insert again
	usedMoves.makeMove(mysys);
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 6; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// insert again
	usedMoves.makeMove(mysys);
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[6] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[7] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[8] - 0.0) < 1.0e-9);
		
	// all the rest should be 0
	for (unsigned int i = 9; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
}

TEST_F (testComputeBiasExpanded, testInSituTMMCMultiComponent) {
	std::vector < double > mu (2, std::numeric_limits<double>::max()); // force an insertion to an empty system
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (2, 3), nmin (2, 0);
	simSystem mysys (2, 1.0, ib, mu, nmax, nmin, Mtot);
	mysys.startTMMC (tmmcSweepSize, Mtot);
	EXPECT_TRUE (mysys.useTMMC);
		
	hardCore hc11, hc12, hc22;
	std::vector < double > params (2, Mtot);
	hc11.setParameters (params);
	params[0] = 0.0;
	hc12.setParameters (params);
	params[0] = 0.0;
	hc22.setParameters (params);
	mysys.addPotential (0, 0, &hc11, false);
	mysys.addPotential (0, 1, &hc12, false);
	mysys.addPotential (1, 1, &hc22, false);
	
	moves usedMoves;
	insertParticle newIns (0, "insert");
	usedMoves.addMove(&newIns, 1.0);
	
	// will insert first species
	usedMoves.makeMove(mysys);
	
	// check TMMC properties - should have incremented where the system started from (N = 0)
	std::vector < double > C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 2; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// insert same species again
	usedMoves.makeMove(mysys);
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 6; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// finish inserting the same species again
	usedMoves.makeMove(mysys);
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[6] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[7] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[8] - 0.0) < 1.0e-9);
		
	// all the rest should be 0
	for (unsigned int i = 9; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
		
	// do insertions with the other species
	moves usedMoves2;
	insertParticle newIns2 (1, "insert");
	usedMoves2.addMove(&newIns2, 1.0);
	
	// insert 1 atom from second species
	usedMoves2.makeMove(mysys);
	
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[6] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[7] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[8] - 0.0) < 1.0e-9);
	
	// N = 1, M = 1 now
	EXPECT_TRUE (fabs(C[9] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[10] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[11] - 0.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 12; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// second level of insertion
	usedMoves2.makeMove(mysys);
	
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[6] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[7] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[8] - 0.0) < 1.0e-9);
	
	// N = 1, M = 2 now
	EXPECT_TRUE (fabs(C[9] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[10] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[11] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[12] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[13] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[14] - 0.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 15; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
	
	// final level of insertion
	usedMoves2.makeMove(mysys);
	
	C = mysys.getTMMCBias()->getC();
	EXPECT_TRUE (fabs(C[0] - 0.0) < 1.0e-9); // infinite mu, implies p_u = 1, so 1-1 = 0
	EXPECT_TRUE (fabs(C[1] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[2] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[3] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[4] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[5] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[6] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[7] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[8] - 0.0) < 1.0e-9);
	
	// N = 2, M = 0 now
	EXPECT_TRUE (fabs(C[9] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[10] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[11] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[12] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[13] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[14] - 0.0) < 1.0e-9);
	
	EXPECT_TRUE (fabs(C[15] - 0.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[16] - 1.0) < 1.0e-9);
	EXPECT_TRUE (fabs(C[17] - 0.0) < 1.0e-9);
	
	// all the rest should be 0
	for (unsigned int i = 18; i < C.size(); ++i) {
		EXPECT_TRUE (fabs(C[i] - 0.0) < 1.0e-9);
	}
}

/* Check potentials in the expanded ensemble */
class lennardJonesExpandedTest : public ::testing::Test {
protected:
	pairPotential* lj;
	
	int Mtot;
	double eps, sigma, rcut, ushift, tol;
	std::vector < double > params;
	
	virtual void SetUp() {
		Mtot = 3;
		lj = new lennardJones;
		params.resize(5, 0);
		params[4] = Mtot;
		tol = 1.0e-9;
	}
};

TEST_F (lennardJonesExpandedTest, badParams0) {
	eps = -1.0;
	sigma = 1.0;
	rcut = 2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	bool caught = false;
	try {
		lj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete lj;
}

TEST_F (lennardJonesExpandedTest, badParams1) {
	eps = 1.0;
	sigma = -1.0;
	rcut = 2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	bool caught = false;
	try {
		lj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete lj;
}

TEST_F (lennardJonesExpandedTest, badParams2) {
	eps = 1.0;
	sigma = 1.0;
	rcut = -2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	bool caught = false;
	try {
		lj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete lj;
}

TEST_F (lennardJonesExpandedTest, StandardParams_M0) {
	eps = 1.0;
	sigma = 1.0;
	rcut = 2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	lj->setParameters(params);
	atom a1, a2;
	std::vector < double > box(3, 2.1*rcut);
	
	a2.pos[2] = pow(2.0, 1./6.);	
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - -eps) < tol );
	
	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 0.9999*rcut;
	EXPECT_TRUE ( lj->energy(&a1, &a2, box) < 0.0 );
	
	a2.pos[2] = sigma;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 1.234;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - -0.812008901) < tol );
	
	delete lj;
}

TEST_F (lennardJonesExpandedTest, StandardParams_M1) {
	eps = 1.0;
	sigma = 1.0;
	rcut = 2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	lj->setParameters(params);
	atom a1, a2;
	
	a1.mState = 1; // M = 1/3
	
	std::vector < double > box(3, 2.1*rcut);
	
	// Energy should be linearly scaled, sigma gets scaled volumetrically
	a2.pos[2] = pow(2.0, 1./6.)*pow(params[1]*params[1]*params[1]/(8.0*Mtot), 1./3.);	
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - -eps/Mtot) < tol );
	
	// rcut unaffected
	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 0.9999*rcut;
	EXPECT_TRUE ( lj->energy(&a1, &a2, box) < 0.0 );
	
	// at tangency
	a2.pos[2] = pow(params[1]*params[1]*params[1]/(8.0*Mtot), 1./3.);
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	delete lj;
}

TEST_F (lennardJonesExpandedTest, StandardParams_M2) {
	eps = 1.0;
	sigma = 1.0;
	rcut = 2.5;
	ushift = 0.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	lj->setParameters(params);
	atom a1, a2;
	
	a1.mState = 2; // M = 2/3
	
	std::vector < double > box(3, 2.1*rcut);
	
	// Energy should be linearly scaled, sigma gets scaled volumetrically
	double r1_cub = params[1]*params[1]*params[1]/(8.0*Mtot);
	a2.pos[2] = pow(2.0, 1./6.)*pow(params[1]*params[1]*params[1]/(8.0*Mtot)+r1_cub, 1./3.);	
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - -2.0*eps/Mtot) < tol );
	
	// rcut unaffected
	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 0.9999*rcut;
	EXPECT_TRUE ( lj->energy(&a1, &a2, box) < 0.0 );
	
	// at tangency
	a2.pos[2] = pow(params[1]*params[1]*params[1]/(8.0*Mtot)+r1_cub, 1./3.);
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	delete lj;
}

TEST_F (lennardJonesExpandedTest, WCAparams_M0) {
	eps = 1.0;
	sigma = 1.0;
	rcut = pow(2.0, 1./6.);
	ushift = 1.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	lj->setParameters(params);
	atom a1, a2;
	std::vector < double > box(3, 2.1*rcut);
	
	a1.mState = 0;
	
	a2.pos[2] = pow(2.0, 1./6.);
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0.0) < tol );
	
	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 0.9999*rcut;
	EXPECT_TRUE ( lj->energy(&a1, &a2, box) > 0.0 );
	
	a2.pos[2] = 0.987;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 1.353386603) < tol );
	
	delete lj;
}

TEST_F (lennardJonesExpandedTest, WCAparams_M1) {
	eps = 1.0;
	sigma = 1.0;
	rcut = pow(2.0, 1./6.);
	ushift = 1.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	lj->setParameters(params);
	atom a1, a2;
	std::vector < double > box(3, 2.1*rcut);
	
	a1.mState = 1;
	
	a2.pos[2] = pow(2.0, 1./6.)*pow(params[1]*params[1]*params[1]/(8.0*Mtot), 1./3.);
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0.0) < tol );
	
	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 0.9999*rcut;
	EXPECT_TRUE ( lj->energy(&a1, &a2, box) > 0.0 );
	
	a2.pos[2] = 0.987*pow(params[1]*params[1]*params[1]/(8.0*Mtot), 1./3.);
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 1.353386603/Mtot) < tol );
	
	delete lj;
}

TEST_F (lennardJonesExpandedTest, WCAparams_M2) {
	eps = 1.0;
	sigma = 1.0;
	rcut = pow(2.0, 1./6.);
	ushift = 1.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	params[3] = ushift;
	lj->setParameters(params);
	atom a1, a2;
	std::vector < double > box(3, 2.1*rcut);
	
	a1.mState = 2;
	
	double r1_cub = params[1]*params[1]*params[1]/(8.0*Mtot);
	a2.pos[2] = pow(2.0, 1./6.)*pow(params[1]*params[1]*params[1]/(8.0*Mtot) + r1_cub, 1./3.);
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0.0) < tol );
	
	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 0) < tol );
	
	a2.pos[2] = 0.9999*rcut;
	EXPECT_TRUE ( lj->energy(&a1, &a2, box) > 0.0 );
	
	a2.pos[2] = 0.987*pow(params[1]*params[1]*params[1]/(8.0*Mtot) + r1_cub, 1./3.);
	EXPECT_TRUE ( fabs(lj->energy(&a1, &a2, box) - 2.0*1.353386603/Mtot) < tol );
	
	delete lj;
}

class squareWellExpandedTest : public ::testing::Test {
protected:
	squareWell *sw;
	double eps, sigma, width, tol;
	std::vector < double > params;
	int Mtot;
	
	virtual void SetUp() {
		sw = new squareWell;
		eps = 1.234;
		sigma = 2.345;
		width = 0.12;
		Mtot = 3;
		params.resize(4, 0);
		params[0] = sigma;
		params[1] = width;
		params[2] = eps;
		params[3] = Mtot;
		tol = 1.0e-9;
	}
};

TEST_F (squareWellExpandedTest, badParams0) {
	params[0] = -1.0;
	bool caught = false;
	try {
		sw->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete sw;
}

TEST_F (squareWellExpandedTest, badParams1) {
	params[1] = -1.0;
	bool caught = false;
	try {
		sw->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete sw;
}

TEST_F (squareWellExpandedTest, badParams2) {
	params[2] = -1.0;
	bool caught = false;
	try {
		sw->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete sw;
}

TEST_F (squareWellExpandedTest, testInRangeLower_M0) {	
	sw->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a1.mState = 0;
	
	a2.pos[2] = 1.0001*sigma;

	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - -eps) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testInRangeLower_M1) {	
	sw->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a1.mState = 1;
	
	a2.pos[2] = 1.0001*pow(params[0]*params[0]*params[0]/(8.0*Mtot), 1./3.);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - -eps/Mtot) < 1.0e-9);

	delete sw;
}

TEST_F (squareWellExpandedTest, testInRangeLower_M2) {	
	sw->setParameters(params);
	atom a1, a2;
	std::vector < double > box (3, 10);
	
	a1.mState = 2;
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 1.0001*pow(params[0]*params[0]*params[0]/(8.0*Mtot)+r1_cub, 1./3.); 

	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - -2.0*eps/Mtot) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testInRangeUpper_M0) {
	sw->setParameters(params);
	atom a1, a2;
	a1.mState = 0;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 0.9999*(sigma+width);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - -eps) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testInRangeUpper_M1) {
	sw->setParameters(params);
	atom a1, a2;
	a1.mState = 1;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 0.9999*(pow(params[0]*params[0]*params[0]/(8.0*Mtot), 1./3.)+width);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - -eps/Mtot) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testInRangeUpper_M2) {
	sw->setParameters(params);
	atom a1, a2;
	a1.mState = 2;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 0.9999*(pow(params[0]*params[0]*params[0]/(8.0*Mtot)+r1_cub, 1./3.)+width);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - -2.0*eps/Mtot) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testOutOfRange_M0) {
	sw->setParameters(params);
	atom a1, a2;
	a2.mState = 0;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 1.0001*(sigma+width);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - 0) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testOutOfRange_M1) {
	sw->setParameters(params);
	atom a1, a2;
	a2.mState = 1;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 1.0001*(pow(r1_cub, 1./3.)+width);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - 0) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testOutOfRange_M2) {
	sw->setParameters(params);
	atom a1, a2;
	a2.mState = 2;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 1.0001*(pow(params[0]*params[0]*params[0]/(8.0*Mtot)+r1_cub, 1./3.)+width);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - 0) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testBadRange_M0) {
	sw->setParameters(params);
	atom a1, a2;
	a1.mState = 0;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 0.9999*(sigma);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - NUM_INFINITY) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testBadRange_M1) {
	sw->setParameters(params);
	atom a1, a2;
	a1.mState = 1;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 0.9999*pow(r1_cub, 1./3.);
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - NUM_INFINITY) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testBadRange_M2) {
	sw->setParameters(params);
	atom a1, a2;
	a1.mState = 2;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 0.9999*(pow(params[0]*params[0]*params[0]/(8.0*Mtot)+r1_cub, 1./3.));
	EXPECT_TRUE (fabs(sw->energy(&a1, &a2, box) - NUM_INFINITY) < tol);
	delete sw;
}

TEST_F (squareWellExpandedTest, testRcut) {
	sw->setParameters(params);
	EXPECT_TRUE (fabs(sw->rcut() - (sigma+width)) < tol); // rcut unaffected by Mtot > 0
	delete sw;
}

TEST_F (squareWellExpandedTest, testTailCorrection) {
	sw->setParameters(params);
	EXPECT_TRUE (fabs(sw->tailCorrection(0.1234) - 0) < tol); // I have chosen not to scale this
	delete sw;
}

class hardCoreExpandedTest : public ::testing::Test {
protected:
	hardCore *hc;
	double sigma, tol;
	std::vector < double > params;
	int Mtot;
	
	virtual void SetUp() {
		hc = new hardCore;
		Mtot = 3;
		sigma = 2.345;
		params.resize(2, 0);
		params[0] = sigma;
		params[1] = Mtot;
		tol = 1.0e-9;
	}
};

TEST_F (hardCoreExpandedTest, badParams0) {
	params[0] = -1.0;
	bool caught = false;
	try {
		hc->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete hc;
}

TEST_F (hardCoreExpandedTest, testRange_M0) {
	hc->setParameters(params);
	atom a1, a2;
	a1.mState = 0;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 1.0001*sigma;
	EXPECT_TRUE (fabs(hc->energy(&a1, &a2, box) - 0) < tol);
	delete hc;
}

TEST_F (hardCoreExpandedTest, testRange_M1) {
	hc->setParameters(params);
	atom a1, a2;
	a1.mState = 1;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 1.0001*pow(r1_cub, 1./3.);
	EXPECT_TRUE (fabs(hc->energy(&a1, &a2, box) - 0) < tol);
	delete hc;
}

TEST_F (hardCoreExpandedTest, testRange_M2) {
	hc->setParameters(params);
	atom a1, a2;
	a1.mState = 2;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 1.0001*pow( params[0]*params[0]*params[0]/(8.0*Mtot)+r1_cub, 1./3.);
	EXPECT_TRUE (fabs(hc->energy(&a1, &a2, box) - 0) < tol);
	delete hc;
}

TEST_F (hardCoreExpandedTest, testBadRange_M0) {
	hc->setParameters(params);
	atom a1, a2;
	a1.mState = 0;
	std::vector < double > box (3, 10);
	
	a2.pos[2] = 0.9999*sigma;
	EXPECT_TRUE (fabs(hc->energy(&a1, &a2, box) - NUM_INFINITY) < tol);
	delete hc;
}

TEST_F (hardCoreExpandedTest, testBadRange_M1) {
	hc->setParameters(params);
	atom a1, a2;
	a1.mState = 1;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 0.9999*pow(r1_cub, 1./3.);
	EXPECT_TRUE (fabs(hc->energy(&a1, &a2, box) - NUM_INFINITY) < tol);
	delete hc;
}

TEST_F (hardCoreExpandedTest, testBadRange_M2) {
	hc->setParameters(params);
	atom a1, a2;
	a1.mState = 2;
	std::vector < double > box (3, 10);
	
	double r1_cub = params[0]*params[0]*params[0]/(8.0*Mtot);
	a2.pos[2] = 0.9999*pow(params[0]*params[0]*params[0]/(8.0*Mtot)+r1_cub, 1./3.);
	EXPECT_TRUE (fabs(hc->energy(&a1, &a2, box) - NUM_INFINITY) < tol);
	delete hc;
}

TEST_F (hardCoreExpandedTest, testRcut) {
	hc->setParameters(params);
	EXPECT_TRUE (fabs(hc->rcut() - sigma) < tol); // unaffected by Mtot > 0
	delete hc;
}

TEST_F (hardCoreExpandedTest, testTailCorrection) {
	hc->setParameters(params);
	EXPECT_TRUE (fabs(hc->tailCorrection(0.1234) - 0) < tol); // unscaled with Mtot
	delete hc;
}

// add expanded ensemble, 2 component SWAP move testing