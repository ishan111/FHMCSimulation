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
		tmmc test (Nmax, Nmin, sweepSize, 1, dummyBox);
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
		tmmc test (Nmax, Nmin, sweepSize, 1, dummyBox);
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
		tmmc test (Nmax, Nmin, sweepSize, 1, dummyBox);
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
		std::cout << ce.what() << std::endl;
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

// test swap move on very specific cases, also consider doing this for inserts and deletes too with/without cell lists

// review and optimize use of __BIAS_TYPE_INT__
