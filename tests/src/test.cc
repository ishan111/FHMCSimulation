#include "gtest-1.7.0/include/gtest/gtest.h"
#include <vector>
#include <iostream>
#include <limits>
#include "../../src/fhmc.h"

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
	test.addInsert(0, prob1);
	std::vector < double > pr = test.reportProbabilities();
	EXPECT_TRUE (fabs(pr[0] - 1.0) < tol);
	test.addInsert(0, prob2);
	std::vector < double > pr2 = test.reportProbabilities();
	EXPECT_TRUE (fabs(pr2[0] - prob1) < tol);
	EXPECT_TRUE (fabs(pr2[1] - (prob1+prob2)) < tol);
}

/* periodic boundary distances */
TEST (PBC, TwoVectors) {
	double tol = 1.0e-9, L = 10;
	std::vector < double > p1(3, 0), p2 (3, 1), p3 (3, L/2.0), p4 (3, 0.75*L), box(3, L);
	EXPECT_TRUE (fabs(pbcDist2(p1, p2, box) - 3.0) < tol );
	EXPECT_TRUE (fabs(pbcDist2(p1, p3, box) - (3.0*L*L/4.0)) < tol );
	EXPECT_TRUE (fabs(pbcDist2(p1, p4, box) - (0.25*L*0.25*L*3.0)) < tol );
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

class fsLennardJonesTest : public ::testing::Test {
protected:
	pairPotential* fslj;

	double eps, sigma, rcut, tol;
	std::vector < double > params;

	virtual void SetUp() {
		fslj = new fsLennardJones;
		params.resize(4, 0);
		params[3] = 1;
		tol = 1.0e-9;
	}
};

TEST_F (fsLennardJonesTest, badParams0) {
	eps = -1.0;
	sigma = 1.0;
	rcut = 3.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	bool caught = false;
	try {
		fslj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete fslj;
}

TEST_F (fsLennardJonesTest, badParams1) {
	eps = 1.0;
	sigma = -1.0;
	rcut = 3.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	bool caught = false;
	try {
		fslj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete fslj;
}

TEST_F (fsLennardJonesTest, badParams2) {
	eps = 1.0;
	sigma = 1.0;
	rcut = -2.5;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	bool caught = false;
	try {
		fslj->setParameters(params);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete fslj;
}

TEST_F (fsLennardJonesTest, StandardParams) {
	eps = 1.0;
	sigma = 1.0;
	rcut = 3.0;
	params[0] = eps;
	params[1] = sigma;
	params[2] = rcut;
	fslj->setParameters(params);
	atom a1, a2;
	std::vector < double > box(3, 2.1*rcut);

	a2.pos[2] = pow(2.0, 1./6.);
	EXPECT_TRUE ( fabs(fslj->energy(&a1, &a2, box) - -0.973973102) < tol );

	a2.pos[2] = 1.0001*rcut;
	EXPECT_TRUE ( fabs(fslj->energy(&a1, &a2, box) - 0) < tol );

	a2.pos[2] = 2.99;
	EXPECT_TRUE ( fabs(fslj->energy(&a1, &a2, box) - -0.0000012851) < tol );

	a2.pos[2] = sigma;
	EXPECT_TRUE ( fabs(fslj->energy(&a1, &a2, box) - 0.0273671019) < tol );

	a2.pos[2] = 1.4;
	EXPECT_TRUE ( fabs(fslj->energy(&a1, &a2, box) - -0.4376973523) < tol );

	a2.pos[2] = 0.9;
	EXPECT_TRUE ( fabs(fslj->energy(&a1, &a2, box) - 6.6645804382) < tol );

	delete fslj;
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
	std::vector < double > params(5, 1);
	mysys.addPotential(0, 0, "lennard_jones", params, true);
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
	std::vector < double > params (5), p1(3, 0.0), p2(3, 0.0);

	// wca
	params[0] = 1.0;
	params[1] = 1.0;
	params[2] = pow(2.0, 1./6.);
	params[3] = 1.0;
	params[4] = 1; // Mtot
	mysys.addPotential(0, 0, "lennard_jones", params, true);

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
	mysys.ppot[0][0]->setParameters(params);

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

TEST (testTMMC, checkReadHCGetHC) {
	std::vector < double > dummyBox (3, 0);
	tmmc* tmmcBias2 = new tmmc (20, 0, 1, 100, dummyBox);
	tmmcBias2->readHC("../data/tmmc_HC.dat");
	const std::vector < double > hc = tmmcBias2->getHC();
	EXPECT_EQ (hc.size(), 3*21);
	for (unsigned int i = 0; i < 5; ++i) {
		EXPECT_EQ ((1+i)*10, hc[i]);
	}
	for (unsigned int i = 5; i < hc.size(); ++i) {
		EXPECT_EQ (hc[i], 0);
	}
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
	tmmcBias->readC("tmmBiaslnPI_checkPrint_C.dat");
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

TEST (testWala, checkReadHGetH) {
	std::vector < double > dummyBox (3, 0);
	wala walaBias (1.0, 0.5, 0.8, 20, 0, 1, dummyBox);
	walaBias.readH("../data/wala_H.dat");
	const std::vector < double > h = walaBias.getH();
	EXPECT_EQ (h.size(), 21);
	for (unsigned int i = 0; i < 5; ++i) {
		EXPECT_EQ (h[i], (1+i)*10);
	}
	for (unsigned int i = 5; i < h.size(); ++i) {
		EXPECT_EQ (h[i], 0);
	}
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
	delete walaBias;
}

TEST_F (testWalaBias, testUpdateSame) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin, 0);
	std::vector < double > H = walaBias->getH(), lnPI = walaBias->getlnPI();
	EXPECT_TRUE (fabs(H[0] - 3.0) < 1.0e-9);
	EXPECT_TRUE (fabs(lnPI[0] - 3.0*lnF) < 1.0e-9);
	delete walaBias;
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
	delete walaBias;
}

TEST_F (testWalaBias, GetGoodAddress) {
	EXPECT_EQ (walaBias->getAddress(Nmin, 0), 0);
	EXPECT_EQ (walaBias->getAddress(Nmin+1, 0), 1);
	EXPECT_EQ (walaBias->getAddress(Nmin+2, 0), 2);
	delete walaBias;
}

TEST_F (testWalaBias, GetBadAddress) {
	bool caught = false;
	try {
		EXPECT_EQ (walaBias->getAddress(Nmin+3, 0), 0);
	} catch (customException &ce) {
		caught = true;
	}
	EXPECT_TRUE (caught);
	delete walaBias;
}

TEST_F (testWalaBias, getBias) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin, 0)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin+1, 0)) - -lnF) < 1.0e-9);
	EXPECT_TRUE (fabs(walaBias->getBias(walaBias->getAddress(Nmin+2, 0)) - -lnF) < 1.0e-9);
	delete walaBias;
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
	delete walaBias;
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

	walaBias->readlnPI("walaBiaslnPI_checkPrint_lnPI.dat");

	std::vector < double > lnPI_new = walaBias->getlnPI();
	for (unsigned int i = 0; i < lnPI_new.size(); ++i) {
		EXPECT_TRUE (fabs(lnPI_new[i] - lnPI_ref[i]) < 1.0e-6); // read loses precision
	}
	delete walaBias;
}

TEST_F (testWalaBias, checkEvaluateFlatnessNo) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	// average of [1, 2, 1] = 1.333, average * s (=0.8) = 1.06 > min (=1)
	bool flat = walaBias->evaluateFlatness();
	EXPECT_TRUE (!flat);
	delete walaBias;
}

TEST_F (testWalaBias, checkEvaluateFlatnessYes) {
	walaBias->update(Nmin, 0);
	walaBias->update(Nmin+1, 0);
	walaBias->update(Nmin+2, 0);
	// average of [1, 1, 1] = 1, average * s (=0.8) = 0.8 < min (=1)
	bool flat = walaBias->evaluateFlatness();
	EXPECT_TRUE (flat);
	delete walaBias;
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

	std::vector < double > params (2, 1.0);
	mysys.addPotential (0, 0, "hard_sphere", params, false);

	moves usedMoves;
	usedMoves.addInsert(0, 1.0);

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

	//hardCore hc11, hc12, hc22;
	std::vector < double > params (2, 1.0);
	//hc11.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);
	params[0] = 0.0; // ergo 1 and 2 can sit on top of each other
	//hc12.setParameters (params);
	mysys.addPotential (0, 1, "hard_sphere", params, false);
	params[1] = 2.0;
	//hc22.setParameters (params);
	mysys.addPotential (1, 1, "hard_sphere", params, false);

	moves usedMoves;
	usedMoves.addInsert(0, 1.0);

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
	usedMoves2.addInsert(1, 1.0);
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

	//hardCore hc;
	std::vector < double > params (2, 1.0);
	//hc.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);

	moves usedMoves;
	usedMoves.addInsert(0, 1.0);

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

	//hardCore hc11, hc12, hc22;
	std::vector < double > params (2, 1.0);
	//hc11.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);
	params[0] = 0.0;
	//hc12.setParameters (params);
	mysys.addPotential (0, 1, "hard_sphere", params, false);
	params[0] = 0.0;
	//hc22.setParameters (params);
	mysys.addPotential (1, 1, "hard_sphere", params, false);

	moves usedMoves;
	usedMoves.addInsert(0, 1.0);

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
	usedMoves2.addInsert(1, 1.0);

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

	//squareWell sw12, sw13, sw23; // only bother to set cross interactions here
	std::vector < double > params (4, 0.0);
	params[0] = 1.0; params[1] = 0.5; params[2] = 1.0, params[3] = 1;
	//sw12.setParameters (params);
	mysys.addPotential (0, 1, "square_well", params, true);
	params[0] = 1.0; params[1] = 0.5; params[2] = 0.3, params[3] = 1;
	//sw13.setParameters (params);
	mysys.addPotential (0, 2, "square_well", params, true);
	params[0] = 1.0; params[1] = 0.5; params[2] = 2.0, params[3] = 1;
	//sw23.setParameters (params);
	mysys.addPotential (1, 2, "square_well", params, true);

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
	usedMoves.addSwap(0,1,1.0); // only can swap a1 and a2

	// swap a1 & a2 should result in no change since they are an isolated pair
	bool noMove = true;
	double lastAns = 0, U_save = mysys.scratchEnergy();
	while (noMove) {
		usedMoves.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves.reportMoveStatistics();
		if (ans[0][0] > lastAns) {
			noMove = false;
		}
		lastAns = ans[0][0];
	}
	double dU1 = mysys.scratchEnergy();
	EXPECT_TRUE (fabs(U_save - dU1) < 1.0e-9);

	// now try to swap a1 and a3 after swapping a1 and a2
	moves usedMoves2;
	usedMoves2.addSwap(1,2,1.0); // swap a2 and a3
	noMove = true;
	lastAns = 0;
	U_save = mysys.scratchEnergy();
	while (noMove) {
		usedMoves2.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves2.reportMoveStatistics();
		if (ans[0][0] > lastAns) {
			noMove = false;
		}
		lastAns = ans[0][0];
	}
	EXPECT_TRUE (fabs(mysys.scratchEnergy() - -0.3) < 1.0e-9); // only 1-3 interaction remains

	// swap a1 and a3 which are now an isolated pair
	moves usedMoves3;
	usedMoves3.addSwap(0,2,1.0); // only can swap a1 and a3
	noMove = true;
	lastAns = 0;
	U_save = mysys.scratchEnergy();
	while (noMove) {
		usedMoves3.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves3.reportMoveStatistics();
		if (ans[0][0] > lastAns) {
			noMove = false;
		}
		lastAns = ans[0][0];
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

	tmmcBias->readC("tmmBiasExpandedlnPI_checkPrint_C.dat");

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
	EXPECT_TRUE (fabs(H[1] - 3.0) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[1] - 3.0*lnF) < 1.0e-8);
}

TEST_F (testExpandedWalaBias, testUpdateDiff) {
	walaBias->update(Nmin, 1);
	walaBias->update(Nmin+1, 2);
	walaBias->update(Nmin+2, 0);
	std::vector < double > H = walaBias->getH(), lnPI = walaBias->getlnPI();
	EXPECT_TRUE (fabs(H[1] - 1.0) < 1.0e-8);
	EXPECT_TRUE (fabs(H[5] - 1.0) < 1.0e-8);
	EXPECT_TRUE (fabs(H[6] - 1.0) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[0] - 0.0) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[1] - lnF) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[2] - 0.0) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[3] - 0.0) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[4] - 0.0) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[5] - lnF) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[6] - lnF) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[7] - 0.0) < 1.0e-8);
	EXPECT_TRUE (fabs(lnPI[8] - 0.0) < 1.0e-8);
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

	walaBias->readlnPI("walaBiasExpandedlnPI_checkPrint_lnPI.dat");

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

	//hardCore hc11, hc12, hc22;
	std::vector < double > params (2, Mtot);
	//hc11.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);
	params[0] = 1.0;
	//hc12.setParameters (params);
	mysys.addPotential (0, 1, "hard_sphere", params, false);
	params[0] = 2.0;
	//hc22.setParameters (params);
	mysys.addPotential (1, 1, "hard_sphere", params, false);

	moves mover (Mtot);
	mover.addInsert(0,1.0);
	mover.addInsert(1,1.0);
	mover.addDelete(0,1.0);
	mover.addDelete(1,1.0);

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

	//hardCore hc11, hc12, hc22;
	std::vector < double > params (2, Mtot);
	//hc11.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);
	params[0] = 1.0;
	//hc12.setParameters (params);
	mysys.addPotential (0, 1, "hard_sphere", params, false);
	params[0] = 2.0;
	//hc22.setParameters (params);
	mysys.addPotential (1, 1, "hard_sphere", params, false);

	moves mover (Mtot);
	mover.addInsert(0,1.0);
	mover.addInsert(1,1.0);
	mover.addDelete(0,1.0);
	mover.addDelete(1,1.0);

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

	//hardCore hc;
	std::vector < double > params (2, 1.0);
	//hc.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);

	moves usedMoves (Mtot);
	usedMoves.addInsert(0, 1.0);

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

	//hardCore hc11, hc12, hc22;
	std::vector < double > params (2, Mtot);
	//hc11.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);
	params[0] = 0.0; // ergo 1 and 2 can sit on top of each other
	//hc12.setParameters (params);
	mysys.addPotential (0, 1, "hard_sphere", params, false);
	params[1] = 2.0;
	//hc22.setParameters (params);
	mysys.addPotential (1, 1, "hard_sphere", params, false);

	moves usedMoves (Mtot);
	usedMoves.addInsert(0,1.0);

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
	moves usedMoves2 (Mtot);
	usedMoves2.addInsert(1,1.0);

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

	//hardCore hc;
	std::vector < double > params (2, Mtot);
	//hc.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);

	moves usedMoves (Mtot);
	usedMoves.addInsert(0, 1.0);

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

	//hardCore hc11, hc12, hc22;
	std::vector < double > params (2, Mtot);
	//hc11.setParameters (params);
	mysys.addPotential (0, 0, "hard_sphere", params, false);
	params[0] = 0.0;
	//hc12.setParameters (params);
	mysys.addPotential (0, 1, "hard_sphere", params, false);
	params[0] = 0.0;
	//hc22.setParameters (params);
	mysys.addPotential (1, 1, "hard_sphere", params, false);

	moves usedMoves (Mtot);
	usedMoves.addInsert(0, 1.0);

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
	moves usedMoves2 (Mtot);
	usedMoves2.addInsert(1,1.0);

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

TEST (testExpandedSwapMove, multicomponentNoSwapTwoFullyInserted) {
	const int Mtot = 3;
	std::vector < double > mu (2, 0);
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (2, 3), nmin (2, 0);
	simSystem mysys (2, 1.0, ib, mu, nmax, nmin, Mtot);

	//squareWell sw11, sw12, sw22;
	std::vector < double > params (4, 0.0);
	params[0] = 1.0; params[1] = 0.1; params[2] = 1.0, params[3] = Mtot;
	//sw11.setParameters (params);
	mysys.addPotential (0, 0, "square_well", params, true);
	params[0] = 1.5; params[1] = 0.1; params[2] = 1.5, params[3] = Mtot;
	//sw12.setParameters (params);
	mysys.addPotential (0, 1, "square_well", params, true);
	params[0] = 2.0; params[1] = 0.1; params[2] = 2.0, params[3] = Mtot;
	//sw22.setParameters (params);
	mysys.addPotential (1, 1, "square_well", params, true);

	atom a1, a2, a3, a4;

	a1.pos[0] = 7;
	a1.pos[1] = 0;
	a1.pos[2] = 5;

	a2.pos[0] = 8.01;
	a2.pos[1] = 0;
	a2.pos[2] = 5;

	a3.pos[0] = 0;
	a3.pos[1] = 0;
	a3.pos[2] = 5;

	a4.pos[0] = 2.01;
	a4.pos[1] = 0;
	a4.pos[2] = 5;

	// Fully insert 4 atoms
	mysys.insertAtom(0, &a1);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);

	mysys.insertAtom(0, &a2);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);

	mysys.insertAtom(1, &a3);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &mysys.atoms[1][0]);

	mysys.insertAtom(1, &a4);
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	mysys.insertAtom(1, &mysys.atoms[1][1]);

	double U_save = mysys.scratchEnergy();
	EXPECT_TRUE (fabs(U_save - -(1.0+2.0)) < 1.0e-9);

	moves usedMoves (Mtot);
	usedMoves.addSwap(0,1,1.0); // only can swap a1 and a2

	bool noMove = true;
	double lastAns = 0;
	int iterMax = 1000, iter = 0;
	while (noMove && iter < iterMax) {
		usedMoves.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves.reportMoveStatistics();
		if (ans[0][0] > lastAns) {
			noMove = false;
		}
		lastAns = ans[0][0];
		iter++;
	}

	// this should nevery succeed since two large and two small particles form different clusters
	EXPECT_TRUE(noMove);
}

TEST (testExpandedSwapMove, multicomponentAllowSwapTwoFullyInserted) {
	const int Mtot = 3;
	std::vector < double > mu (2, 0);
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (2, 3), nmin (2, 0);
	simSystem mysys (2, 1.0, ib, mu, nmax, nmin, Mtot);

	//squareWell sw11, sw12, sw22;
	std::vector < double > params (4, 0.0);
	params[0] = 1.0; params[1] = 0.1; params[2] = 1.0, params[3] = Mtot;
	//sw11.setParameters (params);
	mysys.addPotential (0, 0, "square_well", params, true);
	params[0] = 1.5; params[1] = 0.1; params[2] = 1.5, params[3] = Mtot;
	//sw12.setParameters (params);
	mysys.addPotential (0, 1, "square_well", params, true);
	params[0] = 2.0; params[1] = 0.1; params[2] = 2.0, params[3] = Mtot;
	//sw22.setParameters (params);
	mysys.addPotential (1, 1, "square_well", params, true);

	atom a1, a2, a3, a4;

	a1.pos[0] = 6.5;
	a1.pos[1] = 0;
	a1.pos[2] = 5;

	a2.pos[0] = 8.01;
	a2.pos[1] = 0;
	a2.pos[2] = 5;

	a3.pos[0] = 0;
	a3.pos[1] = 0;
	a3.pos[2] = 5;

	a4.pos[0] = 2.01;
	a4.pos[1] = 0;
	a4.pos[2] = 5;

	// Fully insert 4 atoms
	mysys.insertAtom(0, &a1);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);

	mysys.insertAtom(0, &a2);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);

	mysys.insertAtom(1, &a3);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &mysys.atoms[1][0]);

	mysys.insertAtom(1, &a4);
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	mysys.insertAtom(1, &mysys.atoms[1][1]);

	double U_save = mysys.scratchEnergy();
	EXPECT_TRUE (fabs(U_save - -(0.0+2.0)) < 1.0e-9);

	moves usedMoves (Mtot);
	usedMoves.addSwap(0,1,1.0); // only can swap a1 and a2

	bool noMove = true;
	double lastAns = 0;
	int iterMax = 1000, iter = 0;
	while (noMove && iter < iterMax) {
		usedMoves.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves.reportMoveStatistics();
		if (ans[0][0] > lastAns) {
			noMove = false;
		}
		lastAns = ans[0][0];
		iter++;
	}

	// this should succeed at somepoint since 0 < dU < infinity
	EXPECT_TRUE(!noMove);
}

TEST (testExpandedSwapMove, multicomponentAllowSingleSwapTwoFullyInserted) {
	const int Mtot = 3;
	std::vector < double > mu (2, 0);
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (2, 3), nmin (2, 0);
	simSystem mysys (2, 1.0, ib, mu, nmax, nmin, Mtot);

	//squareWell sw11, sw12, sw22;
	std::vector < double > params (4, 0.0);
	params[0] = 1.0; params[1] = 0.1; params[2] = 1.0, params[3] = Mtot;
	//sw11.setParameters (params);
	mysys.addPotential (0, 0, "square_well", params, true);
	params[0] = 1.5; params[1] = 0.1; params[2] = NUM_INFINITY, params[3] = Mtot; // once 1-2 come together they won't separate
	//sw12.setParameters (params);
	mysys.addPotential (0, 1, "square_well", params, true);
	params[0] = 2.0; params[1] = 0.1; params[2] = 2.0, params[3] = Mtot;
	//sw22.setParameters (params);
	mysys.addPotential (1, 1, "square_well", params, true);

	atom a1, a2, a3, a4;

	a1.pos[0] = 6.5;
	a1.pos[1] = 0;
	a1.pos[2] = 5;

	a2.pos[0] = 8.01;
	a2.pos[1] = 0;
	a2.pos[2] = 5;

	a3.pos[0] = 0;
	a3.pos[1] = 0;
	a3.pos[2] = 5;

	a4.pos[0] = 2.01;
	a4.pos[1] = 0;
	a4.pos[2] = 5;

	// Fully insert 4 atoms
	mysys.insertAtom(0, &a1);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);

	mysys.insertAtom(0, &a2);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);

	mysys.insertAtom(1, &a3);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &mysys.atoms[1][0]);

	mysys.insertAtom(1, &a4);
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	mysys.insertAtom(1, &mysys.atoms[1][1]);

	double U_save = mysys.scratchEnergy();
	EXPECT_TRUE (fabs(U_save - -(0.0+2.0)) < 1.0e-9);

	moves usedMoves (Mtot);
	usedMoves.addSwap(0,1,1.0); // only can swap a1 and a2

	bool done = false, badSwap = true;
	double lastAns = 0;
	int iterMax = 1000, iter = 0;
	while (iter < iterMax && !done) {
		usedMoves.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves.reportMoveStatistics();

		if (ans[0][0] > lastAns) {
			if ((mysys.atoms[0][0].pos[0] == 8.01 || mysys.atoms[0][1].pos[0] == 8.01) && (mysys.atoms[1][0].pos[0] == 6.5 || mysys.atoms[1][1].pos[0] == 6.5)) {
				badSwap = false;
			} else if ((mysys.atoms[0][0].pos[0] == 6.5 || mysys.atoms[0][1].pos[0] == 6.5) && (mysys.atoms[1][0].pos[0] == 8.01 || mysys.atoms[1][1].pos[0] == 8.01)) {
				badSwap = false;
			} else {
				badSwap = true;
			}
			done = true;
		}
		lastAns = ans[0][0];
		iter++;
	}

	// this should only succeed once since breaking like clusters has finite energy, but after that swap, have 1-2 pairs which have infinite attraction
	EXPECT_TRUE (done);
	EXPECT_TRUE (!badSwap);
}

TEST (testExpandedSwapMove, multicomponentAllowSwapsNotFullyInserted) {
	const int Mtot = 3;
	std::vector < double > mu (2, 0);
	std::vector < double > ib (3, 10);
	std::vector <int> nmax (2, 3), nmin (2, 0);
	simSystem mysys (2, 1.0, ib, mu, nmax, nmin, Mtot);

	//squareWell sw11, sw12, sw22;
	std::vector < double > params (4, 0.0);
	params[0] = 1.0; params[1] = 0.1; params[2] = 1.0, params[3] = Mtot;
	//sw11.setParameters (params);
	mysys.addPotential (0, 0, "square_well", params, true);
	params[0] = 1.5; params[1] = 0.1; params[2] = NUM_INFINITY, params[3] = Mtot; // once 1-2 come together they won't separate
	//sw12.setParameters (params);
	mysys.addPotential (0, 1, "square_well", params, true);
	params[0] = 2.0; params[1] = 0.1; params[2] = 2.0, params[3] = Mtot;
	//sw22.setParameters (params);
	mysys.addPotential (1, 1, "square_well", params, true);

	atom a1, a2, a3, a4;

	a1.pos[0] = 6.5;
	a1.pos[1] = 0;
	a1.pos[2] = 5;

	a2.pos[0] = 8.01;
	a2.pos[1] = 0;
	a2.pos[2] = 5;

	a3.pos[0] = 0;
	a3.pos[1] = 0;
	a3.pos[2] = 5;

	a4.pos[0] = 2.01;
	a4.pos[1] = 0;
	a4.pos[2] = 5;

	// Fully insert 4 atoms
	mysys.insertAtom(0, &a1);
	mysys.insertAtom(0, &mysys.atoms[0][0]);
	mysys.insertAtom(0, &mysys.atoms[0][0]);

	mysys.insertAtom(0, &a2);
	mysys.insertAtom(0, &mysys.atoms[0][1]);
	mysys.insertAtom(0, &mysys.atoms[0][1]);

	mysys.insertAtom(1, &a3);
	mysys.insertAtom(1, &mysys.atoms[1][0]);
	mysys.insertAtom(1, &mysys.atoms[1][0]);

	mysys.insertAtom(1, &a4);
	mysys.insertAtom(1, &mysys.atoms[1][1]);
	mysys.insertAtom(1, &mysys.atoms[1][1]);

	// partially remove a1 to set it as the partially inserted atom
	mysys.deleteAtom(0, 0);

	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][0]);
	EXPECT_EQ (mysys.getTotN(), 3);
	EXPECT_EQ (mysys.getCurrentM(), 2);

	double U_save = mysys.scratchEnergy();
	EXPECT_TRUE (fabs(U_save - -(0.0+2.0)) < 1.0e-9);

	moves usedMoves (Mtot);
	usedMoves.addSwap(0,1,1.0); // only can swap a1 and a2

	double lastAns = 0;
	int iterMax = 1000, iter = 0;
	while (iter < iterMax) {
		usedMoves.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves.reportMoveStatistics();

		if (ans[0][0] > lastAns) {
			// each time a move succeeds, check that M, N are the same
			EXPECT_EQ (mysys.getTotN(), 3);
			EXPECT_EQ (mysys.getCurrentM(), 2);

			// although the insert/deletes move fractionalAtom around, check that the pointer still points to the one with mState == 2
			EXPECT_TRUE (mysys.getFractionalAtom()->mState == 2);

			// and that no other atoms have that property
			for (unsigned int i = 0; i < 2; ++i) {
				for (unsigned int j = 0; j < 2; ++j) {
					if (&mysys.atoms[i][j] != mysys.getFractionalAtom()) {
						EXPECT_TRUE (mysys.atoms[i][j].mState == 0);
					}
				}
			}
		}
		lastAns = ans[0][0];
		iter++;
	}

	// now insert that atom again and set to a2
	mysys.insertAtom(0, mysys.getFractionalAtom());
	EXPECT_EQ (mysys.getTotN(), 4);
	EXPECT_EQ (mysys.getCurrentM(), 0);

	mysys.deleteAtom(0, 1);
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[0][1]);
	EXPECT_EQ (mysys.getTotN(), 3);
	EXPECT_EQ (mysys.getCurrentM(), 2);

	lastAns = 0;
	iterMax = 1000;
	iter = 0;
	while (iter < iterMax) {
		usedMoves.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves.reportMoveStatistics();

		if (ans[0][0] > lastAns) {
			// each time a move succeeds, check that M, N are the same
			EXPECT_EQ (mysys.getTotN(), 3);
			EXPECT_EQ (mysys.getCurrentM(), 2);

			// although the insert/deletes move fractionalAtom around, check that the pointer still points to the one with mState == 2
			EXPECT_TRUE (mysys.getFractionalAtom()->mState == 2);

			// and that no other atoms have that property
			for (unsigned int i = 0; i < 2; ++i) {
				for (unsigned int j = 0; j < 2; ++j) {
					if (&mysys.atoms[i][j] != mysys.getFractionalAtom()) {
						EXPECT_TRUE (mysys.atoms[i][j].mState == 0);
					}
				}
			}
		}
		lastAns = ans[0][0];
		iter++;
	}

	// now insert that atom and switch to second species
	mysys.insertAtom(0, mysys.getFractionalAtom());
	EXPECT_EQ (mysys.getTotN(), 4);
	EXPECT_EQ (mysys.getCurrentM(), 0);

	mysys.deleteAtom(1, 0);
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[1][0]);
	EXPECT_EQ (mysys.getTotN(), 3);
	EXPECT_EQ (mysys.getCurrentM(), 2);

	lastAns = 0;
	iterMax = 1000;
	iter = 0;
	while (iter < iterMax) {
		usedMoves.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves.reportMoveStatistics();
		if (ans[0][0] > lastAns) {
			// each time a move succeeds, check that M, N are the same
			EXPECT_EQ (mysys.getTotN(), 3);
			EXPECT_EQ (mysys.getCurrentM(), 2);

			// although the insert/deletes move fractionalAtom around, check that the pointer still points to the one with mState == 2
			EXPECT_TRUE (mysys.getFractionalAtom()->mState == 2);

			// and that no other atoms have that property
			for (unsigned int i = 0; i < 2; ++i) {
				for (unsigned int j = 0; j < 2; ++j) {
					if (&mysys.atoms[i][j] != mysys.getFractionalAtom()) {
						EXPECT_TRUE (mysys.atoms[i][j].mState == 0);
					}
				}
			}
		}
		lastAns = ans[0][0];
		iter++;
	}

	// now insert that atom and switch to other of the second species
	mysys.insertAtom(1, mysys.getFractionalAtom());
	EXPECT_EQ (mysys.getTotN(), 4);
	EXPECT_EQ (mysys.getCurrentM(), 0);

	mysys.deleteAtom(1, 1);
	EXPECT_TRUE (mysys.getFractionalAtom() == &mysys.atoms[1][1]);
	EXPECT_EQ (mysys.getTotN(), 3);
	EXPECT_EQ (mysys.getCurrentM(), 2);

	lastAns = 0;
	iterMax = 1000;
	iter = 0;
	while (iter < iterMax) {
		usedMoves.makeMove (mysys);
		std::vector < std::vector < double > > ans = usedMoves.reportMoveStatistics();
		if (ans[0][0] > lastAns) {
			// each time a move succeeds, check that M, N are the same
			EXPECT_EQ (mysys.getTotN(), 3);
			EXPECT_EQ (mysys.getCurrentM(), 2);

			// although the insert/deletes move fractionalAtom around, check that the pointer still points to the one with mState == 2
			EXPECT_TRUE (mysys.getFractionalAtom()->mState == 2);

			// and that no other atoms have that property
			for (unsigned int i = 0; i < 2; ++i) {
				for (unsigned int j = 0; j < 2; ++j) {
					if (&mysys.atoms[i][j] != mysys.getFractionalAtom()) {
						EXPECT_TRUE (mysys.atoms[i][j].mState == 0);
					}
				}
			}
		}
		lastAns = ans[0][0];
		iter++;
	}
}

class testHistogram : public ::testing::Test {
protected:
	histogram *h;
	std::vector < double > lb, ub, coords;
	std::vector < unsigned long long int > dx;

	virtual void SetUp() {
		coords.resize(3, 0);
		lb.resize(3, 10);
		ub.resize(3, 20);
		dx.resize(3, 11);
		h = new histogram (lb, ub, dx);
	}
};

TEST_F (testHistogram, getAddress) {
	coords[0] = 10;
	coords[1] = 10;
	coords[2] = 10;
	EXPECT_EQ(h->getAddress (coords), 0);

	coords[0] = 20;
	coords[1] = 20;
	coords[2] = 20;
	EXPECT_EQ(h->getAddress (coords), 11*11*11-1);

	coords[0] = 11;
	coords[1] = 10;
	coords[2] = 10;
	EXPECT_EQ(h->getAddress (coords), 1);

	coords[0] = 20;
	coords[1] = 10;
	coords[2] = 10;
	EXPECT_EQ(h->getAddress (coords), 10);

	coords[0] = 10;
	coords[1] = 11;
	coords[2] = 10;
	EXPECT_EQ(h->getAddress (coords), 11);

	coords[0] = 20;
	coords[1] = 20;
	coords[2] = 10;
	EXPECT_EQ(h->getAddress (coords), 11*11-1);

	coords[0] = 20;
	coords[1] = 20;
	coords[2] = 11;
	EXPECT_EQ(h->getAddress (coords), 2*11*11-1);

	coords[0] = 15;
	coords[1] = 15;
	coords[2] = 15;
	EXPECT_EQ(h->getAddress (coords), 11*11*5 + 11*5 + 5);

	delete h;
}

TEST_F (testHistogram, getCoords) {
	std::vector <double> coords (3, 0);

	coords = h->getCoords(0);
	EXPECT_EQ (coords[0], 10);
	EXPECT_EQ (coords[1], 10);
	EXPECT_EQ (coords[2], 10);

	coords = h->getCoords(11*11*11-1);
	EXPECT_EQ (coords[0], 20);
	EXPECT_EQ (coords[1], 20);
	EXPECT_EQ (coords[2], 20);

	coords = h->getCoords(1);
	EXPECT_EQ (coords[0], 11);
	EXPECT_EQ (coords[1], 10);
	EXPECT_EQ (coords[2], 10);

	coords = h->getCoords(10);
	EXPECT_EQ (coords[0], 20);
	EXPECT_EQ (coords[1], 10);
	EXPECT_EQ (coords[2], 10);

	coords = h->getCoords(11);
	EXPECT_EQ (coords[0], 10);
	EXPECT_EQ (coords[1], 11);
	EXPECT_EQ (coords[2], 10);

	coords = h->getCoords(11*11-1);
	EXPECT_EQ (coords[0], 20);
	EXPECT_EQ (coords[1], 20);
	EXPECT_EQ (coords[2], 10);

	coords = h->getCoords(2*11*11-1);
	EXPECT_EQ (coords[0], 20);
	EXPECT_EQ (coords[1], 20);
	EXPECT_EQ (coords[2], 11);

	coords = h->getCoords(11*11*5 + 11*5 + 5);
	EXPECT_EQ (coords[0], 15);
	EXPECT_EQ (coords[1], 15);
	EXPECT_EQ (coords[2], 15);

	delete h;
}

TEST_F (testHistogram, incrementCoords) {
	std::vector <double> coords (3, 0), hist1, hist2;

	hist1 = h->getRawHistogram();
	coords = h->getCoords(0);
	h->increment(coords, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[0], 0);
	EXPECT_EQ (hist2[0], 3);

	hist1 = hist2;
	coords = h->getCoords(11*11*11-1);
	h->increment(coords, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[11*11*11-1], 0);
	EXPECT_EQ (hist2[11*11*11-1], 3);

	hist1 = hist2;
	coords = h->getCoords(1);
	h->increment(coords, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[1], 0);
	EXPECT_EQ (hist2[1], 3);

	hist1 = hist2;
	coords = h->getCoords(10);
	h->increment(coords, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[10], 0);
	EXPECT_EQ (hist2[10], 3);

	hist1 = hist2;
	coords = h->getCoords(11);
	h->increment(coords, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[11], 0);
	EXPECT_EQ (hist2[11], 3);

	hist1 = hist2;
	coords = h->getCoords(11*11-1);
	h->increment(coords, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[11*11-1], 0);
	EXPECT_EQ (hist2[11*11-1], 3);

	hist1 = hist2;
	coords = h->getCoords(2*11*11-1);
	h->increment(coords, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[2*11*11-1], 0);
	EXPECT_EQ (hist2[2*11*11-1], 3);

	hist1 = hist2;
	coords = h->getCoords(11*11*5 + 11*5 + 5);
	h->increment(coords, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[11*11*5 + 11*5 + 5], 0);
	EXPECT_EQ (hist2[11*11*5 + 11*5 + 5], 3);

	delete h;
}

TEST_F (testHistogram, incrementAddress) {
	std::vector <double> hist1, hist2;

	hist1 = h->getRawHistogram();
	h->increment(0, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[0], 0);
	EXPECT_EQ (hist2[0], 3);

	hist1 = hist2;
	h->increment(11*11*11-1, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[11*11*11-1], 0);
	EXPECT_EQ (hist2[11*11*11-1], 3);

	hist1 = hist2;
	h->increment(1, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[1], 0);
	EXPECT_EQ (hist2[1], 3);

	hist1 = hist2;;
	h->increment(10, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[10], 0);
	EXPECT_EQ (hist2[10], 3);

	hist1 = hist2;
	h->increment(11, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[11], 0);
	EXPECT_EQ (hist2[11], 3);

	hist1 = hist2;
	h->increment(11*11-1, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[11*11-1], 0);
	EXPECT_EQ (hist2[11*11-1], 3);

	hist1 = hist2;
	h->increment(2*11*11-1, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[2*11*11-1], 0);
	EXPECT_EQ (hist2[2*11*11-1], 3);

	hist1 = hist2;
	h->increment(11*11*5 + 11*5 + 5, 3);
	hist2 = h->getRawHistogram();
	EXPECT_EQ (hist1[11*11*5 + 11*5 + 5], 0);
	EXPECT_EQ (hist2[11*11*5 + 11*5 + 5], 3);

	delete h;
}

class testHardWallZ : public ::testing::Test {
protected:
    atom a1;
    double sigma, H;
    int M;
    std::vector < double > box;

    virtual void SetUp() {
        std::vector < double > coords (3, 0);
        sigma = 1.234;
        coords[2] = sigma/2.0;
        a1.pos = coords;
        H = 3*sigma;
        M = 3;
        box.resize(3, 4*sigma);
    }
};

TEST_F (testHardWallZ, badInit) {
    hardWallZ *hwz;

    // ub > lb
    bool caught = false;
    try {
        hwz = new hardWallZ (H, 0, sigma, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete hwz;
	}

    // bad sigma
    caught = false;
    try {
        hwz = new hardWallZ (0, H, -sigma, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete hwz;
	}

    // bad M
    caught = false;
    try {
        hwz = new hardWallZ (0, H, sigma, 0);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete hwz;
	}
}

TEST_F (testHardWallZ, inside) {
    hardWallZ hwz (0, H, sigma); // M defaults to 1

    double zpos = a1.pos[2]; // at the border

    // test bottom wall
    a1.pos[2] = zpos*1.01; // inside (not overlapping wall)
    bool inside = hwz.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = zpos*0.99; // outside (overlapping wall)
    inside = hwz.inside (&a1, box);
    EXPECT_TRUE (!inside);

    // test top wall
    a1.pos[2] = H-zpos*1.01; // inside (not overlapping wall)
    inside = hwz.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = H-zpos*0.99; // outside (overlapping wall)
    inside = hwz.inside (&a1, box);
    EXPECT_TRUE (!inside);
}

TEST_F (testHardWallZ, energy) {
    hardWallZ hwz (0, H, sigma); // M defaults to 1

    double zpos = a1.pos[2]; // at the border

    // test bottom wall
    a1.pos[2] = zpos*1.01; // inside (not overlapping wall)
    double U = hwz.energy (&a1, box);
    EXPECT_EQ (U, 0);

    a1.pos[2] = zpos*0.99; // outside (overlapping wall)
    U = hwz.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    // test top wall
    a1.pos[2] = H-zpos*1.01; // inside (not overlapping wall)
    U = hwz.energy (&a1, box);
    EXPECT_EQ (U, 0);

    a1.pos[2] = H-zpos*0.99; // outside (overlapping wall)
    U = hwz.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);
}

TEST_F (testHardWallZ, insideM) {
    hardWallZ hwz (0, H, sigma, M);
    a1.mState = M-1;

    // test bottom wall
    a1.pos[2] = sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    bool inside = hwz.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    inside = hwz.inside (&a1, box);
    EXPECT_TRUE (!inside);

    // test top wall
    a1.pos[2] = H-sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    inside = hwz.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = H-sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    inside = hwz.inside (&a1, box);
    EXPECT_TRUE (!inside);
}

TEST_F (testHardWallZ, energyM) {
    hardWallZ hwz (0, H, sigma, M);
    a1.mState = M-1;

    // test bottom wall
    a1.pos[2] = sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    double U = hwz.energy (&a1, box);
    EXPECT_EQ (U, 0);

    a1.pos[2] = sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    U = hwz.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    // test top wall
    a1.pos[2] = H-sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    U = hwz.energy (&a1, box);
    EXPECT_EQ (U, 0);

    a1.pos[2] = H-sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    U = hwz.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);
}

class testSquareWellWallZ : public ::testing::Test {
protected:
    atom a1;
    double eps, sigma, H, range;
    int M;
    std::vector < double > box;

    virtual void SetUp() {
        std::vector < double > coords (3, 0);
        sigma = 1.234;
        coords[2] = sigma/2.0;
        a1.pos = coords;
        eps = 2.345;
        H = 3*sigma;
        M = 3;
        range = sigma;
        box.resize(3, 4*sigma);
    }
};

TEST_F (testSquareWellWallZ, badInit) {
    squareWellWallZ *sqz;

    // ub > lb
    bool caught = false;
    try {
        sqz = new squareWellWallZ (H, 0, sigma, range, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete sqz;
	}

    // bad sigma
    caught = false;
    try {
        sqz = new squareWellWallZ (0, H, -sigma, range, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete sqz;
	}

    caught = false;
    try {
        sqz = new squareWellWallZ (0, H, range*2.0, range, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete sqz;
	}

    // bad M
    caught = false;
    try {
        sqz = new squareWellWallZ (0, H, sigma, range, eps, 0);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete sqz;
	}

    // bad range
    caught = false;
    try {
        sqz = new squareWellWallZ (0, H, sigma, -range, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete sqz;
	}

    // bad eps
    caught = false;
    try {
        sqz = new squareWellWallZ (0, H, sigma, range, -eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete sqz;
	}
}

TEST_F (testSquareWellWallZ, inside) {
    squareWellWallZ sqz (0, H, sigma, range, eps); // M defaults to 1

    double zpos = a1.pos[2]; // at the border

    // test bottom wall
    a1.pos[2] = zpos*1.01; // inside (not overlapping wall)
    bool inside = sqz.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = zpos*0.99; // outside (overlapping wall)
    inside = sqz.inside (&a1, box);
    EXPECT_TRUE (!inside);

    // test top wall
    a1.pos[2] = H-zpos*1.01; // inside (not overlapping wall)
    inside = sqz.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = H-zpos*0.99; // outside (overlapping wall)
    inside = sqz.inside (&a1, box);
    EXPECT_TRUE (!inside);
}

TEST_F (testSquareWellWallZ, energy) {
    squareWellWallZ sqz (0, H, sigma, range, eps); // M defaults to 1

    double zpos = a1.pos[2]; // at the border

    // test bottom wall
    a1.pos[2] = zpos*1.01; // inside (not overlapping wall)
    double U = sqz.energy (&a1, box);
    EXPECT_EQ (U, -eps);

    a1.pos[2] = zpos*0.99; // outside (overlapping wall)
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = range*0.99; // inside range
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, -eps);

    a1.pos[2] = range*1.01; // outside range
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, 0);

    // test top wall
    a1.pos[2] = H-zpos*1.01; // inside (not overlapping wall)
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, -eps);

    a1.pos[2] = H-zpos*0.99; // outside (overlapping wall)
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = H-range*0.99; // inside range
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, -eps);

    a1.pos[2] = H-range*1.01; // outside range
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, 0);
}

TEST_F (testSquareWellWallZ, insideM) {
    squareWellWallZ sqz (0, H, sigma, range, eps, M);
    a1.mState = M-1;

    // test bottom wall
    a1.pos[2] = sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    bool inside = sqz.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    inside = sqz.inside (&a1, box);
    EXPECT_TRUE (!inside);

    // test top wall
    a1.pos[2] = H-sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    inside = sqz.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = H-sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    inside = sqz.inside (&a1, box);
    EXPECT_TRUE (!inside);
}

TEST_F (testSquareWellWallZ, energyM) {
    squareWellWallZ sqz (0, H, sigma, range, eps, M);
    a1.mState = M-1;

    // test bottom wall
    a1.pos[2] = sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    double U = sqz.energy (&a1, box);
    EXPECT_EQ (U, -eps/M*a1.mState);

    a1.pos[2] = sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = range*0.99; // inside range
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, -eps/M*a1.mState);

    a1.pos[2] = range*1.01; // outside range
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, 0);

    // test top wall
    a1.pos[2] = H-sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, -eps/M*a1.mState);

    a1.pos[2] = H-sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = H-range*0.99; // inside range
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, -eps/M*a1.mState);

    a1.pos[2] = H-range*1.01; // outside range
    U = sqz.energy (&a1, box);
    EXPECT_EQ (U, 0);
}

//  then add 2 square well walls with diff ranges
// also test M

class testCompositeBarrier : public ::testing::Test {
protected:
    compositeBarrier cB;
    atom a1;
    double eps, sigma, H, range1, range2;
    int M;
    std::vector < double > box;

    virtual void SetUp() {
        std::vector < double > coords (3, 0);
        sigma = 1.234;
        coords[2] = sigma/2.0;
        a1.pos = coords;
        eps = 2.345;
        H = 3*sigma;
        M = 3;
        range1 = sigma;
        range2 = 0.75*sigma;
        box.resize(3, 4*sigma);
    }
};

TEST_F (testCompositeBarrier, empty) {
    const int nincrs = 5;
    for (unsigned int i = 0; i < nincrs; ++i) {
        double zpos = H/(nincrs-1)*i;
        a1.pos[2] = zpos;
        bool inside = cB.inside(&a1, box);
        EXPECT_TRUE (inside);
    }
}

TEST_F (testCompositeBarrier, singleHardWall) {
    cB.addHardWallZ (0, H, sigma); // M defaults to 1

    double zpos = a1.pos[2]; // at the border

    // test bottom wall
    a1.pos[2] = zpos*1.01; // inside (not overlapping wall)
    bool inside = cB.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = zpos*0.99; // outside (overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (!inside);

    // test top wall
    a1.pos[2] = H-zpos*1.01; // inside (not overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = H-zpos*0.99; // outside (overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (!inside);
}

TEST_F (testCompositeBarrier, singleHardWallM) {
    cB.addHardWallZ (0, H, sigma, M);

    a1.mState = M-1;

    // test bottom wall
    a1.pos[2] = sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    bool inside = cB.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (!inside);

    // test top wall
    a1.pos[2] = H-sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = H-sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (!inside);
}

TEST_F (testCompositeBarrier, singleSquareWellWallZ) {
    cB.addSquareWellWallZ (0, H, sigma, range1, eps); // M defaults to 1

    double zpos = a1.pos[2]; // at the border

    // test bottom wall
    a1.pos[2] = zpos*1.01; // inside (not overlapping wall)
    bool inside = cB.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = zpos*0.99; // outside (overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (!inside);

    // test top wall
    a1.pos[2] = H-zpos*1.01; // inside (not overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = H-zpos*0.99; // outside (overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (!inside);
}

TEST_F (testCompositeBarrier, singleSquareWellWallZM) {
    cB.addSquareWellWallZ (0, H, sigma, range1, eps, M);
    a1.mState = M-1;

    // test bottom wall
    a1.pos[2] = sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    bool inside = cB.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (!inside);

    // test top wall
    a1.pos[2] = H-sigma/2.0*a1.mState/M*1.01; // inside (not overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (inside);

    a1.pos[2] = H-sigma/2.0*a1.mState/M*0.99; // outside (overlapping wall)
    inside = cB.inside (&a1, box);
    EXPECT_TRUE (!inside);
}

TEST_F (testCompositeBarrier, pairSquareWellWallZ) {
    cB.addSquareWellWallZ (0, H, sigma, range1, eps); // M defaults to 1
    cB.addSquareWellWallZ (0, H, sigma, range2, eps); // M defaults to 1

    double zpos = a1.pos[2]; // at the border

    // test bottom wall
    a1.pos[2] = zpos*1.01; // inside (not overlapping wall)
    double U = cB.energy (&a1, box);
    EXPECT_EQ (U, -2.0*eps);

    a1.pos[2] = zpos*0.99; // outside (overlapping wall)
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = range1*0.99; // inside single range
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, -eps);

    a1.pos[2] = range1*1.01; // outside ranges
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, 0);

    a1.pos[2] = range2*0.99; // inside both ranges
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, -2.0*eps);

    a1.pos[2] = range2*1.01; // outside single range, but inside other
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, -eps);

    // test top wall
    a1.pos[2] = H-zpos*1.01; // inside (not overlapping wall)
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, -2.0*eps);

    a1.pos[2] = H-zpos*0.99; // outside (overlapping wall)
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = H-range1*0.99; // inside single range
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, -eps);

    a1.pos[2] = H-range1*1.01; // outside ranges
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, 0);

    a1.pos[2] = H-range2*0.99; // inside both ranges
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, -2.0*eps);

    a1.pos[2] = H-range2*1.01; // outside single range, but inside other
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, -eps);
}

TEST_F (testCompositeBarrier, hardWallInsideSquareWell) {
    cB.addSquareWellWallZ (0, H, sigma, range1, eps); // M defaults to 1
    cB.addHardWallZ (H/2-sigma/2*1.01, H/2+sigma/2*1.01, sigma); // M defaults to 1

    // test bottom wall
    a1.pos[2] = 0;
    double U = cB.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = H/3.0;
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = H/2.0;
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, 0);

    a1.pos[2] = H*2.0/3.0;
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);

    a1.pos[2] = H;
    U = cB.energy (&a1, box);
    EXPECT_EQ (U, NUM_INFINITY);
}

class testRightTriangleXZ : public ::testing::Test {
protected:
    rightTriangleXZ *rtz;
    atom a1;

    double width, theta, lamW, eps, sigma, sep, offset, L, zbase;
    std::vector < double > box;
    bool top, inside;
    int M;

    virtual void SetUp() {
        sigma = 1.234;
        eps = 1.234;
        lamW = 2.345;
        L = 10;
        box.resize(3, L);
        sep = 0;
        width = (L/2.0 - sep);
        M = 3;
        top = false;
        offset = 0;
        theta = PI/4.0;
        zbase = 0;
    }
};

TEST_F (testRightTriangleXZ, badSep) {
    // sep < 0
    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, -0.01, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
}

TEST_F (testRightTriangleXZ, badPeriodicity) {
    // bad periodicity
    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep+0.01, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
}

TEST_F (testRightTriangleXZ, badZbase) {
    // zbase < 0
    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, -0.01, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);

    // zbase > Z-width of box
    caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, box[2]+0.01, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
}

TEST_F (testRightTriangleXZ, badTheta) {
    // theta <= 0
    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, 0.0, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);

    caught = false;
    try {
        rtz = new rightTriangleXZ (width, -0.01, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);

    // theta >= 90 degrees (PI/2)
    caught = false;
    try {
        rtz = new rightTriangleXZ (width, PI/2.0, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);

    caught = false;
    try {
        rtz = new rightTriangleXZ (width, PI/2.0*1.1, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
}

TEST_F (testRightTriangleXZ, badSigma) {
    // sigma <= 0
    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, 0.0, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
}

TEST_F (testRightTriangleXZ, badLamW) {
    // lamW < 1
    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, 0.99, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
}

TEST_F (testRightTriangleXZ, badEps) {
    // eps < 0
    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, -0.01, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
}

TEST_F (testRightTriangleXZ, badM) {
    // M <= 0
    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, 0);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);

    caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, -1);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_0_M_0) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.7234325 + image*box[2];

        a1.pos[2] = 2.47564552135*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.47564552135*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.92371764031*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.92371764031*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.14232755646*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.14232755646*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_0_M_1) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.241144166667 + image*box[2];

        a1.pos[2] = 0.829007281949*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.829007281949*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 0.97457254677*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.97457254677*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.34682728195*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.34682728195*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_0_M_2) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.482288333333 + image*box[2];

        a1.pos[2] = 1.6580145639*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.6580145639*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 1.94914509354*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.94914509354*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.37159404395*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.37159404395*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_1_M_0) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.875 + image*box[2];

        a1.pos[2] = 1.79498198691*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.79498198691*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 4.42392334365*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.42392334365*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_1_M_1) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false, inside;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.875 + image*box[2];

        a1.pos[2] = 1.32001516546*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.32001516546*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.19632895104*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.19632895104*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_1_M_2) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.875 + image*box[2];

        a1.pos[2] = 1.55749857618*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.55749857618*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.31012614735*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.31012614735*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_2_M_0) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.7875 + image*box[2];

        a1.pos[2] = 2.60669076301*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.60669076301*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 4.27337279813*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.27337279813*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_2_M_1) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.7875 + image*box[2];

        a1.pos[2] = 2.30557386417*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.30557386417*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.86113454254*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.86113454254*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_2_M_2) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.7875 + image*box[2];

        a1.pos[2] = 2.45613231359*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.45613231359*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.56725367034*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.56725367034*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_3_M_0) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 4.375 + image*box[2];

        a1.pos[2] = 2.31653175473*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.31653175473*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 4.11595264758*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.11595264758*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_3_M_1) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 4.375 + image*box[2];

        a1.pos[2] = 1.49386508806*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.49386508806*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.70371439199*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70371439199*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_3_M_2) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 4.375 + image*box[2];

        a1.pos[2] = 1.9051984214*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.9051984214*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.40983351978*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.40983351978*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_4_M_0) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 6.25302184585 + image*box[2];

        a1.pos[2] = 1.43588273218*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.43588273218*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.61273821062*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.61273821062*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.06482408892*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.06482408892*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_4_M_1) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.41767394862 + image*box[2];

        a1.pos[2] = 0.478627577393*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.478627577393*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 1.20572083333*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.20572083333*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.35494136297*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.35494136297*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_4_M_2) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.83534789723 + image*box[2];

        a1.pos[2] = 0.957255154787*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.957255154787*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.41144166667*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.41144166667*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70988272595*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70988272595*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_0) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -1.25302184585 + image*box[2];

        a1.pos[2] = 1.43588273218*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.43588273218*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.61273821062*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.61273821062*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.06482408892*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.06482408892*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_1) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.417673948616 + image*box[2];

        a1.pos[2] = 0.478627577393*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.478627577393*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 1.20572083333*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.20572083333*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.35494136297*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.35494136297*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_2) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.835347897231 + image*box[2];

        a1.pos[2] = 0.957255154787*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.957255154787*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.41144166667*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.41144166667*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70988272595*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70988272595*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_0) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2];

        a1.pos[2] = 2.31653175473*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.31653175473*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.98054767068*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.98054767068*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.11595264758*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.11595264758*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_1) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2];

        a1.pos[2] = 1.49386508806*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.49386508806*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.70371439199*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70371439199*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_2) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2];

        a1.pos[2] = 1.9051984214*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.9051984214*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.40983351978*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.40983351978*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_0) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2];

        a1.pos[2] = 2.62008822264*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.62008822264*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 4.28677025775*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.28677025775*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_1) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2];

        a1.pos[2] = 2.31897132379*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.31897132379*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.87453200216*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.87453200216*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_2) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2];

        a1.pos[2] = 2.46952977321*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.46952977321*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.58065112996*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.58065112996*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_0) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2];

        a1.pos[2] = 1.79498198691*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.79498198691*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.53986473581*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.53986473581*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.42392334365*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.42392334365*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_1) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2];

        a1.pos[2] = 1.32001516546*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.32001516546*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.19632895104*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.19632895104*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_2) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2];

        a1.pos[2] = 1.55749857618*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.55749857618*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.31012614735*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.31012614735*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_0) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.7234325 + image*box[2];

        a1.pos[2] = 2.47564552135*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.47564552135*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.92371764031*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.92371764031*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.14232755646*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.14232755646*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_1) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.24114416667 + image*box[2];

        a1.pos[2] = 0.829007281949*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.829007281949*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 0.97457254677*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.97457254677*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.34682728195*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.34682728195*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_2) {
    theta = PI/3.0;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.48228833333 + image*box[2];

        a1.pos[2] = 1.6580145639*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.6580145639*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 1.94914509354*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.94914509354*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.37159404395*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.37159404395*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_0_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -1.25302184585 + image*box[2];

        a1.pos[2] = box[2] - 1.43588273218*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.43588273218*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 3.61273821062*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 3.61273821062*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.06482408892*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.06482408892*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_1_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.417673948616 + image*box[2];

        a1.pos[2] = box[2] - 0.478627577393*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 0.478627577393*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 1.20572083333*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.20572083333*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.35494136297*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.35494136297*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_2_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.835347897231 + image*box[2];

        a1.pos[2] = box[2] - 0.957255154787*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 0.957255154787*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 2.41144166667*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.41144166667*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.70988272595*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.70988272595*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_0_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2];

        a1.pos[2] = box[2] - 2.31653175473*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.31653175473*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 2.98054767068*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.98054767068*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.11595264758*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.11595264758*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_1_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2];

        a1.pos[2] = box[2] - 1.49386508806*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.49386508806*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 2.70371439199*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.70371439199*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_2_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2];

        a1.pos[2] = box[2] - 1.9051984214*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.9051984214*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 3.40983351978*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 3.40983351978*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_0_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2];

        a1.pos[2] = box[2] - 2.62008822264*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.62008822264*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 4.28677025775*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.28677025775*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_1_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2];

        a1.pos[2] = box[2] - 2.31897132379*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.31897132379*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 2.87453200216*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.87453200216*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_2_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2];

        a1.pos[2] = box[2] - 2.46952977321*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.46952977321*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 3.58065112996*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 3.58065112996*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_0_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2];

        a1.pos[2] = box[2] - 1.79498198691*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.79498198691*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 2.53986473581*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.53986473581*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.42392334365*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.42392334365*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_1_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2];

        a1.pos[2] = box[2] - 1.32001516546*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.32001516546*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 2.19632895104*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.19632895104*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_2_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2];

        a1.pos[2] = box[2] - 1.55749857618*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.55749857618*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 3.31012614735*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 3.31012614735*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_0_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.7234325 + image*box[2];

        a1.pos[2] = box[2] - 2.47564552135*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.47564552135*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 2.92371764031*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.92371764031*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.14232755646*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 4.14232755646*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_1_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.24114416667 + image*box[2];

        a1.pos[2] = box[2] - 0.829007281949*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 0.829007281949*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 0.97457254677*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 0.97457254677*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.34682728195*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 2.34682728195*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_2_top) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.48228833333 + image*box[2];

        a1.pos[2] = box[2] - 1.6580145639*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.6580145639*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2] - 1.94914509354*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.94914509354*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 3.37159404395*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 3.37159404395*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_0_top_z_periodicity) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.7234325 + image*box[2];

        a1.pos[2] = 3*box[2] - 2.47564552135*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3*box[2] - 2.47564552135*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3*box[2] - 2.92371764031*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3*box[2] - 2.92371764031*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3*box[2] - 4.14232755646*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3*box[2] - 4.14232755646*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_1_top_z_periodicity) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.24114416667 + image*box[2];

        a1.pos[2] = -3*box[2] - 0.829007281949*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = -3*box[2] - 0.829007281949*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = -3*box[2] - 0.97457254677*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = -3*box[2] - 0.97457254677*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = -3*box[2] - 2.34682728195*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = -3*box[2] - 2.34682728195*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_2_top_z_periodicity) {
    theta = PI/3.0;
    top = true;
    zbase = box[2];
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.48228833333 + image*box[2];

        a1.pos[2] = 3*box[2] - 1.6580145639*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = -3*box[2] - 1.6580145639*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3*box[2] - 1.94914509354*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = -3*box[2] - 1.94914509354*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3*box[2] - 3.37159404395*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = -3*box[2] - 3.37159404395*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_0_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -1.25302184585 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 1.43588273218*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.43588273218*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 3.61273821062*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 3.61273821062*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.06482408892*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.06482408892*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_1_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.417673948616 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 0.478627577393*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 0.478627577393*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 1.20572083333*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.20572083333*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.35494136297*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.35494136297*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_0_M_2_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.835347897231 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 0.957255154787*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 0.957255154787*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 2.41144166667*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.41144166667*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.70988272595*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.70988272595*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_0_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 2.31653175473*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.31653175473*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 2.98054767068*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.98054767068*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.11595264758*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.11595264758*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_1_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 1.49386508806*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.49386508806*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 2.70371439199*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.70371439199*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_1_M_2_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 0.625 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 1.9051984214*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.9051984214*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 3.40983351978*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 3.40983351978*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_0_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 2.62008822264*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.62008822264*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 4.28677025775*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.28677025775*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_1_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 2.31897132379*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.31897132379*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 2.87453200216*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.87453200216*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_2_M_2_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.2625 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 2.46952977321*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.46952977321*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 3.58065112996*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 3.58065112996*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_0_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 1.79498198691*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2] - 1.79498198691*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 2.53986473581*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.53986473581*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.42392334365*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.42392334365*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_1_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 1.32001516546*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.32001516546*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 2.19632895104*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.19632895104*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_3_M_2_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.125 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 1.55749857618*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.55749857618*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 3.31012614735*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 3.31012614735*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_0_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.7234325 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 2.47564552135*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.47564552135*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 2.92371764031*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.92371764031*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.14232755646*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 4.14232755646*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_1_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.24114416667 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 0.829007281949*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 0.829007281949*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 0.97457254677*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 0.97457254677*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.34682728195*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 2.34682728195*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 60_degrees_energy_pos_4_M_2_top_zbase_offset) {
    theta = PI/3.0;
    top = true;
    zbase = box[2]*0.75;
    offset = 2.345;
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.48228833333 + image*box[2] + offset;

        a1.pos[2] = box[2]*0.75 - 1.6580145639*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.6580145639*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = box[2]*0.75 - 1.94914509354*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 1.94914509354*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 3.37159404395*1.01;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = box[2]*0.75 - 3.37159404395*0.99;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_0_M_0_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.7234325 + image*box[2] + offset;

        a1.pos[2] = 2.47564552135*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.47564552135*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.92371764031*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.92371764031*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.14232755646*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.14232755646*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_0_M_1_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.241144166667 + image*box[2] + offset;

        a1.pos[2] = 0.829007281949*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.829007281949*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 0.97457254677*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.97457254677*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.34682728195*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.34682728195*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_0_M_2_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = -0.482288333333 + image*box[2] + offset;

        a1.pos[2] = 1.6580145639*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.6580145639*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 1.94914509354*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.94914509354*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.37159404395*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.37159404395*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_1_M_0_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.875 + image*box[2] + offset;

        a1.pos[2] = 1.79498198691*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.79498198691*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 4.42392334365*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.42392334365*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_1_M_1_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false, inside;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.875 + image*box[2] + offset;

        a1.pos[2] = 1.32001516546*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.32001516546*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.19632895104*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.19632895104*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_1_M_2_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 1.875 + image*box[2] + offset;

        a1.pos[2] = 1.55749857618*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.55749857618*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.31012614735*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.31012614735*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_2_M_0_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.7875 + image*box[2] + offset;

        a1.pos[2] = 2.60669076301*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.60669076301*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 4.27337279813*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.27337279813*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_2_M_1_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.7875 + image*box[2] + offset;

        a1.pos[2] = 2.30557386417*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.30557386417*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.86113454254*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.86113454254*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_2_M_2_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 3.7875 + image*box[2] + offset;

        a1.pos[2] = 2.45613231359*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.45613231359*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.56725367034*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.56725367034*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_3_M_0_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 2; ++image) {
        a1.pos[0] = 4.375 + image*box[2] + offset;

        a1.pos[2] = 2.31653175473*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.31653175473*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 4.11595264758*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.11595264758*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_3_M_1_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 4.375 + image*box[2] + offset;

        a1.pos[2] = 1.49386508806*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.49386508806*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.70371439199*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70371439199*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_3_M_2_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 4.375 + image*box[2] + offset;

        a1.pos[2] = 1.9051984214*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.9051984214*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.40983351978*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.40983351978*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_4_M_0_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 6.25302184585 + image*box[2] + offset;

        a1.pos[2] = 1.43588273218*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.43588273218*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 3.61273821062*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 3.61273821062*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.06482408892*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 4.06482408892*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_4_M_1_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 1;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.41767394862 + image*box[2] + offset;

        a1.pos[2] = 0.478627577393*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.478627577393*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 1.20572083333*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.20572083333*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps/3.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.35494136297*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 1.35494136297*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps/3.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

TEST_F (testRightTriangleXZ, 30_degrees_energy_pos_4_M_2_zbase_offset) {
    theta = PI/6.0;
    double U = 0;
    a1.pos[1] = 0;
    zbase = 1.234;
    offset = 2.345;

    bool caught = false;
    try {
        rtz = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);

    a1.mState = 2;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {
        a1.pos[0] = 5.83534789723 + image*box[2] + offset;

        a1.pos[2] = 0.957255154787*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 0.957255154787*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (!inside);

        a1.pos[2] = 2.41144166667*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.41144166667*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -2*eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70988272595*1.01 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, 0.0);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);

        a1.pos[2] = 2.70988272595*0.99 + zbase;
        U = rtz->energy(&a1, box);
        EXPECT_EQ (U, -eps*2./3.);
        inside = rtz->inside(&a1, box);
        EXPECT_TRUE (inside);
    }

    delete rtz;
}

class testCompositeRightTriangleXZ : public ::testing::Test {
protected:
    compositeBarrier cB;
    atom a1;

    double width, theta, lamW, eps, sigma, sep, offset, L, zbase;
    std::vector < double > box;
    bool top, inside;
    int M;

    virtual void SetUp() {
        sigma = 1.234;
        eps = 1.234;
        lamW = 2.345;
        L = 10;
        box.resize(3, L);
        sep = 0;
        width = (L/2.0 - sep);
        M = 3;
        top = false;
        offset = 0;
        theta = PI/3.0;
        zbase = 0;
    }
};

TEST_F (testCompositeRightTriangleXZ, top_bottom) {
    double U = 0;
    a1.pos[1] = 0;

    bool caught = false;
    try {
        cB.addRightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE(!caught);

    caught = false;
    try {
        cB.addRightTriangleXZ (width, theta, lamW, eps, sigma, sep, 2.0, box, box[2], true, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE(!caught);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {

        // chosen from visualization
        a1.pos[0] = 1.447489 + image*box[2];
        a1.pos[2] = 2.230522;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 3.940259 + image*box[2];
        a1.pos[2] = 0.741173;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 6.719604 + image*box[2];
        a1.pos[2] = 2.580019;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 0.405736 + image*box[2];
        a1.pos[2] = 9.447467;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 4.329733 + image*box[2];
        a1.pos[2] = 8.129570;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 8.543460 + image*box[2];
        a1.pos[2] = 7.317294;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 1.314948 + image*box[2];
        a1.pos[2] = 5.024857;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 4.379501 + image*box[2];
        a1.pos[2] = 4.112680;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 6.185803 + image*box[2];
        a1.pos[2] = 5.944372;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 0.513676 + image*box[2];
        a1.pos[2] = 3.656887;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 3.623591 + image*box[2];
        a1.pos[2] = 3.435689;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 9.915614 + image*box[2];
        a1.pos[2] = 3.731327;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 0.169984 + image*box[2];
        a1.pos[2] = 7.152836;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 4.721392 + image*box[2];
        a1.pos[2] = 5.679895;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 8.310266 + image*box[2];
        a1.pos[2] = 7.118589;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 1.142391 + image*box[2];
        a1.pos[2] = 6.530009;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 6.791347 + image*box[2];
        a1.pos[2] = 8.781039;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 0.269557 + image*box[2];
        a1.pos[2] = 2.405319;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 4.214802 + image*box[2];
        a1.pos[2] = 3.517564;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 9.671770 + image*box[2];
        a1.pos[2] = 1.179469;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

    }
}

TEST_F (testCompositeRightTriangleXZ, bottom_bottom) {
    double U = 0;
    a1.pos[1] = 0;

    cB.addRightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    cB.addRightTriangleXZ (width, theta, lamW, eps, sigma, sep, 2.0, box, zbase, top, M);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {

        // chosen from visualization
        a1.pos[0] = 0.378401 + image*box[2];
        a1.pos[2] = 1.753589;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 1.186937 + image*box[2];
        a1.pos[2] = 2.398154;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 2.644928 + image*box[2];
        a1.pos[2] = 2.069054;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 5.183863 + image*box[2];
        a1.pos[2] = 1.699569;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 6.974709 + image*box[2];
        a1.pos[2] = 2.412507;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 7.946108 + image*box[2];
        a1.pos[2] = 2.479060;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 0.791921 + image*box[2];
        a1.pos[2] = 4.314692;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 2.892660 + image*box[2];
        a1.pos[2] = 4.943779;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 3.489028 + image*box[2];
        a1.pos[2] = 4.528201;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 6.885559 + image*box[2];
        a1.pos[2] = 4.685033;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 9.849988 + image*box[2];
        a1.pos[2] = 4.700240;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 0.212752 + image*box[2];
        a1.pos[2] = 4.289830;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 1.264140 + image*box[2];
        a1.pos[2] = 4.131850;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 3.076959 + image*box[2];
        a1.pos[2] = 4.441474;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 5.832169 + image*box[2];
        a1.pos[2] = 4.127904;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 9.559871 + image*box[2];
        a1.pos[2] = 4.602983;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 9.181351 + image*box[2];
        a1.pos[2] = 3.913560;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 0.553184 + image*box[2];
        a1.pos[2] = 3.360072;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 1.144858 + image*box[2];
        a1.pos[2] = 3.826946;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 0.581205 + image*box[2];
        a1.pos[2] = 3.949698;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 3.257825 + image*box[2];
        a1.pos[2] = 2.978000;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 4.902641 + image*box[2];
        a1.pos[2] = 3.822539;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 6.712595 + image*box[2];
        a1.pos[2] = 3.839948;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 8.125311 + image*box[2];
        a1.pos[2] = 4.196345;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 9.710505 + image*box[2];
        a1.pos[2] = 3.775819;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -2*eps);

        a1.pos[0] = 0.034884 + image*box[2];
        a1.pos[2] = 3.223501;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -3*eps);

        a1.pos[0] = 0.786501 + image*box[2];
        a1.pos[2] = 3.577915;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -3*eps);

        a1.pos[0] = 1.262904 + image*box[2];
        a1.pos[2] = 3.626476;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -3*eps);

        a1.pos[0] = 2.524760 + image*box[2];
        a1.pos[2] = 2.283315;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -3*eps);

        a1.pos[0] = 8.008121 + image*box[2];
        a1.pos[2] = 2.609876;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -3*eps);

        a1.pos[0] = 0.523231 + image*box[2];
        a1.pos[2] = 2.871053;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -4*eps);

        a1.pos[0] = 0.358148 + image*box[2];
        a1.pos[2] = 1.888662;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -4*eps);

        a1.pos[0] = 1.007053 + image*box[2];
        a1.pos[2] = 2.714060;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -4*eps);

        a1.pos[0] = 5.029243 + image*box[2];
        a1.pos[2] = 2.038298;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -4*eps);

        a1.pos[0] = 5.407055 + image*box[2];
        a1.pos[2] = 2.999161;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -4*eps);

        a1.pos[0] = 5.937607 + image*box[2];
        a1.pos[2] = 2.724622;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -4*eps);
    }
}


TEST_F (testCompositeRightTriangleXZ, bottom_bottom_sep) {
    double U = 0;
    a1.pos[1] = 0;

    width = 2.0;
    sep = 8.0;
    box[0] = 20;
    box[1] = 20;
    box[2] = 20;
    lamW = 1.1;

    cB.addRightTriangleXZ (width, theta, lamW, eps, sigma, sep, 0.0, box, zbase, top, M);
    cB.addRightTriangleXZ (width, theta, lamW, eps, sigma, sep, 4.5, box, zbase, top, M);

    a1.mState = 0;

    // check periodicity outside of box
    for (int image = -2; image < 3; ++image) {

        // chosen from visualization
        a1.pos[0] = 0.371768 + image*box[2];
        a1.pos[2] = 0.287454;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 1.054832 + image*box[2];
        a1.pos[2] = 1.250588;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 5.460271 + image*box[2];
        a1.pos[2] = 1.241101;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 6.768400 + image*box[2];
        a1.pos[2] = 0.550438;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 9.822715 + image*box[2];
        a1.pos[2] = 0.408006;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 10.068356 + image*box[2];
        a1.pos[2] = 1.038385;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 15.347902 + image*box[2];
        a1.pos[2] = 0.616920;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 11.825713 + image*box[2];
        a1.pos[2] = 0.434620;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, NUM_INFINITY);

        a1.pos[0] = 0.312175 + image*box[2];
        a1.pos[2] = 2.088335;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 1.973930 + image*box[2];
        a1.pos[2] = 1.971786;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 2.978860 + image*box[2];
        a1.pos[2] = 0.385523;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 4.367439 + image*box[2];
        a1.pos[2] = 2.207314;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 7.392719 + image*box[2];
        a1.pos[2] = 0.647471;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 8.752916 + image*box[2];
        a1.pos[2] = 0.788170;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 12.863974 + image*box[2];
        a1.pos[2] = 0.359594;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 17.864019 + image*box[2];
        a1.pos[2] = 0.729772;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, 0.0);

        a1.pos[0] = 0.450586 + image*box[2];
        a1.pos[2] = 1.564356;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 2.612290 + image*box[2];
        a1.pos[2] = 0.274607;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 4.022984 + image*box[2];
        a1.pos[2] = 0.859342;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 6.083365 + image*box[2];
        a1.pos[2] = 1.475007;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 9.491282 + image*box[2];
        a1.pos[2] = 0.731473;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);

        a1.pos[0] = 18.987976 + image*box[2];
        a1.pos[2] = 0.490000;
        U = cB.energy(&a1, box);
        EXPECT_EQ (U, -eps);
    }
}

// cylinderZ on its own
class testCylinderZ : public ::testing::Test {
protected:
    atom a1;
    double eps, sigma, radius, width, xc, yc, L;
    int M;
    std::vector < double > box;

    virtual void SetUp() {
		sigma = 1.234;
		radius = 2.0*sigma;
		L = 2*radius+sigma;
		xc = L/2.0; // center of box
		yc = L/2.0; // center of box
        std::vector < double > coords (3, L/2.0);
        a1.pos = coords;
        eps = 2.345;
        width = 1.5*sigma;
        box.resize(3, L);
		M = 1;
    }
};

TEST_F (testCylinderZ, badInit) {
    cylinderZ *cz;

    // safe init
    bool caught = false;
    try {
        cz = new cylinderZ (xc, yc, radius, width, sigma, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (!caught);
	if (!caught) {
		delete cz;
	}

    // bad sigma
    caught = false;
    try {
        cz = new cylinderZ (xc, yc, radius, width, -1.0, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete cz;
	}

    caught = false;
    try {
        cz = new cylinderZ (xc, yc, radius, width, width*2.0, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete cz;
	}

	// bad radius
	caught = false;
    try {
        cz = new cylinderZ (xc, yc, sigma, width, sigma, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete cz;
	}

    // bad M
    caught = false;
    try {
        cz = new cylinderZ (xc, yc, radius, width, sigma, eps, 0);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete cz;
	}

    // bad width
    caught = false;
    try {
        cz = new cylinderZ (xc, yc, radius, -0.001, sigma, eps, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete cz;
	}

    // bad eps
    caught = false;
    try {
        cz = new cylinderZ (xc, yc, radius, width, sigma, -0.001, M);
    } catch (customException &ce) {
        caught = true;
    }
    EXPECT_TRUE (caught);
	if (!caught) {
		delete cz;
	}
}

TEST_F (testCylinderZ, inside) {
 	cylinderZ cz (xc, yc, radius, width, sigma, eps, M);
	bool inside;

	// test valid along entire z-axis
	const double dz = box[2]/10.0;
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		inside = cz.inside(&a1, box);
		EXPECT_TRUE (inside);
	}
	a1.pos[2] = xc;

	// test x
	a1.pos[0] = radius+xc;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = (radius+xc - sigma/2.0);
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = (radius+xc - sigma/2.0)*0.99;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.pos[0] = xc;

	// test y
	a1.pos[1] = radius+yc;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[1] = (radius+yc - sigma/2.0);
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[1] = (radius+yc - sigma/2.0)*0.99;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.pos[1] = yc;

	// test 45 degrees
	a1.pos[0] = xc + (radius - sigma/2.0)/sqrt(2)*1.0001;
	a1.pos[1] = yc + (radius - sigma/2.0)/sqrt(2)*1.0001;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc + (radius - sigma/2.0)/sqrt(2)*0.99;
	a1.pos[1] = yc + (radius - sigma/2.0)/sqrt(2)*0.99;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.pos[0] = xc - (radius - sigma/2.0)/sqrt(2)*1.0001;
	a1.pos[1] = yc - (radius - sigma/2.0)/sqrt(2)*1.0001;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc - (radius - sigma/2.0)/sqrt(2)*0.99;
	a1.pos[1] = yc - (radius - sigma/2.0)/sqrt(2)*0.99;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	// test pbc inside cylinder
	a1.pos[0] = xc + L;
	a1.pos[1] = yc;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.pos[0] = xc + (radius - sigma/2.0)*0.9999 + L;
	a1.pos[1] = yc;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.pos[0] = xc - (radius - sigma/2.0)*0.9999 - L;
	a1.pos[1] = yc;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.pos[0] = xc;
	a1.pos[1] = yc + (radius - sigma/2.0)*0.9999 + L;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.pos[0] = xc;
	a1.pos[1] = yc - (radius - sigma/2.0)*0.9999 - L;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	// test pbc outside of cylinder
	a1.pos[0] = 2*L;
	a1.pos[1] = yc;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc + (radius - sigma/2.0)*1.0001 + L;
	a1.pos[1] = yc;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc - (radius - sigma/2.0)*1.0001 - L;
	a1.pos[1] = yc;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc;
	a1.pos[1] = yc + (radius - sigma/2.0)*1.0001 + L;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc;
	a1.pos[1] = yc - (radius - sigma/2.0)*1.0001 - L;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);
}

TEST_F (testCylinderZ, energy) {
    cylinderZ cz (xc, yc, radius, width, sigma, eps, M);
	double U = 0.0;

	// test U = 0 along entire z-axis at center of cylinder
	const double dz = box[2]/10.0;
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, 0);
	}
	a1.pos[2] = xc;

	// test U = -eps along entire z-axis
	a1.pos[0] = xc + (radius - width*0.99);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, -eps);
	}
	a1.pos[0] = xc;
	a1.pos[2] = xc;

	// test U = 0 along entire z-axis
	a1.pos[0] = xc + (radius - width*1.001);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, 0);
	}
	a1.pos[0] = xc;
	a1.pos[2] = xc;

	// test U = INF along entire z-axis
	a1.pos[0] = xc + (radius - sigma/2.0*0.99);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, NUM_INFINITY);
	}
	a1.pos[0] = xc;
	a1.pos[2] = xc;

	// test U = -eps along entire z-axis
	a1.pos[0] = xc - (radius - width*0.99);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, -eps);
	}
	a1.pos[0] = xc;
	a1.pos[2] = xc;

	// test U = 0 along entire z-axis
	a1.pos[0] = xc - (radius - width*1.001);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, 0);
	}
	a1.pos[0] = xc;
	a1.pos[2] = xc;

	// test U = INF along entire z-axis
	a1.pos[0] = xc - (radius - sigma/2.0*0.99);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, NUM_INFINITY);
	}
	a1.pos[0] = xc;
	a1.pos[2] = xc;

	// test U = -eps along entire z-axis
	a1.pos[1] = yc + (radius - width*0.99);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, -eps);
	}
	a1.pos[1] = yc;
	a1.pos[2] = xc;

	// test U = 0 along entire z-axis
	a1.pos[1] = yc + (radius - width*1.001);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, 0);
	}
	a1.pos[1] = yc;
	a1.pos[2] = xc;

	// test U = INF along entire z-axis
	a1.pos[1] = yc + (radius - sigma/2.0*0.99);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, NUM_INFINITY);
	}
	a1.pos[1] = yc;
	a1.pos[2] = xc;

	// test U = -eps along entire z-axis
	a1.pos[1] = yc - (radius - width*0.99);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, -eps);
	}
	a1.pos[1] = yc;
	a1.pos[2] = xc;

	// test U = 0 along entire z-axis
	a1.pos[1] = yc - (radius - width*1.001);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, 0);
	}
	a1.pos[1] = yc;
	a1.pos[2] = xc;

	// test U = INF along entire z-axis
	a1.pos[1] = yc - (radius - sigma/2.0*0.99);
	for (double z = 0; z <= box[2]; z += dz) {
		a1.pos[2] = z;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, NUM_INFINITY);
	}
	a1.pos[1] = yc;
	a1.pos[2] = xc;
}

TEST_F (testCylinderZ, insideM) {
	int mm = 3;
	bool inside = false;
    cylinderZ cz (xc, yc, radius, width, sigma, eps, mm);
    a1.mState = mm-1;

    // test edge
	a1.pos[0] = xc + (radius - (a1.mState*sigma)/mm/2.0)*1.0001;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc + (radius - (a1.mState*sigma)/mm/2.0)*0.9999;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.mState = mm-2;
	a1.pos[0] = xc + (radius - (a1.mState*sigma)/mm/2.0)*1.0001;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc + (radius - (a1.mState*sigma)/mm/2.0)*0.9999;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);

	a1.mState = 0;
	a1.pos[0] = xc + (radius - sigma/2.0)*1.0001;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (!inside);

	a1.pos[0] = xc + (radius - sigma/2.0)*0.9999;
	inside = cz.inside(&a1, box);
	EXPECT_TRUE (inside);
}

TEST_F (testCylinderZ, energyM) {
	int mm = 3;
	double U = 0.0;
    cylinderZ cz (xc, yc, radius, width, sigma, eps, mm);

	for (int mx = 1; mx <= mm-1; ++mx) {
		a1.mState = mx;

		a1.pos[0] = xc + (radius - (a1.mState*sigma)/mm/2.0)*1.0001;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, NUM_INFINITY);

		a1.pos[0] = xc + (radius - (a1.mState*sigma)/mm/2.0)*0.9999;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, -eps/mm*a1.mState);

		a1.pos[0] = xc;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, 0);

		a1.pos[0] = xc + (radius - width)*0.999;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, 0);

		a1.pos[0] = xc + (radius - width)*1.0001;
		U = cz.energy(&a1, box);
		EXPECT_EQ (U, -eps/mm*a1.mState);
	}

	a1.mState = 0;

	a1.pos[0] = xc + (radius - sigma/2.0)*1.0001;
	U = cz.energy(&a1, box);
	EXPECT_EQ (U, NUM_INFINITY);

	a1.pos[0] = xc + (radius - sigma/2.0)*0.9999;
	U = cz.energy(&a1, box);
	EXPECT_EQ (U, -eps);

	a1.pos[0] = xc;
	U = cz.energy(&a1, box);
	EXPECT_EQ (U, 0);

	a1.pos[0] = xc + (radius - width)*0.999;
	U = cz.energy(&a1, box);
	EXPECT_EQ (U, 0);

	a1.pos[0] = xc + (radius - width)*1.0001;
	U = cz.energy(&a1, box);
	EXPECT_EQ (U, -eps);
}

/* Check potentials */
class quaternionTest : public ::testing::Test {
protected:
	quaternion* q;
	std::vector < double > val;
	double tol;

	virtual void SetUp() {
		q = new quaternion;
		val.resize(4, 0);
		tol = 1.0e-9;
	}
};

TEST_F (quaternionTest, set_get) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;
	q->set(val);

	std::vector < double > is (4);
	is = q->get();

	for (unsigned int i = 0; i < 4; ++i) {
		EXPECT_EQ (is[i], val[i]);
	}
}

TEST_F (quaternionTest, conjugate) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;
	q->set(val);
	q->conjugate();

	std::vector < double > is (4);
	is = q->get();

	EXPECT_EQ (is[0], val[0]);
	for (unsigned int i = 1; i < 4; ++i) {
		EXPECT_EQ (is[i], -val[i]);
	}
}

TEST_F (quaternionTest, getNorm) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;
	std::vector < double > is (4);

	q->set(val);
	is = q->get();

	const double norm = q->getNorm();
	EXPECT_EQ (norm, sqrt(val[0]*val[0]+val[1]*val[1]+val[2]*val[2]+val[3]*val[3]));

	q->conjugate();
	const double norm2 = q->getNorm();
	EXPECT_EQ (norm, norm2);
}

TEST_F (quaternionTest, translate) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;
	std::vector < double > is (4), t (3, -1.0);

	q->set(val);
	q->translate(t);
	is = q->get();

	EXPECT_EQ (is[0], val[0]);
	for (unsigned int i = 1; i < 4; ++i) {
		EXPECT_EQ (is[i], val[i]-1.0);
	}
}

TEST_F (quaternionTest, normalize) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;

	q->set(val);
	q->normalize();
	double new_norm = q->getNorm();

	EXPECT_EQ (new_norm, 1.0);
}

TEST_F (quaternionTest, inverse) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;

	q->set(val);
	q->inverse();
	std::vector< double > is (4);
	is = q->get();

	double n2 = 0.0;
	for (unsigned int i = 0; i < 4; ++i) {
		n2 += val[i]*val[i];
	}
	EXPECT_EQ (is[0], val[0]/n2);
	for (unsigned int i = 1; i < 4; ++i) {
		EXPECT_EQ (is[i], -val[i]/n2);
	}
}

TEST_F (quaternionTest, add1) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;

	q->set(val);

	quaternion p;
	val[0] = -1.234;
	val[1] = -2.345;
	val[2] = -3.456;
	val[3] = -4.567;

	p.set(val);

	quaternion r;
	r = p + (*q);
	std::vector <double> sum = r.get();

	for (unsigned int i = 0; i < 4; ++i) {
		EXPECT_EQ (sum[i], 0.0);
	}
}

TEST_F (quaternionTest, add2) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;

	q->set(val);

	quaternion p;
	p.set(val);

	quaternion r;
	r = p + (*q);
	std::vector <double> sum = r.get();

	for (unsigned int i = 0; i < 4; ++i) {
		EXPECT_EQ (sum[i], val[i]*2);
	}
}

TEST_F (quaternionTest, subtract1) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;

	q->set(val);

	quaternion p;
	val[0] = -1.234;
	val[1] = -2.345;
	val[2] = -3.456;
	val[3] = -4.567;

	p.set(val);

	quaternion r;
	r = p - (*q);
	std::vector <double> sum = r.get();

	for (unsigned int i = 0; i < 4; ++i) {
		EXPECT_EQ (sum[i], 2*val[i]);
	}
}

TEST_F (quaternionTest, subtract2) {
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;

	q->set(val);

	quaternion p;
	p.set(val);

	quaternion r;
	r = p - (*q);
	std::vector <double> sum = r.get();

	for (unsigned int i = 0; i < 4; ++i) {
		EXPECT_EQ (sum[i], 0.0);
	}
}

TEST_F (quaternionTest, quatProd) {
	double w0, w1, x0, x1, y0, y1, z0, z1;
	w0 = 1;
	x0 = 2;
	y0 = 3;
	z0 = 4;
	w1 = 5;
	x1 = 6;
	y1 = 7;
	z1 = 8;

	std::vector < double > v1 (4, 0), v2 (4, 0), res(4, 0);
	v1[0] = w0; v1[1] = x0; v1[2] = y0; v1[3] = z0;
	v2[0] = w1; v2[1] = x1; v2[2] = y1; v2[3] = z1;

	quaternion a, b, c;

	a.set(v1);
	b.set(v2);

	c = a*b;
	res = c.get();

	EXPECT_EQ (res[0], w0*w1-x0*x1-y0*y1-z0*z1);
	EXPECT_EQ (res[1], w0*x1+x0*w1+y0*z1-z0*y1);
	EXPECT_EQ (res[2], w0*y1-x0*z1+y0*w1+z0*x1);
	EXPECT_EQ (res[3], w0*z1+x0*y1-y0*x1+z0*w1);
}

TEST_F (quaternionTest, vecRotatePosX) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	theta = 90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);

	p[1] = 1.0;
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 1.0) < tol);

	// 180
	theta = 180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 270
	theta = 270.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - -1.0) < tol);

	// 360
	theta = 360.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotateNegX) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	theta = -90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);

	p[1] = 1.0;
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - -1.0) < tol);

	// 180
	theta = 180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 270
	theta = -270.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 1.0) < tol);

	// 360
	theta = -360.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotatePosY) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	theta = 90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);

	p[0] = 1.0;
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - -1.0) < tol);

	// 180
	theta = 180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 270
	theta = 270.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 1.0) < tol);

	// 360
	theta = 360.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotateNegY) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	theta = -90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);

	p[0] = 1.0;
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 1.0) < tol);

	// 180
	theta = -180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 270
	theta = -270.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - -1.0) < tol);

	// 360
	theta = -360.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotatePosZ) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	theta = 90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);

	p[0] = 1.0;
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 180
	theta = 180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 270
	theta = 270.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 360
	theta = 360.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotateNegZ) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	theta = -90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);

	p[0] = 1.0;
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 180
	theta = -180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 270
	theta = -270.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);

	// 360
	theta = -360.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotateYZ) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	// (1, 0, 0) about Y then Z
	p[0] = 1.0;

	theta = 180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q->set(a);

	quaternion q2;
	theta = 90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q2.set(a);

	// composite rotation, first q, then q2
	quaternion r = q2*(*q);
	pp = r.rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotateXZ) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	// (0, 1, 0) about X then Z
	p[1] = 1.0;

	theta = 180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);

	quaternion q2;
	theta = 90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q2.set(a);

	// composite rotation, first q, then q2
	quaternion r = q2*(*q);
	pp = r.rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotateXY) {
	std::vector < double > a (4, 0), p (3, 0), pp (3, 0);
	double theta = 0.0, c, s;

	// (0, 0, -1) about X then Y
	p[2] = -1.0;

	theta = 180.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 0.0*s;
	a[3] = 0.0*s;
	q->set(a);

	quaternion q2;
	theta = 90.0*(PI/180.0);
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 1.0*s;
	a[3] = 0.0*s;
	q2.set(a);

	// composite rotation, first q, then q2
	quaternion r = q2*(*q);
	pp = r.rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 0.0) < tol);
}

TEST_F (quaternionTest, vecRotateOn) {
	std::vector < double > a (4, 0), p (3, 3.1415), pp (3, 0);
	double theta = 0.0, c, s;

	// 3.1415*(1, 1, 1) when lies on axis of rotation

	theta = 1.23456*(PI/180.0); // rotation should be irrelevant
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 1.0*s;
	a[2] = 1.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 3.1415) < tol);
	EXPECT_TRUE (fabs(pp[1] - 3.1415) < tol);
	EXPECT_TRUE (fabs(pp[2] - 3.1415) < tol);
}

TEST_F (quaternionTest, vecRotate45) {
	std::vector < double > a (4, 0), p (3, 1), pp (3, 0);
	double theta = 0.0, c, s;

	p[1] = 0;

	theta = 180*(PI/180.0); // rotation should be irrelevant
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 1.0) < tol);

	theta = -90*(PI/180.0); // rotation should be irrelevant
	c = std::cos(theta/2.0); s = std::sin(theta/2.0);
	a[0] = c;
	a[1] = 0.0*s;
	a[2] = 0.0*s;
	a[3] = 1.0*s;
	q->set(a);
	pp = q->rotateVec(p);

	EXPECT_TRUE (fabs(pp[0] - 0.0) < tol);
	EXPECT_TRUE (fabs(pp[1] - -1.0) < tol);
	EXPECT_TRUE (fabs(pp[2] - 1.0) < tol);
}

TEST_F (quaternionTest, axisAngle) {
	quaternion q;
	std::vector < double > axis (3, 1), res (4, 0);
	double angle = 90.0*(PI/180.0);

	q.setAxisAngle (axis, angle);
	res = q.get();

	EXPECT_TRUE (fabs(res[0] - sqrt(2)/2.0) < tol);
	EXPECT_TRUE (fabs(res[1] - sqrt(2)/2.0) < tol);
	EXPECT_TRUE (fabs(res[0] - sqrt(2)/2.0) < tol);
	EXPECT_TRUE (fabs(res[0] - sqrt(2)/2.0) < tol);
}

TEST_F (quaternionTest, randomInit) {
	std::vector < double > last (4, 0);

	// just generate a bunch of trials and make sure each is different
	for (unsigned int i = 0; i < 10; ++i) {
		quaternion q1;
		std::vector < double > v = q1.get();
		for (unsigned int j = 0; j < 4; ++j) {
			EXPECT_TRUE (v[j] != last[j]);
			last[j] = v[j];
		}
	}
}

TEST_F (quaternionTest, equal) {
	quaternion a, b;
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;

	a.set(val);
	b.set(val);
	EXPECT_TRUE (a == b);

	val[3] = 5.678;
	b.set(val);
	EXPECT_TRUE (!(a == b));
}

TEST_F (quaternionTest, notEqual) {
	quaternion a, b;
	val[0] = 1.234;
	val[1] = 2.345;
	val[2] = 3.456;
	val[3] = 4.567;

	a.set(val);
	b.set(val);
	EXPECT_TRUE (!(a != b));

	val[3] = 5.678;
	b.set(val);
	EXPECT_TRUE (a != b);
}

class checkpointTest : public ::testing::Test {
protected:
	moves usedMovesEq, usedMovesPr;
	std::string fname;

	virtual void SetUp() {
		fname = "../data/cpt_input.json";
	}
};

TEST_F (checkpointTest, restarts) {
	// test that system can restart from multiple configurations without error
	simSystem sys = initialize (fname, &usedMovesEq, &usedMovesPr);
	setup (sys, fname);

	bool failed = false;
	try {
		sys.readConfig("../data/res2.xyz");
		sys.readConfig("../data/res3.xyz");
		sys.readConfig("../data/res1.xyz");
	} catch (...) {
		failed = true;
	}
	EXPECT_TRUE (!failed);
}

TEST_F (checkpointTest, load) {
	// test the system can load info from checkpoint
	simSystem sys = initialize (fname, &usedMovesEq, &usedMovesPr);
	setup (sys, fname);
	checkpoint cpt ("../data/checkpt", 900, sys);

	EXPECT_EQ(cpt.hasCheckpoint, true);
	EXPECT_EQ(cpt.takeSnaps, false);
	EXPECT_EQ(cpt.tmmcDone, true);
	EXPECT_EQ(cpt.crossoverDone, true);
	EXPECT_EQ(cpt.walaDone, true);
	EXPECT_EQ(cpt.resFromTMMC, true); // sample was taken from completed sim, but should still flag as in tmmc stage
	EXPECT_EQ(cpt.freq, 900);
	EXPECT_EQ(cpt.moveCounter, 0.0);
	EXPECT_EQ(cpt.sweepCounter, 0.0);
	EXPECT_EQ(cpt.dir, "../data/checkpt");
}

TEST_F (checkpointTest, dump) {
	// test the system can dump
	simSystem sys = initialize (fname, &usedMovesEq, &usedMovesPr);
	setup (sys, fname);
	checkpoint cpt ("../data/checkpt", 900, sys);

	bool failed = false;
	try {
		cpt.dump(sys);
	} catch (...) {
		failed = true;
	}
	EXPECT_TRUE(!failed);
}

TEST_F (checkpointTest, loadDump) {
	// test that system can load,dump same info
	simSystem sys = initialize (fname, &usedMovesEq, &usedMovesPr);
	setup (sys, fname);
	checkpoint cpt ("../data/checkpt", 900, sys);

	// use one system to write a copy of its checkpoint
	cpt.dir = "../data/checkpt2/";
	cpt.chkptName = "../data/checkpt2/state.json"; // dumps a diff checkpoint
	cpt.dump(sys);

	// use a second checkpoint to load
	checkpoint cpt2 ("../data/checkpt2", 900, sys);

	EXPECT_EQ(cpt.hasCheckpoint, cpt2.hasCheckpoint);
	EXPECT_EQ(cpt.takeSnaps, cpt2.takeSnaps);
	EXPECT_EQ(cpt.tmmcDone, cpt2.tmmcDone);
	EXPECT_EQ(cpt.crossoverDone, cpt2.crossoverDone);
	EXPECT_EQ(cpt.walaDone, cpt2.walaDone);
	EXPECT_EQ(cpt.resFromWALA, cpt2.resFromWALA);
	EXPECT_EQ(cpt.resFromCross, cpt2.resFromCross);
	EXPECT_EQ(cpt.resFromTMMC, cpt2.resFromTMMC);
	EXPECT_EQ(cpt.freq, cpt2.freq);
	EXPECT_EQ(cpt.moveCounter, cpt2.moveCounter);
	EXPECT_EQ(cpt.sweepCounter, cpt2.sweepCounter);

	EXPECT_EQ(cpt.dir, cpt2.dir);

	FILE* fp = fopen("../data/checkpt/state.json", "r");
	char readBuffer[65536];
	rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer) );
	rapidjson::Document doc;
	doc.ParseStream(is);
	fclose(fp);

	std::vector < double > ctr1 (doc["extMomCounter"].Size(), 0);
	for (unsigned int i = 0; i < doc["extMomCounter"].Size(); ++i) {
		ctr1[i] = doc["extMomCounter"][i].GetDouble();
	}

	FILE* fp2 = fopen("../data/checkpt2/state.json", "r");
	char readBuffer2[65536];
	rapidjson::FileReadStream is2(fp, readBuffer2, sizeof(readBuffer2) );
	rapidjson::Document doc2;
	doc2.ParseStream(is2);
	fclose(fp2);

	std::vector < double > ctr2 (doc["extMomCounter"].Size(), 0);
	for (unsigned int i = 0; i < doc["extMomCounter"].Size(); ++i) {
		ctr2[i] = doc["extMomCounter"][i].GetDouble();
	}

	EXPECT_EQ (ctr2.size(), ctr1.size());
	const double tol = 1.0e-9;
	for (unsigned int i = 0; i < ctr1.size(); ++i) {
		EXPECT_TRUE(fabs(ctr1[i] - ctr2[i]) < tol);
	}
}

TEST_F (checkpointTest, check) {
	// test the system can check timing correctly
	simSystem sys = initialize (fname, &usedMovesEq, &usedMovesPr);
	setup (sys, fname);
	checkpoint cpt ("../data/checkpt", 2, sys);

	std::cout << "Patience: testing checkpoint timing ..." << std::endl;

	pauseCode(1);
	std::cout << " ... 1/4 ... " << std::endl;
	cpt.dir = "../data/checkpt2/";
	cpt.chkptName = "../data/checkpt2/state.json";
	bool now = cpt.check(sys);
	EXPECT_TRUE (!now);

	pauseCode(1);
	std::cout << " ... 2/4 ... " << std::endl;
	now = cpt.check(sys);
	EXPECT_TRUE(now);

	pauseCode(1);
	std::cout << " ... 3/4 ... " << std::endl;
	now = cpt.check(sys);
	EXPECT_TRUE(!now);

	pauseCode(1);
	std::cout << " ... 4/4 ... " << std::endl;
	now = cpt.check(sys);
	EXPECT_TRUE(now);
}

TEST_F (checkpointTest, restartEnergyHistogram) {
	// test that system can load,dump restartEnergyHistogram info
	simSystem sys = initialize (fname, &usedMovesEq, &usedMovesPr);
	setup (sys, fname);

	sys.restartEnergyHistogram("../data/eHist");
	sys.printEnergyHistogram("eHist_test", false);

	// compare the two files
	int result = system("diff ../data/eHist.dat eHist_test.dat > dmp");
	EXPECT_EQ (result, 0);

	sys.restartEnergyHistogram("../data/eHist");
	sys.printEnergyHistogram("eHist_test", true);

	// normalization changes the result
	result = system("diff ../data/eHist.dat eHist_test.dat > dmp");
	EXPECT_TRUE (result != 0);
}

TEST_F (checkpointTest, restartPkHistogram) {
	// test that system can load,dump restartEnergyHistogram info
	simSystem sys = initialize (fname, &usedMovesEq, &usedMovesPr);
	setup (sys, fname);

	sys.restartPkHistogram("../data/pkHist");
	sys.printPkHistogram("pkHist_test", false);

	// compare the two files
	int result = system("diff ../data/pkHist_1.dat pkHist_test_1.dat > dmp");
	EXPECT_EQ (result, 0);
	result = system("diff ../data/pkHist_2.dat pkHist_test_2.dat > dmp");
	EXPECT_EQ (result, 0);

	sys.restartPkHistogram("../data/pkHist");
	sys.printPkHistogram("pkHist_test", true);

	// normalization changes the result
	result = system("diff ../data/pkHist_1.dat pkHist_test_1.dat > dmp");
	EXPECT_TRUE (result != 0);
	result = system("diff ../data/pkHist_2.dat pkHist_test_2.dat > dmp");
	EXPECT_TRUE (result != 0);
}

TEST_F (checkpointTest, restartExtMoments) {
	// test that system can load,dump restartEnergyHistogram info
	simSystem sys = initialize (fname, &usedMovesEq, &usedMovesPr);
	setup (sys, fname);

	FILE* fp = fopen("../data/checkpt/state.json", "r");
	char readBuffer[65536];
	rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer) );
	rapidjson::Document doc;
	doc.ParseStream(is);
	fclose(fp);

	std::vector < double > ctr (doc["extMomCounter"].Size(), 0);
	for (unsigned int i = 0; i < doc["extMomCounter"].Size(); ++i) {
		ctr[i] = doc["extMomCounter"][i].GetDouble();
	}

	// compare the two files
	sys.restartExtMoments("../data/extMom", ctr);
	sys.printExtMoments("extMom_test", false);
	int result = system("diff ../data/extMom.dat extMom_test.dat > dmp");
	EXPECT_EQ (result, 0);

	// normalization changes result, when ctr != 1 for all members of ctr
	sys.restartExtMoments("../data/extMom", ctr);
	sys.printExtMoments("extMom_test", true);
	result = system("diff ../data/extMom.dat extMom_test.dat > dmp");
	EXPECT_TRUE (result != 0);

	// when ctr = 1 for all, normalization should not not changes
	std::fill(ctr.begin(), ctr.end(), 1);
	sys.restartExtMoments("../data/extMom", ctr);
	sys.printExtMoments("extMom_test", false);
	result = system("diff ../data/extMom.dat extMom_test.dat > dmp");
	EXPECT_EQ (result, 0);
}

TEST (testNone, cleanUp) {
	// cleanup after other tests
	int res;
	res = system("rm ppot_*");
	EXPECT_EQ (res, 0);
	res = system("rm tmmBias*");
	EXPECT_EQ (res, 0);
	res = system("rm walaBias*");
	EXPECT_EQ (res, 0);
	res = system("rm dmp eHist_test.dat");
	EXPECT_EQ (res, 0);
	res = system("rm ../data/checkpt2/e* ../data/checkpt2/pk* ../data/checkpt2/snap.xyz ../data/checkpt2/tmmc* ");
	EXPECT_EQ (res, 0);
	res = system("rm ../data/checkpt/e* ../data/checkpt/pk* ../data/checkpt/tmmc_lnPI.dat ");
	EXPECT_EQ (res, 0);
	res = system("rm pkHist_test_1.dat pkHist_test_2.dat");
	EXPECT_EQ (res, 0);
	res = system("rm extMom_test.dat");
	EXPECT_EQ (res, 0);
}
