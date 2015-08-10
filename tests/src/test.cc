#include "gtest-1.7.0/include/gtest/gtest.h"
#include <vector>
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
#include <iostream>

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
		params.resize(4, 0);
		tol = 1.0e-9;
	}
};

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
	
	EXPECT_TRUE ( fabs(lj->energy(pow(2.0, 1./6.)) - -eps) < tol );
	EXPECT_TRUE ( fabs(lj->energy(1.0001*rcut) - 0) < tol );
	EXPECT_TRUE ( lj->energy(0.9999*rcut) < 0.0 );
	EXPECT_TRUE ( fabs(lj->energy(sigma) - 0) < tol );
	EXPECT_TRUE ( fabs(lj->energy(1.234) - -0.812008901) < tol );
	
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
	
	EXPECT_TRUE ( fabs(lj->energy(pow(2.0, 1./6.)) - 0.0) < tol );
	EXPECT_TRUE ( fabs(lj->energy(1.0001*rcut) - 0) < tol );
	EXPECT_TRUE ( lj->energy(0.9999*rcut) > 0.0 );
	EXPECT_TRUE ( fabs(lj->energy(0.987) - 1.353386603) < tol );
}

/* Check the System */
class InitializeSystem : public ::testing::Test {
protected:
	unsigned int nSpecies;
	double beta, tol;
	std::vector < double > box, mu;
	std::vector < int > maxSpecies, minSpecies;
	
	virtual void SetUp() {
		nSpecies = 3;
		maxSpecies.resize(nSpecies), mu.resize(nSpecies), box.resize(3), minSpecies.resize(nSpecies);
		tol = 1.0e-9;
		beta = 1.234;
		for (unsigned int i = 0; i < nSpecies; ++i) {
			box[i] = 2*i+1.0;
			mu[i] = 2*i+2.0;
			maxSpecies[i] = 10*i+10;
			minSpecies[i] = 0;
		}
	}
};

TEST_F (InitializeSystem, nSpecies) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	EXPECT_EQ (mysys.nSpecies(), nSpecies);
}

TEST_F (InitializeSystem, beta) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	EXPECT_EQ (mysys.beta(), beta);
}

TEST_F (InitializeSystem, box) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	EXPECT_EQ (mysys.box(), box);
}

TEST_F (InitializeSystem, mu) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	for (unsigned int i = 0; i < mysys.nSpecies(); ++i) {
		EXPECT_EQ (mysys.mu(i), mu[i]);
	}
}

TEST_F (InitializeSystem, maxSpecies) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	for (unsigned int i = 0; i < mysys.nSpecies(); ++i) {
		EXPECT_EQ (mysys.maxSpecies(i), maxSpecies[i]);
	}
}

TEST_F (InitializeSystem, minSpecies) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	for (unsigned int i = 0; i < mysys.nSpecies(); ++i) {
		EXPECT_EQ (mysys.minSpecies(i), minSpecies[i]);
	}
}

TEST_F (InitializeSystem, incrementEnergy) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	double dU = -0.1;
	mysys.incrementEnergy(dU);
	EXPECT_TRUE (fabs(mysys.energy() - dU) < tol);
	mysys.incrementEnergy(-dU);
	EXPECT_TRUE (fabs(mysys.energy() - 0) < tol);
}

TEST_F (InitializeSystem, addPotential) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	lennardJones ljtest;
	mysys.addPotential(0, 0, &ljtest);
	EXPECT_EQ (mysys.potentialIsSet (0, 0), true);
	EXPECT_EQ (mysys.potentialIsSet (0, 1), false);
	EXPECT_EQ (mysys.potentialIsSet (1, 1), false);
	EXPECT_EQ (mysys.potentialIsSet (1, 1), false);
}

TEST_F (InitializeSystem, insertAtom) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	atom typeA1, typeA2, typeB1, typeB2;
	std::vector < double > pos (3, 1.234), pos2 (3, 2.345);
	typeA1.pos = pos;
	typeA2.pos = pos2;
	
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
}

TEST_F (InitializeSystem, deleteAtom) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	atom typeA1, typeA2, typeB1, typeB2;
	std::vector < double > pos (3, 1.234), pos2 (3, 2.345);
	typeA1.pos = pos;
	typeA2.pos = pos2;
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
}

TEST_F (InitializeSystem, scratchEnergy) {
	simSystem mysys (nSpecies, beta, box, mu, maxSpecies, minSpecies);
	lennardJones ljtest;
	std::vector < double > params (4), p1(3, 0.0), p2(3, 0.0);
	
	// wca
	params[0] = 1.0;
	params[1] = 1.0;
	params[2] = pow(2.0, 1./6.);
	params[3] = 1.0;
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