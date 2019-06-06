#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POISSONPARAMETERS_H_
#define _POISSONPARAMETERS_H_

#include <iostream>
using std::cerr;

#include "AMRMultiGrid.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "VariableCoeffPoissonOperator.H"

// Boson Star includes
#include "BosonStarParams.hpp"
#include "ComplexPotential.hpp"

/// Class for solver input parameters
/**
 */
class PoissonParameters
{
  public:
    IntVect nCells;
    int maxGridSize;
    int blockFactor;
    int bufferSize;
    Real fillRatio;
    Real refineThresh;
    int coefficient_average_type;
    int verbosity;
    Vector<int> periodic;
    int domBcType;
    int maxLevel;
    int numLevels;
    Vector<int> refRatio;
    ProblemDomain coarsestDomain;
    Real coarsestDx;
    RealVect domainLength;
    RealVect probLo;
    RealVect probHi;
    Real alpha;
    Real beta;
    Real G_Newton;
    BosonStar_params_t boson_star1_params, boson_star2_params;
    bool identical; // whether or not the 2 boson stars have the same profile
    Potential::params_t potential_params;
};

/// get stuff from input file
extern void getPoissonParameters(PoissonParameters &a_params);

#endif
