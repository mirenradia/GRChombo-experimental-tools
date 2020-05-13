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
#include "GRChomboBCs.hpp"
#include "ParmParse.H"
#include "SetBCs.H"
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
    int gridBufferSize;
    int tagBufferSize;
    Real fillRatio;
    Real refineThresh;
    int num_bh_tagging_levels;
    Real horizon_tagging_buffer;
    bool forceGRChomboRegrid;
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
    int num_ghosts;
    SetBCs::params_t bc_params;
    GRChomboBCs::params_t grchombo_boundary_params;
    bool symmetric_boundaries_exist;
    bool nonperiodic_boundaries_exist;
    Real pre_collapsed_lapse_exponent;
    Real alpha;
    Real beta;
    Real G_Newton;
    Real bh_bare_mass;
    RealVect bh_spin;
    RealVect bh_momentum;
    RealVect bh_offset;
    BosonStar_params_t boson_star_params;
    bool rescale_radius;
    Real binary_separation;
    Potential::params_t potential_params;
    int extraction_level; // will ensure this level contains the extraction
                          // radius and a 20% buffer zone
    Real extraction_radius;
};

/// get stuff from input file
extern void getPoissonParameters(PoissonParameters &a_params);

#endif
