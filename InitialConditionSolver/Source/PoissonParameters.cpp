#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GRParmParse.hpp"
#include "PoissonParameters.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// function to read in the key params for solver
void getPoissonParameters(PoissonParameters &a_params)
{
    GRParmParse pp;

    // Set verbosity
    //a_params.verbosity = 3;
    pp.load("verbosity", a_params.verbosity);

    // Chombo grid params
    pp.load("max_level", a_params.maxLevel);
    a_params.numLevels = a_params.maxLevel + 1;
    std::vector<int> nCellsArray(SpaceDim);
    pp.load("N", nCellsArray, SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        a_params.nCells[idir] = nCellsArray[idir];
    }

    // Enforce that dx is same in every directions
    // and that ref_ratio = 2 always as these conditions
    // are required in several places in our code
    a_params.refRatio.resize(a_params.numLevels);
    a_params.refRatio.assign(2);
    Real domain_length;
    pp.load("L", domain_length);
    a_params.coarsestDx = domain_length / a_params.nCells[0];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        a_params.domainLength[idir] =
            a_params.coarsestDx * a_params.nCells[idir];
    }

    // Chombo refinement and load balancing criteria
    pp.load("refine_threshold", a_params.refineThresh);

    // alias the weird chombo names to something more descriptive
    // for these box params, and default to some reasonable values
    if (pp.contains("block_factor"))
    {
        pp.load("block_factor", a_params.blockFactor);
    }
    else
    {
        pp.load("min_box_size", a_params.blockFactor, 8);
    }
    if (pp.contains("max_grid_size"))
    {
        pp.load("max_grid_size", a_params.maxGridSize);
    }
    else
    {
        pp.load("max_box_size", a_params.maxGridSize, 64);
    }
    pp.load("fill_ratio", a_params.fillRatio);
    pp.load("buffer_size", a_params.bufferSize);

    // set average type -
    // set to a bogus default value, so we only break from solver
    // default if it's set to something real
    a_params.coefficient_average_type = -1;
    if (pp.contains("coefficient_average_type"))
    {
        std::string tempString;
        pp.load("coefficient_average_type", tempString);
        if (tempString == "arithmetic")
        {
            a_params.coefficient_average_type = CoarseAverage::arithmetic;
        }
        else if (tempString == "harmonic")
        {
            a_params.coefficient_average_type = CoarseAverage::harmonic;
        }
        else
        {
            MayDay::Error("bad coefficient_average_type in input");
        }
    } // end if an average_type is present in inputs

    // set up coarse domain box
    IntVect lo = IntVect::Zero;
    IntVect hi = a_params.nCells;
    hi -= IntVect::Unit;
    Box crseDomBox(lo, hi);
    a_params.probLo = RealVect::Zero;
    a_params.probHi = RealVect::Zero;
    a_params.probHi += a_params.domainLength;

    // Periodicity - for the moment enforce same in all directions
    ProblemDomain crseDom(crseDomBox);
    int is_periodic;
    pp.load("is_periodic", is_periodic, 0);
    a_params.periodic.resize(SpaceDim);
    a_params.periodic.assign(is_periodic);
    for (int dir = 0; dir < SpaceDim; dir++)
    {
        crseDom.setPeriodic(dir, is_periodic);
    }
    a_params.coarsestDomain = crseDom;

    pout() << "periodicity = " << is_periodic << endl;

    // problem specific params
    //pp.get("alpha", a_params.alpha);
    //pp.get("beta", a_params.beta);
    // hardcode these to the expected values
    a_params.alpha = 1.0;
    a_params.beta = -1.0
    // print out the overall coeffs just to be sure we have selected them
    // correctly
    pout() << "alpha, beta = " << a_params.alpha << ", " << a_params.beta
           << endl;

    // Initial conditions for the scalar field
    pp.load("G_Newton", a_params.G_Newton);

    // Common Boson Star parameters
    pp.load("abs_error", a_params.boson_star1_params.abs_error, 1.0e-14);
    pp.load("rel_error", a_params.boson_star1_params.rel_error, 1.0e-14);
    pp.load("initial_step_size", a_params.boson_star1_params.initial_step_size,
            0.015625);
    pp.load("max_radius", a_params.boson_star1_params.max_radius, 256.);
    pp.load("binary_search_tol", a_params.boson_star1_params.binary_search_tol,
            1.0e-15);
    pp.load("max_binary_search_iter",
            a_params.boson_star1_params.max_binary_search_iter, 1000);

    if (a_params.verbosity)
    {
        pout() << "\nBoson Star Profile Solver Parameters:\n"
               << "abs_error = "
                    << a_params.boson_star1_params.abs_error << "\n"
               << "rel_error = "
                    << a_params.boson_star1_params.rel_error << "\n"
               << "initial_step_size = "
                    << a_params.boson_star1_params.initial_step_size << "\n"
               << "max_radius = "
                    << a_params.boson_star1_params.max_radius << "\n"
               << "binary_search_tol = "
                    << a_params.boson_star1_params.binary_search_tol << "\n"
               << "max_binary_search_iter = "
                    << a_params.boson_star1_params.max_binary_search_iter
               << endl;
    }

    // Boson Star 1 parameters
    pp.load("central_amplitude_CSF1",
            boson_star1_params.central_amplitude_CSF);
    pp.load("phase1", boson_star1_params.phase, 0.0);
    pp.load("star_centre1", boson_star1_params.star_centre);

    pout() << "\nBoson Star 1 parameters:\n"
           << "central_amplitude = "
                << a_params.boson_star1_params.central_amplitude_CSF << "\n"
           << "phase = " << a_params.boson_star1_params.phase << "\n"
           << "centre = (" << a_params.boson_star1_params.centre[0] << ", "
                           << a_params.boson_star1_params.centre[1] << ", "
                           << a_params.boson_star2_params.centre[2] << ")"
           << endl;


    // Initialize values for boson_star2_params to same as boson_star1_params
    // and then assign that ones that should differ below
    boson_star2_params = boson_star1_params;

    // Are the two stars' profiles identical
    pp.load("identical", identical, false);

    // Boson Star 2 parameters
    if (!identical)
    {
        pp.load("central_amplitude_CSF2",
                boson_star2_params.central_amplitude_CSF);
    }
    pp.load("phase2", boson_star2_params.phase, 0.0);
    pp.load("star_centre2", boson_star2_params.star_centre,
            {L - boson_star1_params.star_centre[0],
             L - boson_star1_params.star_centre[1],
             L - boson_star1_params.star_centre[2]});

    pout() << "\nBoson Star 2 parameters:\n"
           << "central_amplitude = "
                << a_params.boson_star2_params.central_amplitude_CSF << "\n"
           << "phase = " << a_params.boson_star2_params.phase << "\n"
           << "centre = (" << a_params.boson_star2_params.centre[0] << ", "
                           << a_params.boson_star2_params.centre[1] << ", "
                           << a_params.boson_star2_params.centre[2] << ")"
           << endl;

    // Potential params
    pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
    pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);

    pout() << "\nPotential parameters:\n"
           << "scalar_mass = " << potential_params.scalar_mass << "\n"
           << "phi^4_coeff = " << potential_params.phi4_coeff << endl;
}
