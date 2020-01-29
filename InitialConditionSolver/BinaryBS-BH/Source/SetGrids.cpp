#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SetGrids.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "PoissonParameters.H"
#include "SetLevelData.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Setup grids
// Note that there is also an option to read in grids, but here we use tagging
// for the refinement
int set_grids(Vector<DisjointBoxLayout> &vectGrids,
                     PoissonParameters &a_params,
                     BosonStar &a_boson_star)
{
    Vector<ProblemDomain> vectDomain;
    Vector<Real> vectDx;
    set_domains_and_dx(vectDomain, vectDx, a_params);

    int numlevels = a_params.numLevels;

    ParmParse pp;

    // grid generation parameters
    vectGrids.resize(numlevels);

    int maxLevel = numlevels - 1;
    Vector<Vector<Box>> newBoxes(numlevels);
    Vector<Vector<Box>> oldBoxes(numlevels);

    // determine grids dynamically, based on grad(RHS)
    // will need temp storage for RHS
    Vector<LevelData<FArrayBox> *> vectRHS(maxLevel + 1, NULL);

    // define base level first
    Vector<Vector<int>> procAssign(maxLevel + 1);
    domainSplit(vectDomain[0], oldBoxes[0], a_params.maxGridSize,
                a_params.blockFactor);
    procAssign[0].resize(oldBoxes[0].size());
    LoadBalance(procAssign[0], oldBoxes[0]);
    vectGrids[0].define(oldBoxes[0], procAssign[0], vectDomain[0]);
    vectRHS[0] = new LevelData<FArrayBox>(vectGrids[0], 1, IntVect::Zero);

    int topLevel = 0;
    bool moreLevels = (maxLevel > 0);

    int nesting_radius = 2;
    // create grid generation object
    BRMeshRefine meshrefine(vectDomain[0], a_params.refRatio,
                            a_params.fillRatio, a_params.blockFactor,
                            nesting_radius, a_params.maxGridSize);

    while (moreLevels)
    {
        // default is moreLevels = false
        // (only repeat loop in the case where a new level
        // is generated which is still less than maxLevel)
        moreLevels = false;

        int baseLevel = 0;
        int oldTopLevel = topLevel;

        // now initialize RHS for this existing hierarchy
        for (int level = 0; level <= topLevel; level++)
        {
            RealVect dxLevel = vectDx[level] * RealVect::Unit;

            LevelData<FArrayBox> *temp_multigrid_vars;
            LevelData<FArrayBox> *temp_dpsi;

            temp_multigrid_vars = new LevelData<FArrayBox>(
                vectGrids[level], NUM_MULTIGRID_VARS, 3 * IntVect::Unit);
            temp_dpsi = new LevelData<FArrayBox>(vectGrids[level], 1,
                                                 3 * IntVect::Unit);

            GRChomboBCs grchombo_boundaries;
            grchombo_boundaries.define(dxLevel[0],
                                       a_params.grchombo_boundary_params,
                                       a_params.coarsestDomain,
                                       a_params.num_ghosts);

            set_initial_conditions(*temp_multigrid_vars, *temp_dpsi,
                                   grchombo_boundaries, dxLevel, a_params,
                                   a_boson_star);

            // set condition for regrid - use the integrability condition
            // integral
            set_regrid_condition(*vectRHS[level], *temp_multigrid_vars, dxLevel,
                                 a_params);

            if (temp_multigrid_vars != NULL)
            {
                delete temp_multigrid_vars;
                temp_multigrid_vars = NULL;
            }
            if (temp_dpsi != NULL)
            {
                delete temp_dpsi;
                temp_dpsi = NULL;
            }
        }

        Vector<IntVectSet> tagVect(topLevel + 1);
        int tags_grow = a_params.bufferSize;
        set_tag_cells(vectRHS, tagVect, vectDx, vectDomain, baseLevel,
                      topLevel + 1, a_params);

        int new_finest =
            meshrefine.regrid(newBoxes, tagVect, baseLevel, topLevel, oldBoxes);

        if (new_finest > topLevel)
        {
            topLevel++;
        }

        oldBoxes = newBoxes;

        //  no need to do this for the base level (already done)
        for (int lev = 1; lev <= topLevel; lev++)
        {
            // do load balancing
            procAssign[lev].resize(newBoxes[lev].size());
            LoadBalance(procAssign[lev], newBoxes[lev]);
            const DisjointBoxLayout newDBL(newBoxes[lev], procAssign[lev],
                                           vectDomain[lev]);
            vectGrids[lev] = newDBL;
            delete vectRHS[lev];
            vectRHS[lev] =
                new LevelData<FArrayBox>(vectGrids[lev], 1, IntVect::Zero);
        } // end loop over levels for initialization

        // figure out whether we need another pass through grid generation
        if ((topLevel < maxLevel) && (topLevel > oldTopLevel))
            moreLevels = true;

    } // end while moreLevels loop

    // clean up temp storage
    for (int ilev = 0; ilev < vectRHS.size(); ilev++)
    {
        if (vectRHS[ilev] != NULL)
        {
            delete vectRHS[ilev];
            vectRHS[ilev] = NULL;
        }
    }

    pout() << "\n----------------------------------------------\n"
           << "\x1b[1mLoad Balancing information:\x1b[0m\n";
    // Get number of boxes on this rank and level and print it
    for (int ilev = 0; ilev < numlevels; ++ilev)
    {
        int nbox = vectGrids[ilev].dataIterator().size();
        pout() << "Number of boxes on level " << ilev << " on this rank: "
               << nbox << "\n";
    }
    pout() << "----------------------------------------------\n" << std::endl;
    return 0;
}

// Set grid hierarchy from input file
void set_domains_and_dx(Vector<ProblemDomain> &vectDomain, Vector<Real> &vectDx,
                        PoissonParameters &a_params)
{

    vectDomain.resize(a_params.numLevels);
    vectDx.resize(a_params.numLevels);
    vectDx[0] = a_params.coarsestDx;
    for (int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
        vectDx[ilev] = vectDx[ilev - 1] / a_params.refRatio[ilev - 1];
    }

    vectDomain[0] = a_params.coarsestDomain;
    for (int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
        vectDomain[ilev] =
            refine(vectDomain[ilev - 1], a_params.refRatio[ilev - 1]);
    }
}

/*
  tag cells for refinement based on magnitude(RHS)
*/
void set_tag_cells(Vector<LevelData<FArrayBox> *> &vectRHS,
                   Vector<IntVectSet> &tagVect, Vector<Real> &vectDx,
                   Vector<ProblemDomain> &vectDomain, const int baseLevel,
                   int numLevels, const PoissonParameters &a_params)
{
    for (int lev = baseLevel; lev != numLevels; lev++)
    {
        IntVectSet local_tags;
        LevelData<FArrayBox> &levelRhs = *vectRHS[lev];
        DisjointBoxLayout level_domain = levelRhs.getBoxes();
        DataIterator dit = levelRhs.dataIterator();

        Real maxRHS = 0;

        maxRHS = norm(levelRhs, levelRhs.interval(), 0);

        Real tagVal = maxRHS * a_params.refineThresh;

        // now loop through grids and tag cells where RHS > tagVal
        for (dit.reset(); dit.ok(); ++dit)
        {
            const Box thisBox = level_domain.get(dit());
            const FArrayBox &thisRhs = levelRhs[dit()];
            BoxIterator bit(thisBox);
            for (bit.begin(); bit.ok(); ++bit)
            {
                const IntVect &iv = bit();
                bool extraction_tag = false;
                if (lev < a_params.extraction_level)
                {
                    RealVect loc;
                    RealVect dx_vect = vectDx[lev] * RealVect::Unit;
                    get_loc(loc, iv, dx_vect, a_params);
                    Real r = loc.vectorLength();
                    extraction_tag = (r < 1.2 * a_params.extraction_radius);
                }
                if (abs(thisRhs(iv)) >= tagVal || extraction_tag)
                    local_tags |= iv;
            }
        } // end loop over grids on this level

        local_tags.grow(a_params.bufferSize);
        const Box &domainBox = vectDomain[lev].domainBox();
        local_tags &= domainBox;

        tagVect[lev] = local_tags;

    } // end loop over levels
}
