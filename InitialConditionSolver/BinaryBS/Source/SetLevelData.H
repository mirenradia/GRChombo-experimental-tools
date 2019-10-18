/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _SETLEVELDATA_H_
#define _SETLEVELDATA_H_

#include <iostream>
using std::cerr;

#include "AMRMultiGrid.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "GRChomboUserVariables.hpp"
#include "MultigridUserVariables.hpp"
#include "GRChomboBCs.hpp"
#include "ParmParse.H"
#include "PoissonParameters.H"
#include "VariableCoeffPoissonOperator.H"
#include "VariableCoeffPoissonOperatorFactory.H"

// Boson Star includes
#include "BosonStar.ICS.hpp"
#include "ComplexPotential.hpp"

// various functions operating on grids to set up data for initial
// conditions, spatially dependent coefficients and rhs etc

// This takes an IntVect and writes the physical coordinates to a RealVect
extern inline void get_loc(RealVect &a_out_loc,
                           const IntVect &a_iv,
                           const RealVect &a_dx,
                           const PoissonParameters &a_params)
{
    a_out_loc = a_iv + 0.5 * RealVect::Unit;
    a_out_loc *= a_dx;
    a_out_loc -= a_params.domainLength / 2.0;
}


extern void set_initial_conditions(LevelData<FArrayBox> &a_multigrid_vars,
                                   LevelData<FArrayBox> &a_dpsi,
                                   GRChomboBCs &a_grchombo_boundaries,
                                   const RealVect &a_dx,
                                   const PoissonParameters &a_params,
                                   BosonStar &a_boson_star1,
                                   BosonStar &a_boson_star2);

inline void set_initial_multigrid_cell(FArrayBox &a_multigrid_vars_box,
                                FArrayBox &a_dpsi_box,
                                const IntVect &a_iv,
                                const RealVect &a_dx,
                                const PoissonParameters &a_params,
                                BosonStar &a_boson_star1,
                                BosonStar &a_boson_star2);

extern void set_rhs(LevelData<FArrayBox> &a_rhs,
                    LevelData<FArrayBox> &a_multigrid_vars,
                    const RealVect &a_dx, const PoissonParameters &a_params);

extern inline void set_laplacian_psi(Real &laplacian_of_psi_0,
                                     const IntVect &a_iv,
                                     const FArrayBox &a_multigrid_vars_box,
                                     const RealVect &a_dx);

extern inline void set_grad_phi_sq(Real &grad_phi_sq,
                                   const IntVect &a_iv,
                                   const FArrayBox &a_multigrid_vars_box,
                                   const RealVect &a_dx);

extern void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                                 LevelData<FArrayBox> &a_multigrid_vars,
                                 const RealVect &a_dx,
                                 const PoissonParameters &a_params);

extern void set_a_coef(LevelData<FArrayBox> &a_aCoef,
                       LevelData<FArrayBox> &a_multigrid_vars,
                       const PoissonParameters &a_params, const RealVect &a_dx);

extern void set_b_coef(LevelData<FArrayBox> &a_bCoef,
                       const PoissonParameters &a_params, const RealVect &a_dx);

extern void set_update_psi0(LevelData<FArrayBox> &a_multigrid_vars,
                            LevelData<FArrayBox> &a_dpsi,
                            const Copier &a_exchange_copier);

extern void set_output_data(LevelData<FArrayBox> &a_vars,
                            LevelData<FArrayBox> &a_multigrid_vars,
                            GRChomboBCs &a_grchombo_boundaries,
                            const PoissonParameters &a_params,
                            const RealVect &a_dx);

inline void set_non_const_output_cell(const FArrayBox &a_multigrid_vars_box,
                                      FArrayBox &a_grchombo_vars_box,
                                      const IntVect &a_iv);

#endif
