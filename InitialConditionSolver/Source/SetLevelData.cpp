/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SetLevelData.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoHiSide.H"
#include "LoadBalance.H"
#include "PoissonParameters.H"
#include "SetBinaryBH.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

// This takes an IntVect and writes the physical coordinates to a RealVect
inline void get_loc(RealVect &a_out_loc, const IntVect &a_iv,
                    const RealVect &a_dx, const PoissonParameters &a_params)
{
    a_out_loc = a_iv + 0.5 * RealVect::Unit;
    a_out_loc *= a_dx;
    a_out_loc -= a_params.domainLength / 2.0;
}

// set initial guess value for the conformal factor psi
// defined by \gamma_ij = \psi^4 \tilde \gamma_ij, and \bar Aij = psi^2 A_ij.
// For now the default setup is 2 Bowen York BHs
void set_initial_conditions(LevelData<FArrayBox> &a_multigrid_vars,
                            LevelData<FArrayBox> &a_dpsi,
                            GRChomboBCs &a_grchombo_boundaries,
                            const RealVect &a_dx,
                            const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &dpsi_box = a_dpsi[dit()];
        Box this_box = multigrid_vars_box.box();
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            set_initial_multigrid_cell(multigrid_vars_box, dpsi_box, bit(),
                                       a_dx, a_params);
        }

        // now fill boundary ghost cells if using nonperiodic boundaries in
        // GRChombo. Note that these cells are unused in the
        IntVect offset_lo, offset_hi;
        a_grchombo_boundaries.get_box_offsets(offset_lo, offset_hi, this_box);

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        a_grchombo_boundaries.remove_outer_ghost_cells(this_box);

        for (int idir = 0; idir < SpaceDim; ++idir)
        {
            if (!a_params.periodic[idir])
            {
                for (SideIterator sit; sit.ok(); ++sit)
                {
                    Box boundary_box = a_grchombo_boundaries.get_boundary_box(
                        sit(), idir, offset_lo, offset_hi, this_box);

                    // now we have the appropriate box, fill it!
                    BoxIterator bbit(boundary_box);
                    for (bbit.begin(); bbit.ok(); ++bbit)
                    {
                        set_initial_multigrid_cell(multigrid_vars_box, dpsi_box,
                                                   bbit(), a_dx, a_params);
                    } // end loop through boundary box
                }     // end loop over sides
            }         // end if (periodic[idir])
        }             // end loop over directions
    }
} // end set_initial_conditions

void set_initial_multigrid_cell(FArrayBox &a_multigrid_vars_box,
                                FArrayBox &a_dpsi_box, const IntVect &a_iv,
                                const RealVect &a_dx,
                                const PoissonParameters &a_params)
{
    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    // set psi to 1.0 and zero dpsi
    // note that we don't include the singular part of psi
    // for the BHs - this is added at the output data stage
    // and when we calculate psi_0 in the rhs etc
    // as it already satisfies Laplacian(psi) = 0
    a_multigrid_vars_box(a_iv, c_psi_reg) = 1.0;
    a_dpsi_box(a_iv, 0) = 0.0;

    // set Aij for spin and momentum according to BH params
    set_binary_bh_Aij(a_multigrid_vars_box, a_iv, loc, a_params);
} // end set set_initial_multigrid_cell

// computes the Laplacian of psi at a point in a box
inline Real get_laplacian_psi(const IntVect &a_iv, const FArrayBox &a_psi_fab,
                              const RealVect &a_dx)
{
    Real laplacian_of_psi = 0.0;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
        IntVect iv_offset1 = a_iv;
        IntVect iv_offset2 = a_iv;
        iv_offset1[idir] -= 1;
        iv_offset2[idir] += 1;

        // 2nd order stencil for now
        Real d2psi_dxdx = 1.0 / (a_dx[idir] * a_dx[idir]) *
                          (+1.0 * a_psi_fab(iv_offset2) -
                           2.0 * a_psi_fab(a_iv) + 1.0 * a_psi_fab(iv_offset1));
        laplacian_of_psi += d2psi_dxdx;
    }
    return laplacian_of_psi;
} // end get_laplacian_psi

// set the rhs source for the poisson eqn
void set_rhs(LevelData<FArrayBox> &a_rhs,
             LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx,
             const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_rhs.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &rhs_box = a_rhs[dit()];
        rhs_box.setVal(0.0, 0);
        Box this_box = rhs_box.box(); // no ghost cells

        FArrayBox psi_fab(Interval(c_psi_reg, c_psi_reg), multigrid_vars_box);

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            Real laplacian_of_psi = get_laplacian_psi(iv, psi_fab, a_dx);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            Real psi_bl = get_psi_brill_lindquist(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bl;

            rhs_box(iv, 0) = -0.125 * A2 * pow(psi_0, -7.0) - laplacian_of_psi;
        }
    }
} // end set_rhs

// set the regrid condition - abs value of this drives AMR regrid
void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                          LevelData<FArrayBox> &a_multigrid_vars,
                          const RealVect &a_dx,
                          const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_condition.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &condition_box = a_condition[dit()];
        condition_box.setVal(0.0, 0);
        Box this_box = condition_box.box(); // no ghost cells

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            // the condition is similar to the rhs but we take abs
            // value of the contributions and add in BH criteria
            Real psi_bl = get_psi_brill_lindquist(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bl;
            condition_box(iv, 0) = 1.5 * A2 * pow(psi_0, -7.0) + log(psi_0);
        }
    }
} // end set_regrid_condition

// Add the correction to psi0 after the solver operates
void set_update_psi0(LevelData<FArrayBox> &a_multigrid_vars,
                     LevelData<FArrayBox> &a_dpsi,
                     const Copier &a_exchange_copier)
{

    // first exchange ghost cells for dpsi so they are filled with the correct
    // values
    a_dpsi.exchange(a_dpsi.interval(), a_exchange_copier);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &dpsi_box = a_dpsi[dit()];

        Box this_box = multigrid_vars_box.box();
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            multigrid_vars_box(iv, c_psi_reg) += dpsi_box(iv, 0);
        }
    }
}

// The coefficient of the I operator on dpsi
void set_a_coef(LevelData<FArrayBox> &a_aCoef,
                LevelData<FArrayBox> &a_multigrid_vars,
                const PoissonParameters &a_params, const RealVect &a_dx)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_aCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &aCoef_box = a_aCoef[dit()];
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        Box this_box = aCoef_box.box();

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            Real psi_bl = get_psi_brill_lindquist(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bl;
            aCoef_box(iv, 0) = -0.875 * A2 * pow(psi_0, -8.0);
        }
    }
}

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
// the rhs source of the Poisson eqn
void set_b_coef(LevelData<FArrayBox> &a_bCoef,
                const PoissonParameters &a_params, const RealVect &a_dx)
{

    CH_assert(a_bCoef.nComp() == 1);
    int comp_number = 0;

    for (DataIterator dit = a_bCoef.dataIterator(); dit.ok(); ++dit)
    {
        FArrayBox &bCoef_box = a_bCoef[dit()];
        bCoef_box.setVal(1.0, comp_number);
    }
}

// used to set output data for all ADM Vars for GRChombo restart
void set_output_data(LevelData<FArrayBox> &a_grchombo_vars,
                     LevelData<FArrayBox> &a_multigrid_vars,
                     GRChomboBCs &a_grchombo_boundaries,
                     const PoissonParameters &a_params, const RealVect &a_dx)
{
    CH_assert(a_grchombo_vars.nComp() == NUM_GRCHOMBO_VARS);
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_grchombo_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &grchombo_vars_box = a_grchombo_vars[dit()];
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        // first set everything to zero
        for (int comp = 0; comp < NUM_GRCHOMBO_VARS; comp++)
        {
            grchombo_vars_box.setVal(0.0, comp);
        }

        // now set non zero terms - const across whole box
        // Conformally flat, and lapse = 1
        grchombo_vars_box.setVal(1.0, c_h11);
        grchombo_vars_box.setVal(1.0, c_h22);
        grchombo_vars_box.setVal(1.0, c_h33);
        grchombo_vars_box.setVal(1.0, c_lapse);

        // now non constant terms by location
        Box this_box = grchombo_vars_box.box();
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            set_non_const_output_cell(multigrid_vars_box, grchombo_vars_box,
                                      bit(), a_dx, a_params);
        }

        // finally non-constant boundary ghosts
        IntVect offset_lo, offset_hi;
        a_grchombo_boundaries.get_box_offsets(offset_lo, offset_hi, this_box);

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        a_grchombo_boundaries.remove_outer_ghost_cells(this_box);

        // get the boundary box (may be Empty)
        for (int idir = 0; idir < SpaceDim; ++idir)
        {
            if (!a_params.periodic[idir])
            {
                for (SideIterator sit; sit.ok(); ++sit)
                {
                    Box boundary_box = a_grchombo_boundaries.get_boundary_box(
                        sit(), idir, offset_lo, offset_hi, this_box);

                    // now we have the appropriate box, fill it!
                    BoxIterator bbit(boundary_box);
                    for (bbit.begin(); bbit.ok(); ++bbit)
                    {
                        set_non_const_output_cell(multigrid_vars_box,
                                                  grchombo_vars_box, bbit(),
                                                  a_dx, a_params);
                    } // end loop through boundary box
                }     // end loop over sides
            }         // end if (periodic[idir])
        }             // end loop over directions
    }                 // end loop over boxes
} // end set_output_data

void set_non_const_output_cell(const FArrayBox &a_multigrid_vars_box,
                               FArrayBox &a_grchombo_vars_box,
                               const IntVect &a_iv, const RealVect &a_dx,
                               const PoissonParameters &a_params)
{
    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    // GRChombo conformal factor chi = psi^-4
    Real psi_bh = get_psi_brill_lindquist(loc, a_params);
    Real chi = pow(a_multigrid_vars_box(a_iv, c_psi_reg) + psi_bh, -4.0);
    a_grchombo_vars_box(a_iv, c_chi) = chi;

    // use pre-collapsed lapse
    a_grchombo_vars_box(a_iv, c_lapse) = sqrt(chi);

    // Copy Aij across - note this is now \tilde Aij not \bar
    // Aij
    Real factor = pow(chi, 1.5);
    a_grchombo_vars_box(a_iv, c_A11) =
        a_multigrid_vars_box(a_iv, c_A11_0) * factor;
    a_grchombo_vars_box(a_iv, c_A12) =
        a_multigrid_vars_box(a_iv, c_A12_0) * factor;
    a_grchombo_vars_box(a_iv, c_A13) =
        a_multigrid_vars_box(a_iv, c_A13_0) * factor;
    a_grchombo_vars_box(a_iv, c_A22) =
        a_multigrid_vars_box(a_iv, c_A22_0) * factor;
    a_grchombo_vars_box(a_iv, c_A23) =
        a_multigrid_vars_box(a_iv, c_A23_0) * factor;
    a_grchombo_vars_box(a_iv, c_A33) =
        a_multigrid_vars_box(a_iv, c_A33_0) * factor;
}
