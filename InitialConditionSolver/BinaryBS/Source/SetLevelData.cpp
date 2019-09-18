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
#include "LoadBalance.H"
#include "LoHiSide.H"
#include "PoissonParameters.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

// This takes an IntVect and writes the physical coordinates to a RealVect
void get_loc(RealVect &a_out_loc, const IntVect &a_iv,
             const RealVect &a_dx, const PoissonParameters &a_params)
{
    a_out_loc = a_iv + 0.5 * RealVect::Unit;
    a_out_loc *= a_dx;
    a_out_loc -= a_params.domainLength / 2.0;
}

// set initial guess value for the conformal factor psi
// defined by \gamma_ij = \psi^4 \tilde \gamma_ij and complex scalar field phi
// This has been/will be modified to take the superposition of two boson stars
// as the initial guess where
// gamma_ij = gamma^(1)_ij + gamma^(2)_ij - delta_ij
void set_initial_conditions(LevelData<FArrayBox> &a_multigrid_vars,
                            LevelData<FArrayBox> &a_dpsi,
                            GRChomboBCs &a_grchombo_boundaries,
                            const RealVect &a_dx,
                            const PoissonParameters &a_params,
                            BosonStar &a_boson_star1,
                            BosonStar &a_boson_star2)
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
            set_initial_multigrid_cell(multigrid_vars_box, dpsi_box,
                bit(), a_dx, a_params, a_boson_star1, a_boson_star2);
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
                            bbit(), a_dx, a_params, a_boson_star1,
                            a_boson_star2);
                    } // end loop through boundary box
                } // end loop over sides
            } // end if (periodic[idir])
        } // end loop over directions
    }
} // end set_initial_conditions

void set_initial_multigrid_cell(FArrayBox &a_multigrid_vars_box,
                                FArrayBox &a_dpsi_box,
                                const IntVect &a_iv,
                                const RealVect &a_dx,
                                const PoissonParameters &a_params,
                                BosonStar &a_boson_star1,
                                BosonStar &a_boson_star2)
{
    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    // work out the displacement from the centre of each star
    RealVect loc1, loc2;
    for(int idir = 0; idir < SpaceDim; ++idir)
    {
        loc1[idir] = loc[idir]
            - a_params.boson_star1_params.star_centre[idir];
        loc2[idir] = loc[idir]
            - a_params.boson_star2_params.star_centre[idir];
    }
    Real r1, r2;
    r1 = loc1.vectorLength();
    r2 = loc2.vectorLength();

    // CCZ4 variables
    Real chi1 = a_boson_star1.m_1d_sol.m_chi(r1);
    Real chi2 = a_boson_star2.m_1d_sol.m_chi(r2);
    Real lapse1 = a_boson_star1.m_1d_sol.m_lapse(r1);
    Real lapse2 = a_boson_star2.m_1d_sol.m_lapse(r2);

    // superposed conformal factor
    // 1/chi = 1/chi1 + 1/chi2 - chi_subtraction_constant

    Real chi = (chi1 * chi2) / (chi1 + chi2
            - a_params.chi_subtraction_constant * chi1 * chi2);
    Real lapse = sqrt(lapse1 * lapse1 + lapse2 * lapse2 - 1.0);
    a_multigrid_vars_box(a_iv, c_psi_0) = pow(chi, -0.25);
    a_multigrid_vars_box(a_iv, c_lapse_0) = lapse;

    // Matter superposition
    Real phase1 = a_params.boson_star1_params.phase;
    Real phase2 = a_params.boson_star2_params.phase;
    Real frequency1 = a_boson_star1.m_1d_sol.m_frequency_over_mass
                        * a_params.potential_params.scalar_mass;
    Real frequency2 = a_boson_star2.m_1d_sol.m_frequency_over_mass
                        * a_params.potential_params.scalar_mass;
    Real mod_phi1 = a_boson_star1.m_1d_sol.m_phi(r1);
    Real mod_phi2 = a_boson_star2.m_1d_sol.m_phi(r2);
    a_multigrid_vars_box(a_iv, c_phi_Re_0) = mod_phi1 * cos(phase1)
                                         + mod_phi2 * cos(phase2);
    a_multigrid_vars_box(a_iv, c_phi_Im_0) = mod_phi1 * sin(phase1)
                                         + mod_phi2 * sin(phase2);
    a_multigrid_vars_box(a_iv, c_Pi_Re_0)
        = frequency1 * mod_phi1 * sin(phase1) / lapse1
        + frequency2 * mod_phi2 * sin(phase2) / lapse2;
    a_multigrid_vars_box(a_iv, c_Pi_Im_0)
        = -frequency1 * mod_phi1 * cos(phase1) / lapse1
         - frequency2 * mod_phi2 * cos(phase2) / lapse2;

    // dpsi is initially zero
    a_dpsi_box(a_iv, 0) = 0.0;
} //end set set_initial_multigrid_cell

// set the rhs source for the poisson eqn
void set_rhs(LevelData<FArrayBox> &a_rhs,
             LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx,
             const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    Potential potential(a_params.potential_params);
    DataIterator dit = a_rhs.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &rhs_box = a_rhs[dit()];
        rhs_box.setVal(0.0, 0);
        Box this_box = rhs_box.box(); // no ghost cells
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();

            // rhs = -laplacian of psi_0 - pi G psi_0^5 (|Pi|^2
            //       + psi_0^{-4} grad_phi_sq + V(|phi|^2))
            Real laplacian_of_psi_0, grad_phi_sq, V_of_mod_phi_sq;
            set_laplacian_psi(laplacian_of_psi_0, iv, multigrid_vars_box, a_dx);
            set_grad_phi_sq(grad_phi_sq, iv, multigrid_vars_box, a_dx);
            potential.compute_potential(V_of_mod_phi_sq,
                                        multigrid_vars_box(iv, c_phi_Re_0),
                                        multigrid_vars_box(iv, c_phi_Im_0));

            Real Pi_Re = multigrid_vars_box(iv, c_Pi_Re_0);
            Real Pi_Im = multigrid_vars_box(iv, c_Pi_Im_0);
            Real mod_Pi_sq = Pi_Re * Pi_Re + Pi_Im * Pi_Im;
            Real psi_0 = multigrid_vars_box(iv, c_psi_0);

            rhs_box(iv, 0) = -laplacian_of_psi_0 - M_PI * a_params.G_Newton *
                             (pow(psi_0, 5.0) * (mod_Pi_sq + V_of_mod_phi_sq) +
                              psi_0 * grad_phi_sq);
        }
    }
} // end set_rhs

// computes the Laplacian of psi at a point in a box
inline void set_laplacian_psi(Real &laplacian_of_psi_0,
                              const IntVect &a_iv,
                              const FArrayBox &a_multigrid_vars_box,
                              const RealVect &a_dx)
{
    laplacian_of_psi_0 = 0.0;
    for(int idir = 0; idir < SpaceDim; ++idir)
    {
        IntVect iv_offset1 = a_iv;
        IntVect iv_offset2 = a_iv;
        iv_offset1[idir] -= 1;
        iv_offset2[idir] += 1;

        // 2nd order stencil for now
        Real d2psi_dxdx = 1.0 / (a_dx[idir] * a_dx[idir]) *
                          (+1.0 * a_multigrid_vars_box(iv_offset2, c_psi_0)
                           -2.0 * a_multigrid_vars_box(a_iv, c_psi_0)
                           +1.0 * a_multigrid_vars_box(iv_offset1, c_psi_0));
        laplacian_of_psi_0 += d2psi_dxdx;
    }
} // end set_laplacian_psi

// computes the gradient of the complex scalar field squared at a point in a box
// i.e. delta^{ij} d_i phi^* d_j phi
// = delta^{ij}(d_i phi_Re d_j phi_Re + d_i phi_Im d_j phi_Im)
inline void set_grad_phi_sq(Real &grad_phi_sq,
                            const IntVect &a_iv,
                            const FArrayBox &a_multigrid_vars_box,
                            const RealVect &a_dx)
{
    grad_phi_sq = 0.0;
    for(int idir = 0; idir < SpaceDim; ++idir)
    {
        IntVect iv_offset1 = a_iv;
        IntVect iv_offset2 = a_iv;
        iv_offset1[idir] -= 1;
        iv_offset2[idir] += 1;
        Real dx_inv = 1 / a_dx[idir];

        // 2nd order stencils for now
        Real dphi_Re_dx = 0.5 * dx_inv *
                          (a_multigrid_vars_box(iv_offset2, c_phi_Re_0) -
                           a_multigrid_vars_box(iv_offset1, c_phi_Re_0));
        Real dphi_Im_dx = 0.5 * dx_inv *
                          (a_multigrid_vars_box(iv_offset2, c_phi_Im_0) -
                           a_multigrid_vars_box(iv_offset1, c_phi_Im_0));
        grad_phi_sq += dphi_Re_dx * dphi_Re_dx + dphi_Im_dx * dphi_Im_dx;
    }
} //end set_grad_phi_sq

// set the regrid condition - abs value of this drives AMR regrid
void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                          LevelData<FArrayBox> &a_multigrid_vars,
                          const RealVect &a_dx,
                          const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    Potential potential(a_params.potential_params);
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

            Real grad_phi_sq, V_of_mod_phi_sq;
            set_grad_phi_sq(grad_phi_sq, iv, multigrid_vars_box, a_dx);
            potential.compute_potential(V_of_mod_phi_sq,
                                        multigrid_vars_box(iv, c_phi_Re_0),
                                        multigrid_vars_box(iv, c_phi_Im_0));

            Real Pi_Re = multigrid_vars_box(iv, c_Pi_Re_0);
            Real Pi_Im = multigrid_vars_box(iv, c_Pi_Im_0);
            Real mod_Pi_sq = Pi_Re * Pi_Re + Pi_Im * Pi_Im;
            Real psi_0 = multigrid_vars_box(iv, c_psi_0);

            // the condition is similar to the 12*rhs but without the laplacian
            condition_box(iv, 0) = 12.0 * M_PI * a_params.G_Newton *
                             (pow(psi_0, 5.0) * (mod_Pi_sq + V_of_mod_phi_sq) +
                              psi_0 * grad_phi_sq);

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
            multigrid_vars_box(iv, c_psi_0) += dpsi_box(iv, 0);
        }
    }
}

// The coefficient of the I operator on dpsi
void set_a_coef(LevelData<FArrayBox> &a_aCoef,
                LevelData<FArrayBox> &a_multigrid_vars,
                const PoissonParameters &a_params, const RealVect &a_dx)
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    Potential potential(a_params.potential_params);
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

            Real grad_phi_sq, V_of_mod_phi_sq;
            set_grad_phi_sq(grad_phi_sq, iv, multigrid_vars_box, a_dx);
            potential.compute_potential(V_of_mod_phi_sq,
                                        multigrid_vars_box(iv, c_phi_Re_0),
                                        multigrid_vars_box(iv, c_phi_Im_0));


            Real Pi_Re = multigrid_vars_box(iv, c_Pi_Re_0);
            Real Pi_Im = multigrid_vars_box(iv, c_Pi_Im_0);
            Real mod_Pi_sq = Pi_Re * Pi_Re + Pi_Im * Pi_Im;
            Real psi_0 = multigrid_vars_box(iv, c_psi_0);

            aCoef_box(iv, 0) = M_PI * a_params.G_Newton *
                               (5 * pow(psi_0, 4.0) * (mod_Pi_sq
                                + V_of_mod_phi_sq) + grad_phi_sq);
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
        // Conformally flat
        grchombo_vars_box.setVal(1.0, c_h11);
        grchombo_vars_box.setVal(1.0, c_h22);
        grchombo_vars_box.setVal(1.0, c_h33);
        //grchombo_vars_box.setVal(1.0, c_lapse);

        // now non constant terms by location
        Box this_box = grchombo_vars_box.box();
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            set_non_const_output_cell(multigrid_vars_box,
                grchombo_vars_box, bit());
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
                            grchombo_vars_box, bbit());
                    } // end loop through boundary box
                } // end loop over sides
            } // end if (periodic[idir])
        } // end loop over directions
    } // end loop over boxes
} // end set_output_data

void set_non_const_output_cell(const FArrayBox &a_multigrid_vars_box,
                               FArrayBox &a_grchombo_vars_box,
                               const IntVect &a_iv)
{
    // GRChombo conformal factor chi = psi^-4
    Real chi = pow(a_multigrid_vars_box(a_iv, c_psi_0), -4.0);
    a_grchombo_vars_box(a_iv, c_chi) = chi;

    // Copy lapse
    a_grchombo_vars_box(a_iv, c_lapse) = a_multigrid_vars_box(a_iv, c_lapse_0);

    // Copy phi and Pi
    a_grchombo_vars_box(a_iv, c_phi_Re)
        = a_multigrid_vars_box(a_iv, c_phi_Re_0);
    a_grchombo_vars_box(a_iv, c_phi_Im)
        = a_multigrid_vars_box(a_iv, c_phi_Im_0);
    a_grchombo_vars_box(a_iv, c_Pi_Re)
        = a_multigrid_vars_box(a_iv, c_Pi_Re_0);
    a_grchombo_vars_box(a_iv, c_Pi_Im)
        = a_multigrid_vars_box(a_iv, c_Pi_Im_0);
}
