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
#include "MyPhiFunction.H"
#include "PoissonParameters.H"
#include "SetLevelDataF_F.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Boson Star includes
#include "BosonStar.ICS.hpp"
#include "ComplexPotential.hpp"

// Set various LevelData functions across the grid

// set initial guess value for the conformal factor psi
// defined by \gamma_ij = \psi^4 \tilde \gamma_ij and complex scalar field phi
// This has been/will be modified to take the superposition of two boson stars
// as the initial guess where
// gamma_ij = gamma^(1)_ij + gamma^(2)_ij - delta_ij
void set_initial_conditions(LevelData<FArrayBox> &a_multigrid_vars,
                            LevelData<FArrayBox> &a_dpsi, const RealVect &a_dx,
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
        Box b = multigrid_vars_box.box();
        BoxIterator bit(b);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();
            RealVect loc(iv + 0.5 * RealVect::Unit);
            loc *= a_dx;
            loc -= a_params.domainLength / 2.0;

            // work out the displacement from the centre of each star
            RealVect loc1, loc2;
            for(int idir = 0; idir < SpaceDim; ++idir)
            {
                loc1[idir] = loc[idir]
                    - a_boson_star1.params_BosonStar.star_centre[idir];
                loc2[idir] = loc[idir]
                    - a_boson_star2.params_BosonStar.star_centre[idir];
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
            // 1/chi = 1/chi1 + 1/chi2 - 1
            Real chi = (chi1 * chi2) / (chi1 + chi2 - chi1 * chi2);
            Real lapse = sqrt(lapse1 * lapse1 + lapse2 * lapse2 - 1.0);
            multigrid_vars_box(iv, c_psi_0) = pow(chi, -0.25);
            multigrid_vars_box(iv, c_lapse_0) = lapse;

            // Matter superposition
            Real phase1 = a_boson_star1.m_params_BosonStar.phase;
            Real phase2 = a_boson_star2.m_params_BosonStar.phase;
            Real frequency1 = a_boson_star1.m_1d_sol.m_frequency_over_mass
                                * a_boson_star1.m_params_potential.scalar_mass;
            Real frequency2 = a_boson_star2.m_1d_sol.m_frequency_over_mass
                                * a_boson_star2.m_params_potential.scalar_mass;
            Real mod_phi1 = a_boson_star1.m_1d_sol.m_phi(r1);
            Real mod_phi2 = a_boson_star2.m_1d_sol.m_phi(r2);
            multigrid_vars_box(iv, c_phi_Re_0) = mod_phi1 * cos(phase1)
                                                 + mod_phi2 * cos(phase2);
            multigrid_vars_box(iv, c_phi_Im_0) = mod_phi1 * sin(phase1)
                                                 + mod_phi2 * sin(phase2);
            multigrid_vars_box(iv, c_Pi_Re_0)
                = frequency1 * mod_phi1 * sin(phase1) / lapse1
                + frequency2 * mod_phi2 * sin(phase2) / lapse2;
            multigrid_vars_box(iv, c_Pi_Im_0)
                = -frequency1 * mod_phi1 * cos(phase1) / lapse1
                 - frequency2 * mod_phi2 * cos(phase2) / lapse2;

            // dpsi is initially zero
            dpsi_box(iv, 0) = 0.0;
        }
    }
} // end set_initial_conditions

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
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();

            // rhs = -laplacian of psi_0 - pi G psi_0^5 (|Pi|^2
            //       + psi_0^{-4} grad_phi_sq + V(|phi|^2))
            Real laplacian_of_psi_0, grad_phi_sq, V_of_mod_phi_sq;
            set_laplacian_psi(laplacian_of_psi_0, iv, multigrid_vars_box, a_dx);
            set_grad_phi_sq(grad_phi_sq, iv, multigrid_vars_box, a_dx);
            compute_potential(V_of_mod_phi_sq,
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
    laplacian_of_psi = 0.0;
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
        laplacian_of_psi += d2psi_dxdx;
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
                           a_multigrid_vars_box(iv_offset1, c_phi_Re_0))
        Real dphi_Im_dx = 0.5 * dx_inv *
                          (a_multigrid_vars_box(iv_offset2, c_phi_Im_0) -
                           a_multigrid_vars_box(iv_offset1, c_phi_Im_0))
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

            Real grad_phi_sq, mod_Pi_sq, V_of_mod_phi_sq;
            set_grad_phi_sq(grad_phi_sq, iv, multigrid_vars_box, a_dx);
            compute_potential(V_of_mod_phi_sq,
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

            Real grad_phi_sq, mod_Pi_sq, V_of_mod_phi_sq;
            set_grad_phi_sq(grad_phi_sq, iv, multigrid_vars_box, a_dx);
            compute_potential(V_of_mod_phi_sq,
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
            IntVect iv = bit();

            // GRChombo conformal factor chi = psi^-4
            Real chi = pow(multigrid_vars_box(iv, c_psi_0), -4.0);
            grchombo_vars_box(iv, c_chi) = chi;

            // Copy lapse
            grchombo_vars_box(iv, c_lapse) = multigrid_vars_box(iv, c_lapse_0);

            // Copy phi and Pi
            grchombo_vars_box(iv, c_phi_Re)
                = multigrid_vars_box(iv, c_phi_Re_0);
            grchombo_vars_box(iv, c_phi_Im)
                = multigrid_vars_box(iv, c_phi_Im_0);
            grchombo_vars_box(iv, c_Pi_Re)
                = multigrid_vars_box(iv, c_Pi_Re_0);
            grchombo_vars_box(iv, c_Pi_Im)
                = multigrid_vars_box(iv, c_Pi_Im_0);
        }
    }
}
