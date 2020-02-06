/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETBLACKHOLE_H_
#define SETBLACKHOLE_H_

#include "PoissonParameters.H"

inline RealVect get_bh_displacement(const RealVect &loc,
                                    const PoissonParameters &a_params)
{
    return (loc - a_params.bh_offset);
}

// Set Aij Bowen York data
// see Alcubierre pg 110 eqn (3.4.22)
Real get_Aij(const int i, const int j, const RealVect &a_bh_displacement,
             const RealVect &a_bh_spin, const RealVect &a_bh_momentum)
{
    std::array<std::array<std::array<double, 3>, 3>, 3> epsilon = {0.};
    epsilon[0][1][2] = 1.0;
    epsilon[1][2][0] = 1.0;
    epsilon[2][0][1] = 1.0;
    epsilon[0][2][1] = -1.0;
    epsilon[2][1][0] = -1.0;
    epsilon[1][0][2] = -1.0;

    // relabel n^i (in Alcubierre) to m^i to avoid confusion with normal to
    // hypersurface
    Real r = a_bh_displacement.vectorLength();
    Real rr = a_bh_displacement.radSquared();
    RealVect m = a_bh_displacement / r;
    const RealVect &P = a_bh_momentum;
    const RealVect &S = a_bh_spin;

    Real Aij = 1.5 / rr * (m[i] * P[j] + m[j] * P[i]);
    for (int k = 0; k < SpaceDim; k++)
    {
        Aij += 1.5 / rr * (m[i] * m[j] - Real(i == j)) * P[k] * m[k];

        for (int l = 0; l < SpaceDim; l++)
        {
            Aij += -3.0 / (rr * r) *
                       (epsilon[i][l][k] * m[j] + epsilon[j][l][k] * m[i]) *
                       m[l] * S[k];
        }
    }

    return Aij;
}

void set_black_hole_Aij(FArrayBox &multigrid_vars_box, const IntVect &iv,
                        const RealVect &loc, const PoissonParameters &a_params)
{
    // set the BH values - location
    RealVect bh_displacement = get_bh_displacement(loc, a_params);

    // set the Aij vars
    multigrid_vars_box(iv, c_A11_0) =
        get_Aij(0, 0, bh_displacement, a_params.bh_spin, a_params.bh_momentum);
    multigrid_vars_box(iv, c_A22_0) =
        get_Aij(1, 1, bh_displacement, a_params.bh_spin, a_params.bh_momentum);
    multigrid_vars_box(iv, c_A33_0) =
        get_Aij(2, 2, bh_displacement, a_params.bh_spin, a_params.bh_momentum);
    multigrid_vars_box(iv, c_A12_0) =
        get_Aij(0, 1, bh_displacement, a_params.bh_spin, a_params.bh_momentum);
    multigrid_vars_box(iv, c_A13_0) =
        get_Aij(0, 2, bh_displacement, a_params.bh_spin, a_params.bh_momentum);
    multigrid_vars_box(iv, c_A23_0) =
        get_Aij(1, 2, bh_displacement, a_params.bh_spin, a_params.bh_momentum);
}

// get the singular Brill-lindquist part of psi
Real get_psi_bl(const RealVect &loc, const PoissonParameters &a_params)
{
    RealVect bh_displacement = get_bh_displacement(loc, a_params);
    Real r = bh_displacement.vectorLength();

    return 0.5 * a_params.bh_bare_mass / r;
}

#endif /* SETBLACKHOLE_H_ */
