/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MULTIGRIDUSERVARIABLES_HPP
#define MULTIGRIDUSERVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_psi_0,

    c_lapse_0,

    c_phi_Re_0, // real part of scalar field
    c_phi_Im_0, // imaginary part of scalar field
    c_Pi_Re_0, // real part of auxiliary variable Pi = -L_n phi
    c_Pi_Im_0, // imaginary part of auxiliary variable Pi = -L_n phi

    NUM_MULTIGRID_VARS
};

namespace MultigridUserVariables
{
static constexpr char const *variable_names[NUM_MULTIGRID_VARS] = {
    "psi_0",

    "lapse_0",

    "phi_Re_0", "phi_Im_0", "Pi_Re_0", "Pi_Im_0"};
}

#endif /* MULTIGRIDUSERVARIABLES_HPP */
