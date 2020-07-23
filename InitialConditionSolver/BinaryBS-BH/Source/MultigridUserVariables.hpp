/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MULTIGRIDUSERVARIABLES_HPP
#define MULTIGRIDUSERVARIABLES_HPP

#include <array>
#include <string>

// assign an enum to each variable
enum
{
    c_psi_reg, // the _reg denotes that this does not contain the singular
               // Brill-Lindquist part of psi which is added on when needed

    c_A11_0,
    c_A12_0,
    c_A13_0,
    c_A22_0,
    c_A23_0,
    c_A33_0,

    // c_lapse_0,

    c_phi_Re_0, // real part of scalar field
    c_phi_Im_0, // imaginary part of scalar field
    c_Pi_Re_0,  // real part of auxiliary variable Pi = -L_n phi
    c_Pi_Im_0,  // imaginary part of auxiliary variable Pi = -L_n phi

    NUM_MULTIGRID_VARS
};

namespace MultigridUserVariables
{
static const std::array<std::string, NUM_MULTIGRID_VARS> variable_names = {
    "psi_reg",

    "A11_0", "A12_0", "A13_0", "A22_0", "A23_0", "A33_0",

    //"lapse_0",

    "phi_Re_0", "phi_Im_0", "Pi_Re_0", "Pi_Im_0"};
}

#endif /* MULTIGRIDUSERVARIABLES_HPP */
