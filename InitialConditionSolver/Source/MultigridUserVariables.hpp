/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MULTIGRIDUSERVARIABLES_HPP
#define MULTIGRIDUSERVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_psi_reg, // the regular part of psi_0 not including the singular
               // Brill-Lindquist part which is added separately
               // See Alcubierre p112 eq. (3.4.29)

    c_A11_0,
    c_A12_0,
    c_A13_0,
    c_A22_0,
    c_A23_0,
    c_A33_0,

    NUM_MULTIGRID_VARS
};

namespace MultigridUserVariables
{
static constexpr char const *variable_names[NUM_MULTIGRID_VARS] = {
    "psi",

    "A11_0", "A12_0", "A13_0", "A22_0", "A23_0", "A33_0"};
}

#endif /* MULTIGRIDUSERVARIABLES_HPP */
