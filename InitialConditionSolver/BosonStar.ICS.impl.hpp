/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTAR_ICS_HPP_)
#error "This file should only be included through BosonStar.hpp"
#endif

#ifndef BOSONSTAR_ICS_IMPL_HPP_
#define BOSONSTAR_ICS_IMPL_HPP_

#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "BosonStarBinarySearch.hpp" //for BosonStarBinarySearch class

inline BosonStar::BosonStar(BosonStar_params_t a_params_BosonStar,
                    Potential::params_t a_params_potential, double a_G_Newton,
                    int a_verbosity)
    : m_1d_sol(a_params_BosonStar, a_params_potential, a_G_Newton, a_verbosity),
    m_G_Newton(a_G_Newton), m_params_BosonStar(a_params_BosonStar),
    m_params_potential(a_params_potential), m_verbosity(a_verbosity)

{
}

void BosonStar::compute_1d_solution(const double a_max_radius)
{
    try
    {
        BosonStarBinarySearch<initial_data_t, initial_state_t> binary_search(
        m_params_BosonStar, m_params_potential, m_G_Newton, m_verbosity);

        binary_search.shoot();
        auto sol = binary_search.getShootedSolution();

        m_1d_sol.makeFromPolarArealSolution(sol, a_max_radius);
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

#endif /* BOSONSTAR_ICS_IMPL_HPP_ */
