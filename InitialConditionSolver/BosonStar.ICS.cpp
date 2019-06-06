/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "REAL.H"
#include "BosonStar.ICS.hpp"
#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "BosonStarBinarySearch.hpp" //for BosonStarBinarySearch class


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
