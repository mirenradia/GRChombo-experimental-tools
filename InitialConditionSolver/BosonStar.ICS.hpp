/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTAR_ICS_HPP_
#define BOSONSTAR_ICS_HPP_

#include "ComplexPotential.hpp"
#include "BosonStarParams.hpp"
#include "BosonStarIsotropicSolution.hpp"
#include "PoissonParameters.H"
#include <vector>
#include "parstream.H" //gives pout

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling. This is copied from GRChombo with some of the
//! GRChombo integration (i.e. compute member function) removed
class BosonStar
{
    typedef std::vector<double> initial_state_t;
    template<class T>
    using initial_data_t = std::vector<T>;
    friend class BinaryBS;

public:
    //! The constructor
    BosonStar(BosonStar_params_t a_params_BosonStar,
        Potential::params_t a_params_potential, double a_G_Newton,
        int a_verbosity)
    : m_1d_sol(a_params_BosonStar, a_params_potential, a_G_Newton, a_verbosity),
    m_G_Newton(a_G_Newton), m_params_BosonStar(a_params_BosonStar),
    m_params_potential(a_params_potential), m_verbosity(a_verbosity)
    {
    }

    //! Computes the 1d solution and stores in m_1d_sol
    void compute_1d_solution(const double a_max_radius);

    BosonStarIsotropicSolution<initial_data_t, initial_state_t> m_1d_sol; /*<
    The object that stores the solution found by the 1d ODE integrator */

    friend void compute_boson_star_profiles(BosonStar &a_boson_star1,
                                            BosonStar &a_boson_star2,
                                            PoissonParameters &a_params);

protected:
    double m_G_Newton;
    BosonStar_params_t m_params_BosonStar; //!< The complex scalar field params
    Potential::params_t m_params_potential; //!< The potential params
    int m_verbosity;
};

#endif /* BOSONSTAR_ICS_HPP_ */
