/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXPOTENTIAL_HPP_
#define COMPLEXPOTENTIAL_HPP_

class Potential
{
  public:
    struct params_t
    {
        double scalar_mass;
        double phi4_coeff;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    void compute_potential(Real &V_of_modulus_phi_squared,
                           const Real &phi_Re_here,
                           const Real &phi_Im_here) const
    {
        // First calculate |phi|^2
        Real modulus_phi_squared =
            phi_Re_here * phi_Re_here + phi_Im_here * phi_Im_here;

        // The potential value at phi (note the convention with factors of 1/2)
        // m^2 |phi|^2 + lambda/2 |phi|^4
        V_of_modulus_phi_squared =
            m_params.scalar_mass * m_params.scalar_mass * modulus_phi_squared +
            0.5 * m_params.phi4_coeff * modulus_phi_squared
            * modulus_phi_squared;
    }
};

#endif /* COMPLEXPOTENTIAL_HPP_ */
