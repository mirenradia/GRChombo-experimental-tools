#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SETBCS_H_
#define _SETBCS_H_

#include <iostream>
using std::cerr;

#include "AMRMultiGrid.H"
#include "BCFunc.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "VariableCoeffPoissonOperator.H"



// This class sets the value for BCs to its member variable m_value
class BCValue : public BCValueFunction
{
  private:
    double m_value;

  public:
    // Constructor
    BCValue(const double a_value) : m_value(a_value) {}

    // The operator actually called for the BC value
    virtual void operator()(Real*           a_pos,
                            int*            a_dir,
                            Side::LoHiSide* a_side,
                            Real*           a_value) override;
};


// BC Class
class SetBCs : public BCFunction
{
  public:
    struct params_t
    {
        std::vector<int> bc_lo;
        std::vector<int> bc_hi;
        double value;
    };

  private:
    int m_NL_iter;
    params_t m_params;
    RefCountedPtr<BCValueFunction> m_bc_value_ptr;
    RefCountedPtr<BCValueFunction> m_bc_zero_ptr;

    inline void printSideBCs(Side::LoHiSide a_side);

    inline void SetSideBCs(FArrayBox &a_state, const Box &a_valid,
                           const ProblemDomain &a_domain, Real a_dx,
                           bool a_homogeneous, Side::LoHiSide a_side);

  public:
    // constructor
    SetBCs(const int a_NL_iter, const params_t &a_params);

    void printBCs();

    virtual void operator() (FArrayBox &a_state, const Box &a_valid,
                             const ProblemDomain &a_domain, Real a_dx,
                             bool a_homogeneous) override;
};


#endif
