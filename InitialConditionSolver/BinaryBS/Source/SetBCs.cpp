#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SetBCs.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

void BCValue::operator()(Real*           a_pos,
                         int*            a_dir,
                         Side::LoHiSide* a_side,
                         Real*           a_value)
{
    *(a_value) = m_value;
}

SetBCs::SetBCs(const int a_NL_iter, const params_t &a_params)
    : m_NL_iter(a_NL_iter), m_params(a_params),
      m_bc_value_ptr(new BCValue(m_params.value)),
      m_bc_zero_ptr(new BCValue(0.0)) {}

void SetBCs::printSideBCs(Side::LoHiSide a_side)
{
    std::string side_str = (a_side == Side::Lo) ? "low" : "high";
    std::vector<int>& bc = (a_side == Side::Lo) ?
        m_params.bc_lo : m_params.bc_hi;
    for(int idir = 0; idir < CH_SPACEDIM; ++idir)
    {
        std::string bc_str;
        switch(bc[idir])
        {
        case 0:
            bc_str = "Constant Dirichlet";
            break;
        case 1:
            bc_str = "Constant Neumann";
            break;
        case 2:
            bc_str = "Periodic";
            break;
        case 3:
            bc_str = "NL Dirichlet";
            break;
        default:
            MayDay::Error("bogus bc flag");
        }
        pout() << bc_str << " bcs imposed for " << side_str
               << " side in direction " << idir << endl;
    }
}

void SetBCs::printBCs()
{
    if (m_NL_iter == 0)
    {
        printSideBCs(Side::Lo);
        printSideBCs(Side::Hi);
    }
}

void SetBCs::SetSideBCs(FArrayBox &a_state, const Box &a_valid,
                        const ProblemDomain &a_domain, Real a_dx,
                        bool a_homogeneous, Side::LoHiSide a_side)
{
    if (!a_domain.domainBox().contains(a_state.box()))
    {
        Box valid = a_valid;
        std::vector<int>& bc = (a_side == Side::Lo) ?
            m_params.bc_lo : m_params.bc_hi;
        BCValueHolder bc_value_holder(m_bc_value_ptr);
        BCValueHolder bc_value_zero_holder(m_bc_zero_ptr);
        for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
            // periodic? If not, check if Dirichlet or Neumann
            if (!a_domain.isPeriodic(idir))
            {
                Box ghostBox = adjCellBox(valid, idir, a_side, 1);
                if (!a_domain.domainBox().contains(ghostBox))
                {
                    switch(bc[idir])
                    {
                    case 0:
                        DiriBC(a_state, valid, a_dx, a_homogeneous,
                               bc_value_holder, idir, a_side);
                        break;
                    case 1:
                        NeumBC(a_state, valid, a_dx, a_homogeneous,
                               bc_value_holder, idir, a_side);
                        break;
                    case 2:
                        break;// periodic case so no need to do anything
                    case 3:
                        if (m_NL_iter == 0)
                        {
                            DiriBC(a_state, valid, a_dx, a_homogeneous,
                                   bc_value_holder, idir, a_side);
                        }
                        else
                        {
                            DiriBC(a_state, valid, a_dx, a_homogeneous,
                                   bc_value_zero_holder, idir, a_side);
                        }
                        break;
                    default:
                        MayDay::Error("bogus bc flag side");
                    }
                }
            } // else - is periodic
        } // close for idir
    }
}

void SetBCs::operator() (FArrayBox &a_state, const Box &a_valid,
                         const ProblemDomain &a_domain, Real a_dx,
                         bool a_homogeneous)
{
    SetSideBCs(a_state, a_valid, a_domain, a_dx, a_homogeneous, Side::Lo);
    SetSideBCs(a_state, a_valid, a_domain, a_dx, a_homogeneous, Side::Hi);
}
