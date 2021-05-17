//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj_H
#define Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj_H

#include "tthAnalysis/tthMEM/interface/mg5/read_slha.h"

#include <complex>

#ifndef MDL_MT_masses
#define MDL_MT_masses
extern "C"
{
  extern struct
  {
    double MDL_MT;
  } masses_;
}
#endif  // MDL_MT_masses

class Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj
{
  public:

    static Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj * getInstance();

    // Define "zero"
    double zero, ZERO;
    // Model parameters independent of aS
    double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_ymb, aS,
        mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB, mdl_cy,
        mdl_ctr, mdl_a2, mdl_a1, mdl_c2, mdl_MZ__exp__2, mdl_MZ__exp__4,
        mdl_sqrt__2, mdl_MH__exp__2, mdl_cy__exp__2, mdl_aEW, mdl_MW,
        mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2,
        mdl_sw, mdl_g1, mdl_gw, mdl_vev, mdl_vev__exp__2, mdl_lam, mdl_yb,
        mdl_yt, mdl_ytau, mdl_muH, mdl_ee__exp__2, mdl_sw__exp__2,
        mdl_cw__exp__2;
    std::complex<double> mdl_complexi;
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G, mdl_G__exp__2;
    // Model couplings independent of aS
    std::complex<double> GC_13, GC_28, GC_29, GC_31;
    // Model couplings dependent on aS
    std::complex<double> GC_38, GC_22, GC_25, GC_34, GC_37, GC_36;

    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha);
    // Set couplings that are unchanged during the run
    void setIndependentCouplings();
    // Set parameters that are changed event by event
    void setDependentParameters();
    // Set couplings that are changed event by event
    void setDependentCouplings(double & kl, double & kt, double & c2, double & cg, double & c2g);

    // Print parameters that are unchanged during the run
    void printIndependentParameters();
    // Print couplings that are unchanged during the run
    void printIndependentCouplings();
    // Print parameters that are changed event by event
    void printDependentParameters();
    // Print couplings that are changed event by event
    void printDependentCouplings();


  private:
    static Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj * instance;
};

#endif  // Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj_H
