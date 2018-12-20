//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "hhAnalysis/bbwwMEM/interface/mg5/parameters/Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"

#include <iomanip>

// Initialize static instance
Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj * Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::instance = nullptr;

// Function to get static instance - only one instance per program
Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj * Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::getInstance()
{
  if (! instance)
    instance = new Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj();

  return instance; 
}

void Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::setIndependentParameters(SLHAReader& slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  mdl_WH = slha.get_block_entry("decay", 25, 6.382339e-03); 
  mdl_WW = slha.get_block_entry("decay", 24, 2.085000e+00); 
  mdl_WZ = slha.get_block_entry("decay", 23, 2.495200e+00); 
  mdl_WT = slha.get_block_entry("decay", 6, 1.508336e+00); 
  mdl_ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  mdl_ymt = slha.get_block_entry("yukawa", 6, 1.730000e+02); 
  mdl_ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00); 
  aS = slha.get_block_entry("sminputs", 3, 1.184000e-01); 
  mdl_Gf = slha.get_block_entry("sminputs", 2, 1.166370e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.279000e+02); 
  mdl_MH = slha.get_block_entry("mass", 25, 1.250000e+02); 
  mdl_MZ = slha.get_block_entry("mass", 23, 9.118760e+01); 
  mdl_MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  mdl_MT = slha.get_block_entry("mass", 6, 1.730000e+02); 
  mdl_MB = slha.get_block_entry("mass", 5, 4.700000e+00); 
  mdl_cy = slha.get_block_entry("bsm", 189, 1.000000e+00); 
  mdl_ctr = slha.get_block_entry("bsm", 188, 1.000000e+00); 
  mdl_a2 = slha.get_block_entry("bsm", 32, 1.000000e+00); 
  mdl_a1 = slha.get_block_entry("bsm", 31, 1.000000e+00); 
  mdl_c2 = slha.get_block_entry("bsm", 30, -1.000000e+00); 
  mdl_MZ__exp__2 = pow(mdl_MZ, 2.); 
  mdl_MZ__exp__4 = pow(mdl_MZ, 4.); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_MH__exp__2 = pow(mdl_MH, 2.); 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_cy__exp__2 = pow(mdl_cy, 2.); 
  mdl_aEW = 1./aEWM1; 
  mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = sqrt(mdl_aEW); 
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
  mdl_MW__exp__2 = pow(mdl_MW, 2.); 
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
  mdl_cw = sqrt(1. - mdl_sw2); 
  mdl_sqrt__sw2 = sqrt(mdl_sw2); 
  mdl_sw = mdl_sqrt__sw2; 
  mdl_g1 = mdl_ee/mdl_cw; 
  mdl_gw = mdl_ee/mdl_sw; 
  mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee; 
  mdl_vev__exp__2 = pow(mdl_vev, 2.); 
  mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2); 
  mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_vev; 
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev; 
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev; 
  mdl_muH = sqrt(mdl_lam * mdl_vev__exp__2); 
  mdl_ee__exp__2 = pow(mdl_ee, 2.); 
  mdl_sw__exp__2 = pow(mdl_sw, 2.); 
  mdl_cw__exp__2 = pow(mdl_cw, 2.); 
}
void Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::setIndependentCouplings()
{
  GC_13 = (mdl_ee * mdl_complexi)/(mdl_sw * mdl_sqrt__2); 
  GC_28 = -6. * mdl_complexi * mdl_lam * mdl_vev * mdl_ctr; 
  GC_29 = (mdl_ee__exp__2 * mdl_complexi * mdl_vev)/(2. * mdl_sw__exp__2); 
  GC_31 = -((mdl_complexi * mdl_yb)/mdl_sqrt__2); 
}
void Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::setDependentParameters()
{
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = pow(G, 2.); 
}
void Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::setDependentCouplings()
{
  GC_38 = -(aS * mdl_complexi)/(3. * M_PI * mdl_vev__exp__2) * mdl_a2; 
  GC_22 = -(aS * mdl_complexi)/(3. * M_PI * mdl_vev__exp__2) * (-3./2.) *
      mdl_cy__exp__2;
  GC_25 = (aS * mdl_complexi)/(3. * M_PI * mdl_vev) * (3./2.) * mdl_cy; 
  GC_34 = -(aS * mdl_complexi)/(3. * M_PI * mdl_vev__exp__2) * (-3./2.) *
      mdl_cy__exp__2;
  GC_37 = (aS * mdl_complexi)/(3. * M_PI * mdl_vev) * mdl_a1; 
  GC_36 = -mdl_c2 * (aS * mdl_complexi)/(3. * M_PI * mdl_vev__exp__2) *
      (-3./2.);
}

// Routines for printing out parameters
void Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::printIndependentParameters()
{
  LOGVRB <<  "BSM_gg_hh model parameters independent of event kinematics:" <<
      endl;
  LOGVRB << setw(20) <<  "mdl_WH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WH << endl;
  LOGVRB << setw(20) <<  "mdl_WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WW << endl;
  LOGVRB << setw(20) <<  "mdl_WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WZ << endl;
  LOGVRB << setw(20) <<  "mdl_WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WT << endl;
  LOGVRB << setw(20) <<  "mdl_ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymtau << endl;
  LOGVRB << setw(20) <<  "mdl_ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymt << endl;
  LOGVRB << setw(20) <<  "mdl_ymb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymb << endl;
  LOGVRB << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  LOGVRB << setw(20) <<  "mdl_Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gf << endl;
  LOGVRB << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  LOGVRB << setw(20) <<  "mdl_MH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MH << endl;
  LOGVRB << setw(20) <<  "mdl_MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MZ << endl;
  LOGVRB << setw(20) <<  "mdl_MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MTA << endl;
  LOGVRB << setw(20) <<  "mdl_MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MT << endl;
  LOGVRB << setw(20) <<  "mdl_MB " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MB << endl;
  LOGVRB << setw(20) <<  "mdl_cy " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cy << endl;
  LOGVRB << setw(20) <<  "mdl_ctr " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ctr << endl;
  LOGVRB << setw(20) <<  "mdl_a2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_a2 << endl;
  LOGVRB << setw(20) <<  "mdl_a1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_a1 << endl;
  LOGVRB << setw(20) <<  "mdl_c2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_c2 << endl;
  LOGVRB << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4 << endl;
  LOGVRB << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2 << endl;
  LOGVRB << setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi << endl;
  LOGVRB << setw(20) <<  "mdl_cy__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cy__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_aEW << endl;
  LOGVRB << setw(20) <<  "mdl_MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MW << endl;
  LOGVRB << setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__aEW << endl;
  LOGVRB << setw(20) <<  "mdl_ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ee << endl;
  LOGVRB << setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw2 << endl;
  LOGVRB << setw(20) <<  "mdl_cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cw << endl;
  LOGVRB << setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__sw2 << endl;
  LOGVRB << setw(20) <<  "mdl_sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw << endl;
  LOGVRB << setw(20) <<  "mdl_g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_g1 << endl;
  LOGVRB << setw(20) <<  "mdl_gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gw << endl;
  LOGVRB << setw(20) <<  "mdl_vev " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_vev << endl;
  LOGVRB << setw(20) <<  "mdl_vev__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_vev__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_lam << endl;
  LOGVRB << setw(20) <<  "mdl_yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yb << endl;
  LOGVRB << setw(20) <<  "mdl_yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yt << endl;
  LOGVRB << setw(20) <<  "mdl_ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ytau << endl;
  LOGVRB << setw(20) <<  "mdl_muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_muH << endl;
  LOGVRB << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2 << endl;
}
void Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::printIndependentCouplings()
{
  LOGVRB <<  "BSM_gg_hh model couplings independent of event kinematics:" <<
      endl;
  LOGVRB << setw(20) <<  "GC_13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_13 << endl;
  LOGVRB << setw(20) <<  "GC_28 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_28 << endl;
  LOGVRB << setw(20) <<  "GC_29 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_29 << endl;
  LOGVRB << setw(20) <<  "GC_31 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_31 << endl;
}
void Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::printDependentParameters()
{
  LOGVRB <<  "BSM_gg_hh model parameters dependent on event kinematics:" << endl;
  LOGVRB << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS << endl;
  LOGVRB << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  LOGVRB << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2 << endl;
}
void Parameters_BSM_gg_hh2bbWW_Wp2lv_Wn2jj::printDependentCouplings()
{
  LOGVRB <<  "BSM_gg_hh model couplings dependent on event kinematics:" << endl;
  LOGVRB << setw(20) <<  "GC_38 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_38 << endl;
  LOGVRB << setw(20) <<  "GC_22 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_22 << endl;
  LOGVRB << setw(20) <<  "GC_25 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_25 << endl;
  LOGVRB << setw(20) <<  "GC_34 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_34 << endl;
  LOGVRB << setw(20) <<  "GC_37 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_37 << endl;
  LOGVRB << setw(20) <<  "GC_36 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_36 << endl;
}


