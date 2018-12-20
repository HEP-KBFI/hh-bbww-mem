//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "hhAnalysis/bbwwMEM/interface/mg5/parameters/Parameters_sm_ttbar2WbWb_WW2lvlv.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"

#include <iomanip>

// Initialize static instance
Parameters_sm_ttbar2WbWb_WW2lvlv * Parameters_sm_ttbar2WbWb_WW2lvlv::instance = nullptr;

// Function to get static instance - only one instance per program
Parameters_sm_ttbar2WbWb_WW2lvlv * Parameters_sm_ttbar2WbWb_WW2lvlv::getInstance()
{
  if (! instance)
    instance = new Parameters_sm_ttbar2WbWb_WW2lvlv();

  return instance; 
}

void Parameters_sm_ttbar2WbWb_WW2lvlv::setIndependentParameters(SLHAReader& slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  mdl_WH = slha.get_block_entry("decay", 25, 6.382339e-03); 
  mdl_WW = slha.get_block_entry("decay", 24, 2.047600e+00); 
  mdl_WZ = slha.get_block_entry("decay", 23, 2.441404e+00); 
  mdl_WT = slha.get_block_entry("decay", 6, 1.491500e+00); 
  mdl_ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  mdl_ymt = slha.get_block_entry("yukawa", 6, 1.730000e+02); 
  mdl_ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00); 
  aS = slha.get_block_entry("sminputs", 3, 1.180000e-01); 
  mdl_Gf = slha.get_block_entry("sminputs", 2, 1.166390e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.325070e+02); 
  mdl_MH = slha.get_block_entry("mass", 25, 1.250000e+02); 
  mdl_MZ = slha.get_block_entry("mass", 23, 9.118800e+01); 
  mdl_MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  mdl_MT = slha.get_block_entry("mass", 6, 1.730000e+02); 
  mdl_MB = slha.get_block_entry("mass", 5, 4.700000e+00); 
  mdl_conjg__CKM3x3 = 1.; 
  mdl_CKM3x3 = 1.; 
  mdl_conjg__CKM1x1 = 1.; 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_MZ__exp__2 = pow(mdl_MZ, 2.); 
  mdl_MZ__exp__4 = pow(mdl_MZ, 4.); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_MH__exp__2 = pow(mdl_MH, 2.); 
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
  mdl_I1x33 = mdl_yb * mdl_conjg__CKM3x3; 
  mdl_I2x33 = mdl_yt * mdl_conjg__CKM3x3; 
  mdl_I3x33 = mdl_CKM3x3 * mdl_yt; 
  mdl_I4x33 = mdl_CKM3x3 * mdl_yb; 
  mdl_ee__exp__2 = pow(mdl_ee, 2.); 
  mdl_sw__exp__2 = pow(mdl_sw, 2.); 
  mdl_cw__exp__2 = pow(mdl_cw, 2.); 
}
void Parameters_sm_ttbar2WbWb_WW2lvlv::setIndependentCouplings()
{
  GC_100 = (mdl_ee * mdl_complexi * mdl_conjg__CKM1x1)/(mdl_sw * mdl_sqrt__2); 
}
void Parameters_sm_ttbar2WbWb_WW2lvlv::setDependentParameters()
{
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = pow(G, 2.); 
}
void Parameters_sm_ttbar2WbWb_WW2lvlv::setDependentCouplings()
{
  GC_11 = mdl_complexi * G; 
  GC_10 = -G; 
}

// Routines for printing out parameters
void Parameters_sm_ttbar2WbWb_WW2lvlv::printIndependentParameters()
{
  LOGVRB <<  "sm model parameters independent of event kinematics:" << endl;
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
  LOGVRB << setw(20) <<  "mdl_conjg__CKM3x3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM3x3 << endl;
  LOGVRB << setw(20) <<  "mdl_CKM3x3 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM3x3 << endl;
  LOGVRB << setw(20) <<  "mdl_conjg__CKM1x1 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM1x1 << endl;
  LOGVRB << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi << endl;
  LOGVRB << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4 << endl;
  LOGVRB << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2 << endl;
  LOGVRB << setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__2 << endl;
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
  LOGVRB << setw(20) <<  "mdl_I1x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1x33 << endl;
  LOGVRB << setw(20) <<  "mdl_I2x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2x33 << endl;
  LOGVRB << setw(20) <<  "mdl_I3x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3x33 << endl;
  LOGVRB << setw(20) <<  "mdl_I4x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4x33 << endl;
  LOGVRB << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2 << endl;
  LOGVRB << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2 << endl;
}
void Parameters_sm_ttbar2WbWb_WW2lvlv::printIndependentCouplings()
{
  LOGVRB <<  "sm model couplings independent of event kinematics:" << endl;
  LOGVRB << setw(20) <<  "GC_100 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_100 << endl;
}
void Parameters_sm_ttbar2WbWb_WW2lvlv::printDependentParameters()
{
  LOGVRB <<  "sm model parameters dependent on event kinematics:" << endl;
  LOGVRB << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS << endl;
  LOGVRB << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  LOGVRB << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2 << endl;
}
void Parameters_sm_ttbar2WbWb_WW2lvlv::printDependentCouplings()
{
  LOGVRB <<  "sm model couplings dependent on event kinematics:" << endl;
  LOGVRB << setw(20) <<  "GC_11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_11 << endl;
  LOGVRB << setw(20) <<  "GC_10 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_10 << endl;
}


