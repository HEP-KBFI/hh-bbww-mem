#ifndef hhAnalysis_bbwwMEM_memAuxFunctions_h
#define hhAnalysis_bbwwMEM_memAuxFunctions_h

#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

#include <string>
#include <vector>

namespace mem
{
  inline double square(double x)
  {
    return x*x;
  }

  inline double cube(double x)
  {
    return x*x*x;
  }

  inline double fourth(double x)
  {
    return x*x*x*x;
  }
  
  inline double fifth(double x)
  {
    return x*x*x*x*x;
  }

  inline double sixth(double x)
  {
    return x*x*x*x*x*x;
  }

  inline double seventh(double x)
  {
    return x*x*x*x*x*x*x;
  }

  inline double eigth(double x)
  {
    return x*x*x*x*x*x*x*x;
  }

  //-----------------------------------------------------------------------------
  // define masses, widths and lifetimes of particles
  // relevant for computing values of likelihood functions in MEM
  //
  // NOTE: the values are taken from
  //        K.A. Olive et al. (Particle Data Group),
  //        Chin. Phys. C 38, 090001 (2014)
  //
  const double electronMass = 0.51100e-3; // GeV
  const double electronMass2 = electronMass*electronMass;
  const double muonMass = 0.10566; // GeV
  const double muonMass2 = muonMass*muonMass; 

  const double bottomQuarkMass = 4.18; // GeV
  const double bottomQuarkMass2 = bottomQuarkMass*bottomQuarkMass; 

  const double topQuarkMass = 173.0; // GeV
  const double topQuarkMass2 = topQuarkMass*topQuarkMass; 
  const double topQuarkWidth = 1.41; // GeV
  const double topQuarkWidth2 = topQuarkWidth*topQuarkWidth; 

  const double wBosonMass = 80.379; // GeV
  const double wBosonMass2 = wBosonMass*wBosonMass; 
  const double wBosonWidth = 2.085; // GeV
  const double wBosonWidth2 = wBosonWidth*wBosonWidth;

  //const double higgsBosonMass = 125.18; // GeV
  const double higgsBosonMass = 125.; // GeV (mH=125 GeV used by MadGraph)
  const double higgsBosonMass2 = higgsBosonMass*higgsBosonMass; 
  const double higgsBosonWidth = 4.07e-3; // GeV (Standard Model prediction, not measured value)
  const double higgsBosonWidth2 = higgsBosonWidth*higgsBosonWidth;

  const double hbar_c = 0.1973; // GeV fm
  //-----------------------------------------------------------------------------

  /// define flags used to switch between associations of lepton+ and lepton- to on-shell and off-shell W bosons
  /// in case of signal hypothesis in dilepton channel
  enum { kPermutationUndefined2L, kOnshellChargedLeptonPlus, kOnshellChargedLeptonMinus };

  /// define flags used to switch between associations of charged lepton and (light quark) jet pair to on-shell and off-shell W bosons
  /// in case of signal hypothesis in single lepton channel
  enum { kPermutationUndefined1L, kOnshellChargedLepton, kOffshellChargedLepton };

  /**
     \typedef mem::Vector
     \brief   spacial momentum vector (equivalent to reco::Candidate::Vector)
  */
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;
  /**
     \typedef mem::LorentzVector
     \brief   lorentz vector (equivalent to reco::Candidate::LorentzVector)
  */
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

  double roundToNdigits(double, int = 3);

  void printLorentzVector(const std::string&, const LorentzVector&);
  void printLorentzVector(const std::string&, const LorentzVector&, const LorentzVector&);
  void printLorentzVector_NA(const std::string&, const LorentzVector&);
  void printVDouble(const std::string&, const double*, unsigned); 

  std::vector<double> compBJetEn_Hbb(const LorentzVector&, double, double);
  double compNuEn_Wlnu(const LorentzVector&, double, double, double = wBosonMass2);
  double compNuStarEn_Hww(const LorentzVector&, double, double);
  std::vector<double> compBJetEn_top(const LorentzVector&, double, double);
  double compHadWJet1En_Wjj(const LorentzVector&, double, double, double = wBosonMass2);
  double compHadWJet2En_Wjj(const LorentzVector&, double, double, double = wBosonMass2);
  double compHadWJet1En_Hww(const LorentzVector&, double, double);
  double compHadWJet2En_Hww(const LorentzVector&, double, double);

  LorentzVector findBJetEn_solution_Hbb(const std::vector<double>&, double, double, const LorentzVector&, bool&);
  LorentzVector findBJetEn_solution_top(const std::vector<double>&, double, double, const LorentzVector&, bool&);

  LorentzVector buildLorentzVector(double, double, double);
  LorentzVector buildLorentzVector(double, double, double, double);

  double compJacobiFactor_Hbb(const LorentzVector&, const LorentzVector&);
  double compJacobiFactor_Wlnu(const LorentzVector&, const LorentzVector&);
  double compJacobiFactor_Hww(const LorentzVector&, const LorentzVector&);
  double compJacobiFactor_top(const LorentzVector&, const LorentzVector&);
  double compJacobiFactor_Wjj(const LorentzVector&, const LorentzVector&);
}

#endif
