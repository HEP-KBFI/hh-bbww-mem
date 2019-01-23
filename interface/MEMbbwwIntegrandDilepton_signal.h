#ifndef hhAnalysis_bbwwMEM_MEMbbwwIntegrandDilepton_signal_h
#define hhAnalysis_bbwwMEM_MEMbbwwIntegrandDilepton_signal_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandDilepton.h"
#include "hhAnalysis/bbwwMEM/interface/mg5/me/mg5_gg_hh2bbWW_WW2lvlv.h"

namespace mem
{

class MEMbbwwIntegrandDilepton_signal : public MEMbbwwIntegrandDilepton
{
 public:
  MEMbbwwIntegrandDilepton_signal(double, const std::string&, int);
  ~MEMbbwwIntegrandDilepton_signal();

  /// fix (flag=true) mass of charged lepton plus neutrino originating from the decay of the "on-shell" W boson to mW,
  /// or allow the mass to vary during the integration (flag=false)
  void applyOnshellWmassConstraint(bool flag);

  /// set measured momenta of charged leptons and b-jets and of missing transverse momentum
  void setInputs(const mem::MeasuredParticle*, const mem::MeasuredParticle*, 
		 const mem::MeasuredParticle*, const mem::MeasuredParticle*, 
		 double, double, const TMatrixD&);

  /// switch between associations of lepton+ and lepton- to on-shell and off-shell W bosons
  //enum { kPermutationUndefined, kOnshellChargedLeptonPlus, kOnshellChargedLeptonMinus };
  void setOnshellChargedLepton(int chargedLeptonPermutation);

  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x) const;

 protected:  
  /// initialize integration variables (grid)
  void initializeIntVars();

  /// leading order (LO) matrix element obtained from MadGraph
  mutable mg5_BSM_gg_hh2bbWW_WW2lvlv me_madgraph_;

  /// flag to either fix (applyOnshellWmassConstraint=true) mass of charged lepton plus neutrino originating from the decay of the "on-shell" W boson to mW,
  /// or allow the mass to vary during the integration (applyOnshellWmassConstraint=false)
  bool applyOnshellWmassConstraint_;

  /// index of integration variables for theta and phi of either b-jet1 or b-jet2 in case one b-jet is "missing", i.e. not reconstructed
  int offsetBJet1Theta_;
  int offsetBJet1Phi_;
  int offsetBJet2Theta_;
  int offsetBJet2Phi_;

  /// flag to switch between associations of lepton+ and lepton- to on-shell and off-shell W bosons
  int chargedLeptonPermutation_;
};

}

#endif
