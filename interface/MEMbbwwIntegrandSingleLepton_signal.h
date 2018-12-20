#ifndef hhAnalysis_bbwwMEM_MEMbbwwIntegrandSingleLepton_signal_h
#define hhAnalysis_bbwwMEM_MEMbbwwIntegrandSingleLepton_signal_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandSingleLepton.h"
//#include "hhAnalysis/bbwwMEM/interface/mg5/me/mg5_gg_hh2bbWW_XXX.h" <-- TO BE UPDATED !!

namespace mem
{

class MEMbbwwIntegrandSingleLepton_signal : public MEMbbwwIntegrandSingleLepton
{
 public:
  MEMbbwwIntegrandSingleLepton_signal(double, const std::string&, int);
  ~MEMbbwwIntegrandSingleLepton_signal();

  /// fix (flag=true) mass of charged lepton plus neutrino originating from the decay of the "on-shell" W boson to mW,
  /// or allow the mass to vary during the integration (flag=false)
  void applyOnshellWmassConstraint(bool flag);

  // set measured momenta of charged lepton, jets from W->jj decay, and b-jets and of missing transverse momentum
  void setInputs(const mem::MeasuredParticle&, const mem::MeasuredParticle&, 
		 const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, 
		 double, double, const TMatrixD&);

  /// switch between associations of charged lepton and (light quark) jet pair to on-shell and off-shell W bosons
  void setOnshellChargedLepton(int chargedLeptonPermutation);

  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x) const;

 protected:  
  /// initialize integration variables (grid)
  void initializeIntVars();

  /// leading order (LO) matrix element obtained from MadGraph
//mutable mg5_BSM_gg_hh2bbWW_WW2lvlv me_madgraph_;

  /// flag to either fix (applyOnshellWmassConstraint=true) mass of charged lepton plus neutrino originating from the decay of the "on-shell" W boson to mW,
  /// or allow the mass to vary during the integration (applyOnshellWmassConstraint=false)
  bool applyOnshellWmassConstraint_;

  /// flag to switch between associations of charged lepton and (light quark) jet pair to on-shell and off-shell W bosons
  int chargedLeptonPermutation_;
};

}

#endif
