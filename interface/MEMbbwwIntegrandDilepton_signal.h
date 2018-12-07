#ifndef hhAnalysis_bbwwMEM_MEMbbwwIntegrandDilepton_signal_h
#define hhAnalysis_bbwwMEM_MEMbbwwIntegrandDilepton_signal_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandBase.h"
#include "hhAnalysis/bbwwMEM/interface/mg5/me/mg5_gg_hh2bbWW_WW2lvlv.h"

namespace mem
{

class MEMbbwwIntegrandDilepton_signal : public MEMbbwwIntegrandBase
{
 public:
  MEMbbwwIntegrandDilepton_signal(double, const std::string&, const std::string&, int);
  ~MEMbbwwIntegrandDilepton_signal();

  /// set measured momenta of charged leptons and b-jets and of missing transverse momentum
  void setInputs(const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, double, double, const TMatrixD&);

  /// switch between associations of lepton+ and lepton- to on-shell and off-shell W bosons
  enum { kPermutationUndefined, kOnshellChargedLeptonPlus, kOnshellChargedLeptonMinus };
  void setOnshellChargedLepton(int chargedLeptonPermutation);

  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x) const;

 protected:  
  /// leading order (LO) matrix element obtained from MadGraph
  mutable mg5_BSM_gg_hh2bbWW_WW2lvlv me_madgraph_;

  int chargedLeptonPermutation_;
};

}

#endif
