#ifndef hhAnalysis_bbwwMEM_MEMbbwwIntegrandSingleLepton_background_h
#define hhAnalysis_bbwwMEM_MEMbbwwIntegrandSingleLepton_background_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandSingleLepton.h"
#include "hhAnalysis/bbwwMEM/interface/mg5/me/mg5_ttbar2WbWb_Wp2jj_Wn2vl.h"
#include "hhAnalysis/bbwwMEM/interface/mg5/me/mg5_ttbar2WbWb_Wp2lv_Wn2jj.h"

namespace mem
{

class MEMbbwwIntegrandSingleLepton_background : public MEMbbwwIntegrandSingleLepton
{
 public:
  MEMbbwwIntegrandSingleLepton_background(double, const std::string&, int);
  ~MEMbbwwIntegrandSingleLepton_background();

  /// set measured momenta of charged lepton, jets from W->jj decay, and b-jets and of missing transverse momentum
  void setInputs(const mem::MeasuredParticle*,
		 const mem::MeasuredParticle*, const mem::MeasuredParticle*,
		 const mem::MeasuredParticle*, const mem::MeasuredParticle*,
		 double, double, const TMatrixD&);

  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x, int & countEval) const;

 protected:
  /// initialize integration variables (grid)
  void initializeIntVars();

  /// leading order (LO) matrix element obtained from MadGraph
  ///
  /// Note: separate matrix elements are used for events with leptons of positive and events with leptons of negative charge
  mutable mg5_sm_ttbar2WbWb_Wp2lv_Wn2jj me_madgraph_chargedLeptonPlus_;
  mutable mg5_sm_ttbar2WbWb_Wp2jj_Wn2vl me_madgraph_chargedLeptonMinus_;

  /// index of integration variables for theta and phi of either jet1 or jet2 from W->jj decay in case one jet from W->jj decay is "missing", i.e. not reconstructed
  int offsetHadWJet1Theta_;
  int offsetHadWJet1Phi_;
  int offsetHadWJet2Theta_;
  int offsetHadWJet2Phi_;
  std::string madgraphFileName_;

  /// index of integration variables for theta and phi of either b-jet1 or b-jet2 in case one b-jet is "missing", i.e. not reconstructed
  int offsetBJet1Theta_;
  int offsetBJet1Phi_;
  int offsetBJet2Theta_;
  int offsetBJet2Phi_;
};

}

#endif
