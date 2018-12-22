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
  void setInputs(const mem::MeasuredParticle&, const mem::MeasuredParticle&, 
		 const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, 
		 double, double, const TMatrixD&);

  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x) const;

 protected:  
  /// leading order (LO) matrix element obtained from MadGraph
  ///
  /// Note: separate matrix elements are used for events with leptons of positive and events with leptons of negative charge
  mutable mg5_sm_ttbar2WbWb_Wp2lv_Wn2jj me_madgraph_chargedLeptonPlus_;
  mutable mg5_sm_ttbar2WbWb_Wp2jj_Wn2vl me_madgraph_chargedLeptonMinus_;
};

}

#endif
