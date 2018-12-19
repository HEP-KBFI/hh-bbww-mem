#ifndef hhAnalysis_bbwwMEM_MEMbbwwIntegrandDilepton_background_h
#define hhAnalysis_bbwwMEM_MEMbbwwIntegrandDilepton_background_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandDilepton.h"
#include "hhAnalysis/bbwwMEM/interface/mg5/me/mg5_ttbar2WbWb_WW2lvlv.h"

namespace mem
{

class MEMbbwwIntegrandDilepton_background : public MEMbbwwIntegrandDilepton
{
 public:
  MEMbbwwIntegrandDilepton_background(double, const std::string&, int);
  ~MEMbbwwIntegrandDilepton_background();

  /// set measured momenta of charged leptons and b-jets and of missing transverse momentum
  void setInputs(const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, 
		 double, double, const TMatrixD&);

  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x) const;

 protected:  
  /// leading order (LO) matrix element obtained from MadGraph
  mutable mg5_sm_ttbar2WbWb_WW2lvlv me_madgraph_;
};

}

#endif
