#ifndef hhAnalysis_bbwwMEM_MEMbbwwIntegrandDilepton_h
#define hhAnalysis_bbwwMEM_MEMbbwwIntegrandDilepton_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandBase.h"
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"

#include <TMatrixD.h>

#include <vector>
#include <string>

namespace mem
{
  
class MEMbbwwIntegrandDilepton : public MEMbbwwIntegrandBase
{
 public:
  MEMbbwwIntegrandDilepton(double, const std::string&, int);
  virtual ~MEMbbwwIntegrandDilepton();
  
  /// set measured momenta of charged leptons and b-jets and of missing transverse momentum
  virtual void setInputs(const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, 
			 double, double, const TMatrixD&);

 protected:  
  /// measured momenta of charged leptons
  MeasuredParticle measuredChargedLeptonPlus_;
  MeasuredParticle measuredChargedLeptonMinus_;  

  /// four-vectors used to evaluate MadGraph matrix element
  double* madgraphChargedLeptonPlusP4_;
  double* madgraphNeutrinoP4_;
  double* madgraphChargedLeptonMinusP4_;
  double* madgraphAntiNeutrinoP4_;
  mutable std::vector<double*> madgraphMomenta_;
};

}

#endif
