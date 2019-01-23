#ifndef hhAnalysis_bbwwMEM_MEMbbwwIntegrandSingleLepton_h
#define hhAnalysis_bbwwMEM_MEMbbwwIntegrandSingleLepton_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandBase.h"
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"
#include "hhAnalysis/bbwwMEM/interface/HadWJetTF.h"

#include <TMatrixD.h>

#include <vector>
#include <string>

namespace mem
{
  
class MEMbbwwIntegrandSingleLepton : public MEMbbwwIntegrandBase
{
 public:
  MEMbbwwIntegrandSingleLepton(double, const std::string&, int);
  virtual ~MEMbbwwIntegrandSingleLepton();
  
  /// set measured momenta of charged leptons and b-jets and of missing transverse momentum
  virtual void setInputs(const mem::MeasuredParticle*, 
			 const mem::MeasuredParticle*, const mem::MeasuredParticle*, 
			 const mem::MeasuredParticle*, const mem::MeasuredParticle*,
			 double, double, const TMatrixD&);

 protected:  
  /// measured momenta of charged lepton and of the two jets from the W->jj decay
  const MeasuredParticle* measuredChargedLepton_;
  const MeasuredParticle* measuredHadWJet1_;
  const MeasuredParticle* measuredHadWJet2_;

  /// transfer functions for jets from W->jj decay
  HadWJetTF* hadWJet1TF_;
  HadWJetTF* hadWJet2TF_;

  /// four-vectors used to evaluate MadGraph matrix element
  double* madgraphChargedLeptonP4_;
  double* madgraphNeutrinoP4_;
  double* madgraphHadWJet1P4_;
  double* madgraphHadWJet2P4_;
  mutable std::vector<double*> madgraphMomenta_chargedLeptonPlus_;
  mutable std::vector<double*> madgraphMomenta_chargedLeptonMinus_;
};

}

#endif
