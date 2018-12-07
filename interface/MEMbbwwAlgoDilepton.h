#ifndef hhAnalysis_bbwwMEM_MEMbbwwAlgoDilepton_h
#define hhAnalysis_bbwwMEM_MEMbbwwAlgoDilepton_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandDilepton_signal.h"
#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandDilepton_background.h"
#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorBase.h"
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"

#include <TBenchmark.h>
#include <TMatrixD.h>
#include <TMath.h>

#include <vector> 

class MEMbbwwAlgoDilepton
{
 public:
  MEMbbwwAlgoDilepton(double, const std::string&, const std::string& = "", int = 0); 
  ~MEMbbwwAlgoDilepton();

  /// number of function calls for VEGAS and VAMP integration (default is 100000)
  void setMaxObjFunctionCalls(unsigned maxObjFunctionCalls) 
  { 
    maxObjFunctionCalls_ = maxObjFunctionCalls;
  }

  /// set integration algorithm (either VEGAS or VAMP algorithm)
  void setIntMode(int intMode)
  {
    intMode_ = intMode;
  }
  
  /// run integration 
  void integrate(const std::vector<mem::MeasuredParticle>&, double, double, const TMatrixD&);

  /// return probabilities for signal and background hypotheses
  struct resultType
  {
    double prob_signal_;
    double probErr_signal_;
    double prob_background_;
    double probErr_background_;
  };
  resultType getResult() const { return result_; }

  /// return computing time (in seconds) spent on last call to integrate method
  double getComputingTime_cpu() const { return numSeconds_cpu_; }
  double getComputingTime_real() const { return numSeconds_real_; }

  /// static pointer to this (needed for interfacing the likelihood function calls to VEGAS and VAMP integration)
  static const mem::MEMbbwwIntegrandBase* gMEMIntegrand;

  enum { kVEGAS, kVAMP }; 

 protected:
  /// initialize integration algorithm (either VEGAS or VAMP)
  void initializeIntAlgo();

  /// run actual integration to compute compatibility of measured particles and MET with either signal or background hypothesis,
  /// and for one particular permutation of measured charged leptons and b-jets
  void runIntAlgo(mem::MEMbbwwIntegrandBase* integrand, double& prob, double& probErr);

  /// pointers to integration classes for signal and background hypotheses
  mem::MEMbbwwIntegrandDilepton_signal* integrand_signal_;
  mem::MEMbbwwIntegrandDilepton_background* integrand_background_;
  double sqrtS_;

  /// measured momenta of charged leptons and b-jets
  std::vector<mem::MeasuredParticle> measuredParticles_;
  const mem::MeasuredParticle* measuredChargedLeptonPlus_;
  const mem::MeasuredParticle* measuredChargedLeptonMinus_;
  const mem::MeasuredParticle* measuredLeadingBJet_; 
  const mem::MeasuredParticle* measuredSubleadingBJet_;

  /// measured missing transverse momentum (MET)
  /// and MET uncertainty matrix
  double measuredMEtPx_;
  double measuredMEtPy_;
  TMatrixD measuredMEtCov_;

  /// interface to integration algorithm (either VEGAS or VAMP)
  int intMode_;
  mem::MEMIntegratorBase* intAlgo_;
  unsigned maxObjFunctionCalls_;
  double precision_;

  /// result of integration (probabilities for signal and background hypotheses)
  resultType result_;

  /// clock for measuring run-time of algorithm
  TBenchmark* clock_;
  double numSeconds_cpu_;
  double numSeconds_real_;

  /// verbosity level
  int verbosity_;
};

#endif
