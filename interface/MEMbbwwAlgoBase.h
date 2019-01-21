#ifndef hhAnalysis_bbwwMEM_MEMbbwwAlgoBase_h
#define hhAnalysis_bbwwMEM_MEMbbwwAlgoBase_h

#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"
#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorBase.h"
#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandBase.h"
#include "hhAnalysis/bbwwMEM/interface/MEMResult.h"

#include <TBenchmark.h>
#include <TMatrixD.h>
#include <TMath.h>

#include <vector> 

class MEMbbwwAlgoBase
{
 public:
  MEMbbwwAlgoBase(double, const std::string&, const std::string&, const std::string& = "", int = 0);
  virtual ~MEMbbwwAlgoBase();

  /// number of function calls for VEGAS and VAMP integration for signal hypothesis (default is 10000)
  void setMaxObjFunctionCalls_signal(unsigned maxObjFunctionCalls) 
  { 
    maxObjFunctionCalls_signal_ = maxObjFunctionCalls;
  }
  /// number of function calls for VEGAS and VAMP integration for background hypothesis (default is 100000)
  void setMaxObjFunctionCalls_background(unsigned maxObjFunctionCalls) 
  { 
    maxObjFunctionCalls_background_ = maxObjFunctionCalls;
  }

  /// set integration algorithm (either VEGAS or VAMP algorithm)
  void setIntMode(int intMode)
  {
    intMode_ = intMode;
  }
  
  /// run integration 
  virtual void integrate(const std::vector<mem::MeasuredParticle>&, double, double, const TMatrixD&) = 0;

  /// return computing time (in seconds) spent on last call to integrate method
  double getComputingTime_cpu() const { return numSeconds_cpu_; }
  double getComputingTime_real() const { return numSeconds_real_; }

  /// return number of times the MadGraph matrix element was evaluated 
  unsigned long getNumMatrixElementEvaluations_signal() const { return numMatrixElementEvaluations_signal_; }
  unsigned long getNumMatrixElementEvaluations_background() const { return numMatrixElementEvaluations_background_; }

  /// return TBenchmark object recording computing time
  const TBenchmark* getClock() const { return clock_; }

  /// static pointer to this (needed for interfacing the likelihood function calls to VEGAS and VAMP integration)
  static const mem::MEMbbwwIntegrandBase* gMEMIntegrand;

  enum { kVEGAS, kVAMP }; 

 protected:
  /// set measured momenta of charged leptons, b-jets, and (light quark) jets originating from hadronic decay of W boson
  virtual void setMeasuredParticles(const std::vector<mem::MeasuredParticle>&) = 0;

  /// set measured missing transverse momentum (MET)
  /// and MET uncertainty matrix
  virtual void setMeasuredMEt_and_Cov(double, double, const TMatrixD&);

  /// initialize integration algorithm (either VEGAS or VAMP)
  void initializeIntAlgo(unsigned maxObjFunctionCalls);

  /// run actual integration to compute compatibility of measured particles and MET with either signal or background hypothesis,
  /// and for one particular permutation of measured charged leptons and b-jets
  void runIntAlgo(mem::MEMbbwwIntegrandBase* integrand, double& prob, double& probErr);

  /// pointers to integration classes for signal and background hypotheses
  static LHAPDF::PDF* pdf_;
  std::string madgraphFileName_signal_;
  std::string madgraphFileName_background_;
  double sqrtS_;

  /// measured missing transverse momentum (MET)
  /// and MET uncertainty matrix
  double measuredMEtPx_;
  double measuredMEtPy_;
  TMatrixD measuredMEtCov_;

  /// interface to integration algorithm (either VEGAS or VAMP)
  int intMode_;
  mem::MEMIntegratorBase* intAlgo_;
  unsigned maxObjFunctionCalls_signal_;
  unsigned maxObjFunctionCalls_background_;
  double precision_;

  /// clock for measuring run-time of algorithm
  TBenchmark* clock_;
  double numSeconds_cpu_;
  double numSeconds_real_;
  unsigned long numMatrixElementEvaluations_signal_;
  unsigned long numMatrixElementEvaluations_background_;

  /// verbosity level
  int verbosity_;
};

#endif
