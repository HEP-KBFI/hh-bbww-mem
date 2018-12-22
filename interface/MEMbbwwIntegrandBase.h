#ifndef hhAnalysis_bbwwMEM_MEMbbwwIntegrandBase_h
#define hhAnalysis_bbwwMEM_MEMbbwwIntegrandBase_h

#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"
#include "hhAnalysis/bbwwMEM/interface/BJetTF.h"
#include "hhAnalysis/bbwwMEM/interface/HadRecoilTF.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
// ignore "-Wmaybe-uninitialized" and "-Wshadow" gcc errors if compiled with -Og -g3 -ggdb3
#if defined(__OPTIMIZE__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <LHAPDF/LHAPDF.h> // LHAPDF::PDF

#pragma GCC diagnostic pop

#include <TMatrixD.h>

#include <vector>
#include <string>

namespace mem
{
  
class MEMbbwwIntegrandBase
{
 public:
  /// error codes
  enum ErrorCodes {
    None            = 0x00000000,
    MatrixInversion = 0x00000001
  };

  MEMbbwwIntegrandBase(double, const std::string&, int);
  virtual ~MEMbbwwIntegrandBase();
  
  /// get dimension of integration region
  unsigned getIntNumDimensions() const { return intNumDimensions_; }

  /// get names of integration variables (to improve readability of debug output)
  const std::vector<std::string>& getIntVarNames() const { return intVarNames_; }

  /// get lower and upper bounds on integration region
  const double* getIntBounds_lower() const { return intIntBounds_lower_; }
  const double* getIntBounds_upper() const { return intIntBounds_upper_; }

  /// reset counter for number of times the MadGraph matrix element was evaluated 
  void resetNumMatrixElementEvaluations() { numMatrixElementEvaluations_ = 0; }

  /// evaluate integrand for given value of integration variables x
  /// (pure virtual function, overwritten by derived classes for signal and background)
  virtual double Eval(const double* x) const = 0;

  /// get number of times the MadGraph matrix element was evaluated 
  unsigned long getNumMatrixElementEvaluations() const { return numMatrixElementEvaluations_; }

  /// set parton distribution function (PDF)
  void setPDF(LHAPDF::PDF* pdf);

 protected:  
  /// pointer to parton-distribution-function (PDF)
  LHAPDF::PDF* pdf_;

  /// print four-vectors passed to MadGraph 
  void printMadGraphMomenta(const std::vector<double*>& ) const;
  
  /// integration variables
  unsigned intNumDimensions_;
  std::vector<std::string> intVarNames_; 
  double* intIntBounds_lower_;
  double* intIntBounds_upper_;

  /// measured momenta of charged leptons and b-jets
  MeasuredParticle measuredBJet1_; 
  MeasuredParticle measuredBJet2_;

  /// measured missing transverse momentum (MET)
  /// and MET uncertainty matrix
  double measuredMEtPx_;
  double measuredMEtPy_;
  TMatrixD measuredMEtCov_;

  /// measured momentum of hadronic recoil
  double measuredHadRecoilPx_;
  double measuredHadRecoilPy_;

  /// transfer functions for b-jets and MET
  BJetTF* bjet1TF_;
  BJetTF* bjet2TF_;
  HadRecoilTF* hadRecoilTF_;

  /// center-of-mass energy
  double sqrtS_;

  /// normalization factor for probability densitity w_1(y|yhat)
  double normFactor_;

  /// file with MadGraph model parameters and four-vectors used to evaluate MadGraph matrix element
  std::string madgraphFileName_;
  bool madgraphIsInitialized_;
  double* madgraphGluon1P4_;
  double* madgraphGluon2P4_;
  double* madgraphBJet1P4_;
  double* madgraphBJet2P4_;
  mutable unsigned long numMatrixElementEvaluations_;

  /// error code that can be passed on
  int errorCode_;

  /// verbosity level
  int verbosity_;
};

}

#endif
