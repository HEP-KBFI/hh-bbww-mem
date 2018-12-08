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
  /// error codes that can be read out by SVfitMEM class
  enum ErrorCodes {
    None            = 0x00000000,
    MatrixInversion = 0x00000001
  };

  MEMbbwwIntegrandBase(double, const std::string&, const std::string&, int);
  virtual ~MEMbbwwIntegrandBase();
  
  /// get dimension of integration region
  unsigned getIntNumDimensions() const { return intNumDimensions_; }

  /// get names of integration variables (to improve readability of debug output)
  const std::vector<std::string>& getIntVarNames() const { return intVarNames_; }

  /// get lower and upper bounds on integration region
  const double* getIntBounds_lower() const { return intIntBounds_lower_; }
  const double* getIntBounds_upper() const { return intIntBounds_upper_; }

  /// set measured momenta of charged leptons and b-jets and of missing transverse momentum
  virtual void setInputs(const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, const mem::MeasuredParticle&, double, double, const TMatrixD&);

  /// evaluate integrand for given value of integration variables x
  /// (pure virtual function, overwritten by derived classes for signal and background)
  virtual double Eval(const double* x) const = 0;

  /// pointer to parton-distribution-functions (PDF)
  /// Note: pointer is static so that same PDF instance gets used for signal and background
  static LHAPDF::PDF* pdf_;
  static std::string pdfName_;
  static bool pdfIsInitialized_;

 protected:  
  void printMadGraphMomenta() const;

  /// integration variables
  unsigned intNumDimensions_;
  std::vector<std::string> intVarNames_; 
  double* intIntBounds_lower_;
  double* intIntBounds_upper_;

  /// measured momenta of charged leptons and b-jets
  MeasuredParticle measuredChargedLeptonPlus_;
  MeasuredParticle measuredChargedLeptonMinus_;
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

  std::string madgraphFileName_;
  bool madgraphIsInitialized_;
  double* madgraphGluon1P4_;
  double* madgraphGluon2P4_;
  double* madgraphChargedLeptonPlusP4_;
  double* madgraphNeutrinoP4_;
  double* madgraphChargedLeptonMinusP4_;
  double* madgraphAntiNeutrinoP4_;
  double* madgraphBJet1P4_;
  double* madgraphBJet2P4_;
  mutable std::vector<double*> madgraphMomenta_;

  /// error code that can be passed on
  int errorCode_;

  /// verbosity level
  int verbosity_;
};

}

#endif
