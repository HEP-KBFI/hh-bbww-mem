#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoBase.h"

#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorVEGAS.h"
#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorVAMP.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>

using namespace mem;
int countEval = -1;

const MEMbbwwIntegrandBase* MEMbbwwAlgoBase::gMEMIntegrand = nullptr;
LHAPDF::PDF * MEMbbwwAlgoBase::pdf_ = nullptr;

namespace
{
  double g_C(double* x, size_t, void*)
  {
    //std::cout << "<g_C>:" << std::endl;
    double retVal = MEMbbwwAlgoBase::gMEMIntegrand->Eval(x, countEval);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }

  double g_Fortran(double** x, size_t, void**)
  {
    //std::cout << "<g_Fortran>:" << std::endl;
    //int countEval;
    double retVal = MEMbbwwAlgoBase::gMEMIntegrand->Eval(*x, countEval);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }
}

MEMbbwwAlgoBase::MEMbbwwAlgoBase(double sqrtS,
				 const std::string& pdfName,
				 const std::string& madgraphFileName_signal, const std::string& madgraphFileName_background,
				 int verbosity)
  : madgraphFileName_signal_(madgraphFileName_signal)
  , madgraphFileName_background_(madgraphFileName_background)
  , sqrtS_(sqrtS)
  , intMode_(kVAMP)
  , intAlgo_(nullptr)
  , maxObjFunctionCalls_signal_(10000)
  , maxObjFunctionCalls_background_(100000)
  , precision_(1.e-3)
  , clock_(nullptr)
  , numSeconds_cpu_(-1.)
  , numSeconds_real_(-1.)
  , numMatrixElementEvaluations_signal_(0)
  , numMatrixElementEvaluations_background_(0)
  , verbosity_(verbosity)
{
  if ( !pdf_ )
  {
    pdf_ = LHAPDF::mkPDF(pdfName);
  }

  clock_ = new TBenchmark();
}

MEMbbwwAlgoBase::~MEMbbwwAlgoBase()
{
  if ( pdf_ )
  {
    delete pdf_;
    pdf_ = nullptr;
  }

  delete clock_;
}

void MEMbbwwAlgoBase::setMeasuredMEt_and_Cov(double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  measuredMEtPx_ = roundToNdigits(measuredMEtPx);
  measuredMEtPy_ = roundToNdigits(measuredMEtPy);
  measuredMEtCov_.ResizeTo(2,2);
  measuredMEtCov_[0][0] = roundToNdigits(measuredMEtCov[0][0]);
  measuredMEtCov_[1][0] = roundToNdigits(measuredMEtCov[1][0]);
  measuredMEtCov_[0][1] = roundToNdigits(measuredMEtCov[0][1]);
  measuredMEtCov_[1][1] = roundToNdigits(measuredMEtCov[1][1]);
  if ( verbosity_ >= 1 ) {
    std::cout << "MET: Px = " << measuredMEtPx_ << ", Py = " << measuredMEtPy_ << std::endl;
    std::cout << "covMET:" << std::endl;
    measuredMEtCov_.Print();
    TMatrixDSym metCov_sym(2);
    metCov_sym(0,0) = measuredMEtCov[0][0];
    metCov_sym(0,1) = measuredMEtCov[0][1];
    metCov_sym(1,0) = measuredMEtCov[1][0];
    metCov_sym(1,1) = measuredMEtCov[1][1];
    TMatrixD EigenVectors(2,2);
    EigenVectors = TMatrixDSymEigen(metCov_sym).GetEigenVectors();
    std::cout << "Eigenvectors =  { " << EigenVectors(0,0) << ", " << EigenVectors(1,0) << " (phi = " << TMath::ATan2(EigenVectors(1,0), EigenVectors(0,0)) << ") },"
	      << " { " << EigenVectors(0,1) << ", " << EigenVectors(1,1) << " (phi = " << TMath::ATan2(EigenVectors(1,1), EigenVectors(0,1)) << ") }" << std::endl;
    TVectorD EigenValues(2);
    EigenValues = TMatrixDSymEigen(metCov_sym).GetEigenValues();
    EigenValues(0) = TMath::Sqrt(EigenValues(0));
    EigenValues(1) = TMath::Sqrt(EigenValues(1));
    std::cout << "Eigenvalues = " << EigenValues(0) << ", " << EigenValues(1) << std::endl;
  }
}

void MEMbbwwAlgoBase::initializeIntAlgo(unsigned maxObjFunctionCalls)
{
  if ( intMode_ == kVEGAS ) {
    unsigned numCallsGridOpt = TMath::Nint(0.20*maxObjFunctionCalls);
    unsigned numCallsIntEval = TMath::Nint(0.80*maxObjFunctionCalls);
    intAlgo_ = new MEMIntegratorVEGAS(
      numCallsGridOpt, numCallsIntEval,
      2., 1);
  } else if ( intMode_ == kVAMP ) {
    unsigned numCallsGridOpt = TMath::Nint(0.20*maxObjFunctionCalls);
    unsigned numCallsIntEval = TMath::Nint(0.80*maxObjFunctionCalls);
    intAlgo_ = new MEMIntegratorVAMP(
      numCallsGridOpt, numCallsIntEval);
  } else {
    std::cerr << "<MEMbbwwAlgoBase::initializeIntAlgo>: Invalid configuration parameter 'intMode' = " << intMode_ << " --> ABORTING !!\n";
    assert(0);
  }
}

void MEMbbwwAlgoBase::runIntAlgo(mem::MEMbbwwIntegrandBase* integrand, double& prob, double& probErr)
{
  unsigned numDimensions = integrand->getIntNumDimensions();
  assert(numDimensions >= 1);
  const double* xl = integrand->getIntBounds_lower();
  const double* xu = integrand->getIntBounds_upper();
  if ( verbosity_ >= 1 ) {
    const std::vector<std::string>& intVarNames = integrand->getIntVarNames();
    assert(intVarNames.size() == numDimensions);
    for ( unsigned idxDimension = 0; idxDimension < numDimensions; ++idxDimension ) {
      std::cout << " intVariable #" << idxDimension << " (" << intVarNames[idxDimension] << "): xl = " << xl[idxDimension] << ", xu = " << xu[idxDimension];
      std::cout << std::endl;
    }
  }

  prob = 0.;
  probErr = 0.;
  if ( intMode_ == kVEGAS ) {
    intAlgo_->integrate(&g_C, xl, xu, numDimensions, prob, probErr);
  } else if ( intMode_ == kVAMP ) {
    intAlgo_->integrate(&g_Fortran, xl, xu, numDimensions, prob, probErr);
  } else assert(0);
}
