#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandDilepton.h"

#include <TMath.h>

using namespace mem;

MEMbbwwIntegrandDilepton::MEMbbwwIntegrandDilepton(double sqrtS, const std::string& madgraphFileName, int verbosity)
  : MEMbbwwIntegrandBase(sqrtS, madgraphFileName, verbosity)
  , madgraphChargedLeptonPlusP4_(nullptr)
  , madgraphNeutrinoP4_(nullptr)
  , madgraphChargedLeptonMinusP4_(nullptr)
  , madgraphAntiNeutrinoP4_(nullptr)
{
  madgraphChargedLeptonPlusP4_ = new double[4];
  madgraphNeutrinoP4_ = new double[4];
  madgraphChargedLeptonMinusP4_ = new double[4];
  madgraphAntiNeutrinoP4_ = new double[4];
}

MEMbbwwIntegrandDilepton::~MEMbbwwIntegrandDilepton()
{
  if ( verbosity_ >= 1 ) 
  {
    std::cout << "<MEMbbwwIntegrandDilepton::~MEMbbwwIntegrandDilepton>:" << std::endl;
  }
  
  delete [] madgraphChargedLeptonPlusP4_;
  delete [] madgraphNeutrinoP4_;
  delete [] madgraphChargedLeptonMinusP4_;
  delete [] madgraphAntiNeutrinoP4_;
}

void 
MEMbbwwIntegrandDilepton::setInputs(const MeasuredParticle& measuredChargedLeptonPlus, const MeasuredParticle& measuredChargedLeptonMinus,
				    const MeasuredParticle& measuredBJet1, const MeasuredParticle& measuredBJet2,
				    double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<MEMbbwwIntegrandDilepton::setInputs>:" << std::endl;
  }
  // reset 'MatrixInversion' error code
  errorCode_ &= (errorCode_ ^ MatrixInversion);
  measuredChargedLeptonPlus_ = measuredChargedLeptonPlus;
  measuredChargedLeptonMinus_ = measuredChargedLeptonMinus;
  measuredBJet1_ = measuredBJet1;
  measuredBJet2_ = measuredBJet2;
  measuredMEtPx_ = measuredMEtPx;
  measuredMEtPy_ = measuredMEtPy;
  measuredMEtCov_.ResizeTo(2,2);
  measuredMEtCov_ = measuredMEtCov;
  measuredHadRecoilPx_ = -(measuredChargedLeptonPlus_.px() + measuredChargedLeptonMinus_.px() + measuredBJet1_.px() + measuredBJet2_.px() + measuredMEtPx_);
  measuredHadRecoilPy_ = -(measuredChargedLeptonPlus_.py() + measuredChargedLeptonMinus_.py() + measuredBJet1_.py() + measuredBJet2_.py() + measuredMEtPy_);
  // set measured momenta of b-jets and of missing transverse momentum
  // in transfer function (TF) objects
  bjet1TF_->setInputs(measuredBJet1_.p4());
  bjet2TF_->setInputs(measuredBJet2_.p4());
  hadRecoilTF_->setInputs(measuredHadRecoilPx_, measuredHadRecoilPy_, measuredMEtCov_);
  if ( hadRecoilTF_->getErrorCode() ) {
    errorCode_ |= MatrixInversion;
  }
}
