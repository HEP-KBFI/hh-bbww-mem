#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandBase.h"

#include <TMath.h>

using namespace mem;

MEMbbwwIntegrandBase::MEMbbwwIntegrandBase(double sqrtS, const std::string& madgraphFileName, int verbosity)
  : bjet1TF_(nullptr)
  , bjet2TF_(nullptr)
  , hadRecoilTF_(nullptr)
  , sqrtS_(sqrtS)
  , normFactor_(-1.)
  , madgraphFileName_(madgraphFileName)
  , madgraphIsInitialized_(false)
  , madgraphGluon1P4_(nullptr)
  , madgraphGluon2P4_(nullptr)
  , madgraphChargedLeptonPlusP4_(nullptr)
  , madgraphNeutrinoP4_(nullptr)
  , madgraphChargedLeptonMinusP4_(nullptr)
  , madgraphAntiNeutrinoP4_(nullptr)
  , madgraphBJet1P4_(nullptr)
  , madgraphBJet2P4_(nullptr)
  , errorCode_(0)
  , verbosity_(verbosity)
{
  bjet1TF_ = new BJetTF(verbosity_);
  bjet2TF_ = new BJetTF(verbosity_);
  hadRecoilTF_ = new HadRecoilTF(verbosity_);

  madgraphGluon1P4_ = new double[4];
  madgraphGluon1P4_[1] = 0.;
  madgraphGluon1P4_[2] = 0.;
  madgraphGluon2P4_ = new double[4];
  madgraphGluon2P4_[1] = 0.;
  madgraphGluon2P4_[2] = 0.;  
  madgraphChargedLeptonPlusP4_ = new double[4];
  madgraphNeutrinoP4_ = new double[4];
  madgraphChargedLeptonMinusP4_ = new double[4];
  madgraphAntiNeutrinoP4_ = new double[4];
  madgraphBJet1P4_ = new double[4];
  madgraphBJet2P4_ = new double[4];
}

MEMbbwwIntegrandBase::~MEMbbwwIntegrandBase()
{
  std::cout << "<MEMbbwwIntegrandBase::~MEMbbwwIntegrandBase>:" << std::endl;
  
  delete [] intIntBounds_lower_;
  delete [] intIntBounds_upper_;

  delete bjet1TF_;
  delete bjet2TF_;
  delete hadRecoilTF_;

  delete [] madgraphGluon1P4_;
  delete [] madgraphGluon2P4_;
  delete [] madgraphChargedLeptonPlusP4_;
  delete [] madgraphNeutrinoP4_;
  delete [] madgraphChargedLeptonMinusP4_;
  delete [] madgraphAntiNeutrinoP4_;
  delete [] madgraphBJet1P4_;
  delete [] madgraphBJet2P4_;
}

void
MEMbbwwIntegrandBase::setPDF(LHAPDF::PDF * pdf)
{
  pdf_ = pdf;
}

void 
MEMbbwwIntegrandBase::setInputs(const MeasuredParticle& measuredChargedLeptonPlus, const MeasuredParticle& measuredChargedLeptonMinus,
				const MeasuredParticle& measuredBJet1, const MeasuredParticle& measuredBJet2,
				double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ ) {
    std::cout << "<MEMbbwwIntegrandBase::setInputs>:" << std::endl;
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

void 
MEMbbwwIntegrandBase::printMadGraphMomenta() const
{
  std::cout << "<MEMbbwwIntegrandBase::printMadGraphMomenta>:" << std::endl;
  size_t numParticles = madgraphMomenta_.size();
  for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
    std::cout << "particle #" << idxParticle << ":" 
	      << " En = " << madgraphMomenta_[idxParticle][0] << "," 
	      << " Px = " << madgraphMomenta_[idxParticle][1] << ","
	      << " Py = " << madgraphMomenta_[idxParticle][2] << ","
	      << " Pz = " << madgraphMomenta_[idxParticle][3] << std::endl;
  }
}

