#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandBase.h"

#include <TMath.h>

using namespace mem;

LHAPDF::PDF* MEMbbwwIntegrandBase::pdf_ = nullptr;
std::string MEMbbwwIntegrandBase::pdfName_ = "";
bool MEMbbwwIntegrandBase::pdfIsInitialized_ = false;

MEMbbwwIntegrandBase::MEMbbwwIntegrandBase(double sqrtS, const std::string& pdfName, const std::string& madgraphFileName, int verbosity) 
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

  // initialize PDF set
  if ( !pdfIsInitialized_ ) {
    pdf_ = LHAPDF::mkPDF(pdfName.data(), 0);
    pdfName_ = pdfName;
    pdfIsInitialized_ = true;
  }

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

  delete pdf_;

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
MEMbbwwIntegrandBase::setInputs(const MeasuredParticle& measuredChargedLeptonPlus, const MeasuredParticle& measuredChargedLeptonMinus,
				const MeasuredParticle& measuredBJet1, const MeasuredParticle& measuredBJet2,
				double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ ) {
    std::cout << "<MEMbbwwIntegrandBase::setInputs>:" << std::endl;
  }
std::cout << "break-point 1 reached" << std::endl;
  // reset 'MatrixInversion' error code
  errorCode_ &= (errorCode_ ^ MatrixInversion);
std::cout << "break-point 1.1 reached" << std::endl;
  measuredChargedLeptonPlus_ = measuredChargedLeptonPlus;
std::cout << "break-point 1.2 reached" << std::endl;
  measuredChargedLeptonMinus_ = measuredChargedLeptonMinus;
std::cout << "break-point 1.3 reached" << std::endl;
  measuredBJet1_ = measuredBJet1;
std::cout << "break-point 1.4 reached" << std::endl;
  measuredBJet2_ = measuredBJet2;
std::cout << "break-point 2 reached" << std::endl;
  measuredMEtPx_ = measuredMEtPx;
  measuredMEtPy_ = measuredMEtPy;
  measuredMEtCov_.ResizeTo(2,2);
  measuredMEtCov_ = measuredMEtCov;
std::cout << "break-point 3 reached" << std::endl;
  measuredHadRecoilPx_ = -(measuredChargedLeptonPlus_.px() + measuredChargedLeptonMinus_.px() + measuredBJet1_.px() + measuredBJet2_.px() + measuredMEtPx_);
  measuredHadRecoilPy_ = -(measuredChargedLeptonPlus_.py() + measuredChargedLeptonMinus_.py() + measuredBJet1_.py() + measuredBJet2_.py() + measuredMEtPy_);
std::cout << "break-point 4 reached" << std::endl;
std::cout << "bjet1TF = " << bjet1TF_ << std::endl;
std::cout << "bjet2TF = " << bjet2TF_ << std::endl;
std::cout << "hadRecoilTF = " << hadRecoilTF_ << std::endl;
  // set measured momenta of b-jets and of missing transverse momentum
  // in transfer function (TF) objects
  bjet1TF_->setInputs(measuredBJet1_.p4());
  bjet2TF_->setInputs(measuredBJet2_.p4());
  hadRecoilTF_->setInputs(measuredHadRecoilPx_, measuredHadRecoilPy_, measuredMEtCov_);
  if ( hadRecoilTF_->getErrorCode() ) {
    errorCode_ |= MatrixInversion;
  }
std::cout << "break-point 5 reached" << std::endl;
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

