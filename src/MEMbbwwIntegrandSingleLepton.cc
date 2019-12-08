#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandSingleLepton.h"

#include <TMath.h>

using namespace mem;

MEMbbwwIntegrandSingleLepton::MEMbbwwIntegrandSingleLepton(double sqrtS, const std::string& madgraphFileName, int verbosity)
  : MEMbbwwIntegrandBase(sqrtS, madgraphFileName, verbosity)
  , measuredChargedLepton_(nullptr)
  , measuredHadWJet1_(nullptr)
  , measuredHadWJet2_(nullptr)
  , hadWJet1TF_(nullptr)
  , hadWJet2TF_(nullptr)
  , madgraphChargedLeptonP4_(nullptr)
  , madgraphNeutrinoP4_(nullptr)
  , madgraphHadWJet1P4_(nullptr)
  , madgraphHadWJet2P4_(nullptr)
{
  hadWJet1TF_ = new HadWJetTF(verbosity_);
  hadWJet2TF_ = new HadWJetTF(verbosity_);

  madgraphChargedLeptonP4_ = new double[4];
  madgraphNeutrinoP4_ = new double[4];
  madgraphHadWJet1P4_ = new double[4];
  madgraphHadWJet2P4_ = new double[4];
}

MEMbbwwIntegrandSingleLepton::~MEMbbwwIntegrandSingleLepton()
{
  //if ( verbosity_ >= 1 ) 
  //{
  //  std::cout << "<MEMbbwwIntegrandSingleLepton::~MEMbbwwIntegrandSingleLepton>:" << std::endl;
  //}

  delete hadWJet1TF_;
  delete hadWJet2TF_;
  
  delete [] madgraphChargedLeptonP4_;
  delete [] madgraphNeutrinoP4_;
  delete [] madgraphHadWJet1P4_;
  delete [] madgraphHadWJet2P4_;
}

void 
MEMbbwwIntegrandSingleLepton::setInputs(const MeasuredParticle* measuredChargedLepton, 
					const MeasuredParticle* measuredHadWJet1, const MeasuredParticle* measuredHadWJet2,
					const MeasuredParticle* measuredBJet1, const MeasuredParticle* measuredBJet2,
					double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  //if ( verbosity_ >= 1 ) 
  //{
  //  std::cout << "<MEMbbwwIntegrandSingleLepton::setInputs>:" << std::endl;
  //}
  // reset 'MatrixInversion' error code
  errorCode_ &= (errorCode_ ^ MatrixInversion);
  measuredChargedLepton_ = measuredChargedLepton;
  measuredHadWJet1_ = measuredHadWJet1;
  measuredHadWJet2_ = measuredHadWJet2;
  measuredBJet1_ = measuredBJet1;
  measuredBJet2_ = measuredBJet2;
  measuredMEtPx_ = measuredMEtPx;
  measuredMEtPy_ = measuredMEtPy;
  measuredMEtCov_.ResizeTo(2,2);
  measuredMEtCov_ = measuredMEtCov;
  measuredHadRecoilPx_ = -measuredMEtPx_;
  measuredHadRecoilPy_ = -measuredMEtPy_;
  if ( measuredChargedLepton_ ) 
  {
    measuredHadRecoilPx_ -= measuredChargedLepton_->px();
    measuredHadRecoilPy_ -= measuredChargedLepton_->py();    
  }
  if ( measuredHadWJet1_ ) 
  {
    measuredHadRecoilPx_ -= measuredHadWJet1_->px();
    measuredHadRecoilPy_ -= measuredHadWJet1_->py();    
  }
  if ( measuredHadWJet2_ ) 
  {
    measuredHadRecoilPx_ -= measuredHadWJet2_->px();
    measuredHadRecoilPy_ -= measuredHadWJet2_->py();
  }
  if ( measuredBJet1_ ) 
  {
    measuredHadRecoilPx_ -= measuredBJet1_->px();
    measuredHadRecoilPy_ -= measuredBJet1_->py();    
  }
  if ( measuredBJet2_ ) 
  {
    measuredHadRecoilPx_ -= measuredBJet2_->px();
    measuredHadRecoilPy_ -= measuredBJet2_->py();
  }
  // set measured momenta of jets from W->jj decay, of b-jets, and of missing transverse momentum
  // in transfer function (TF) objects
  if ( measuredHadWJet1_ ) 
  {
    hadWJet1TF_->setInputs(measuredHadWJet1_->p4());
  }
  if ( measuredHadWJet2_ )
  {
    hadWJet2TF_->setInputs(measuredHadWJet2_->p4());
  }
  if ( measuredBJet1_ ) 
  {
    bjet1TF_->setInputs(measuredBJet1_->p4());
  }
  if ( measuredBJet2_ ) 
  {
    bjet2TF_->setInputs(measuredBJet2_->p4());
  }
  hadRecoilTF_->setInputs(measuredHadRecoilPx_, measuredHadRecoilPy_, measuredMEtCov_);
  if ( hadRecoilTF_->getErrorCode() ) 
  {
    errorCode_ |= MatrixInversion;
  }
}
