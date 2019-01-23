#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandBase.h"

#include <TMath.h>

using namespace mem;

MEMbbwwIntegrandBase::MEMbbwwIntegrandBase(double sqrtS, const std::string& madgraphFileName, int verbosity)
  : intBounds_lower_(nullptr)
  , intBounds_upper_(nullptr)
  , measuredBJet1_(nullptr)
  , measuredBJet2_(nullptr)
  , bjet1TF_(nullptr)
  , bjet2TF_(nullptr)
  , hadRecoilTF_(nullptr)
  , sqrtS_(sqrtS)
  , normFactor_(-1.)
  , madgraphFileName_(madgraphFileName)
  , madgraphIsInitialized_(false)
  , madgraphGluon1P4_(nullptr)
  , madgraphGluon2P4_(nullptr)
  , madgraphBJet1P4_(nullptr)
  , madgraphBJet2P4_(nullptr)
  , numMatrixElementEvaluations_(0)
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
  madgraphBJet1P4_ = new double[4];
  madgraphBJet2P4_ = new double[4];
}

MEMbbwwIntegrandBase::~MEMbbwwIntegrandBase()
{
  //if ( verbosity_ >= 1 ) 
  //{
  //  std::cout << "<MEMbbwwIntegrandBase::~MEMbbwwIntegrandBase>:" << std::endl;
  //}
  
  delete [] intBounds_lower_;
  delete [] intBounds_upper_;

  delete bjet1TF_;
  delete bjet2TF_;
  delete hadRecoilTF_;

  delete [] madgraphGluon1P4_;
  delete [] madgraphGluon2P4_;
  delete [] madgraphBJet1P4_;
  delete [] madgraphBJet2P4_;
}

void
MEMbbwwIntegrandBase::setPDF(LHAPDF::PDF * pdf)
{
  pdf_ = pdf;
}

void 
MEMbbwwIntegrandBase::printMadGraphMomenta(const std::vector<double*>& madgraphMomenta) const
{
  std::cout << "<MEMbbwwIntegrandBase::printMadGraphMomenta>:" << std::endl;
  size_t numParticles = madgraphMomenta.size();
  for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
    std::cout << "particle #" << idxParticle << ":" 
	      << " En = " << madgraphMomenta[idxParticle][0] << "," 
	      << " Px = " << madgraphMomenta[idxParticle][1] << ","
	      << " Py = " << madgraphMomenta[idxParticle][2] << ","
	      << " Pz = " << madgraphMomenta[idxParticle][3] << std::endl;
  }
}

