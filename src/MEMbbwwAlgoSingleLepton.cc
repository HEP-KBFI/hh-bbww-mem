#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoSingleLepton.h"

#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorVEGAS.h"
#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorVAMP.h"

#include "DataFormats/Math/interface/deltaR.h" // deltaR

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>

using namespace mem;

MEMbbwwAlgoSingleLepton::MEMbbwwAlgoSingleLepton(double sqrtS, 
						 const std::string& pdfName, 
						 const std::string& madgraphFileName_signal, const std::string& madgraphFileName_background, 
						 int verbosity) 
  : MEMbbwwAlgoBase(sqrtS, pdfName, madgraphFileName_signal, madgraphFileName_background, verbosity)   
  , integrand_signal_(nullptr)
  , integrand_background_(nullptr)
  , maxNumHadWJetPairs_(8)
  , sortHadJetPairOption_(kSortHadWJetPairsByMass)
{ 
  integrand_signal_ = new MEMbbwwIntegrandSingleLepton_signal(sqrtS_, madgraphFileName_signal_, verbosity_);
  integrand_signal_->setPDF(pdf_);
  integrand_background_ = new MEMbbwwIntegrandSingleLepton_background(sqrtS_, madgraphFileName_background_, verbosity_);
  integrand_background_->setPDF(pdf_);
  
  result_.prob_signal_ = -1.;
  result_.probErr_signal_ = -1.;
  result_.prob_background_ = -1.;
  result_.probErr_background_ = -1.;
}

MEMbbwwAlgoSingleLepton::~MEMbbwwAlgoSingleLepton() 
{
  delete integrand_signal_;
  delete integrand_background_;
}

void 
MEMbbwwAlgoSingleLepton::applyOnshellWmassConstraint_signal(bool flag) 
{ 
  integrand_signal_->applyOnshellWmassConstraint(flag);
}

namespace
{
  struct sortMeasuredParticles 
  {
    bool operator() (const MeasuredParticle& measuredParticle1, const MeasuredParticle& measuredParticle2)
    {
      return ( measuredParticle1.pt() > measuredParticle2.pt() );
    }
  };

  struct sortMeasuredHadWJetPairsByMass
  {
    bool operator() (const MEMbbwwAlgoSingleLepton::MeasuredHadWJetPair& measuredHadWJetPair1, const MEMbbwwAlgoSingleLepton::MeasuredHadWJetPair& measuredHadWJetPair2)
    {
      double deltaMass1 = TMath::Abs(measuredHadWJetPair1.p4().mass() - wBosonMass);
      double deltaMass2 = TMath::Abs(measuredHadWJetPair2.p4().mass() - wBosonMass);
      return ( deltaMass1 < deltaMass2 );
    }
  };
  struct sortMeasuredHadWJetPairsByDeltaR
  {
    bool operator() (const MEMbbwwAlgoSingleLepton::MeasuredHadWJetPair& measuredHadWJetPair1, const MEMbbwwAlgoSingleLepton::MeasuredHadWJetPair& measuredHadWJetPair2)
    {
      double deltaR1 = deltaR(measuredHadWJetPair1.jet1()->p4(), measuredHadWJetPair1.jet2()->p4());
      double deltaR2 = deltaR(measuredHadWJetPair2.jet1()->p4(), measuredHadWJetPair2.jet2()->p4());
      return ( deltaR1 < deltaR2 );
    }
  };
  struct sortMeasuredHadWJetPairsByPt
  {
    bool operator() (const MEMbbwwAlgoSingleLepton::MeasuredHadWJetPair& measuredHadWJetPair1, const MEMbbwwAlgoSingleLepton::MeasuredHadWJetPair& measuredHadWJetPair2)
    {
      double pt1 = measuredHadWJetPair1.p4().pt();
      double pt2 = measuredHadWJetPair2.p4().pt();
      return ( pt1 > pt2 );
    }
  };
  struct sortMeasuredHadWJetPairsByScalarPt
  {
    bool operator() (const MEMbbwwAlgoSingleLepton::MeasuredHadWJetPair& measuredHadWJetPair1, const MEMbbwwAlgoSingleLepton::MeasuredHadWJetPair& measuredHadWJetPair2)
    {
      double sumPt1 = measuredHadWJetPair1.jet1()->pt() + measuredHadWJetPair1.jet2()->pt();
      double sumPt2 = measuredHadWJetPair2.jet1()->pt() + measuredHadWJetPair2.jet2()->pt();
      return ( sumPt1 > sumPt2 );
    }
  };
}

void
MEMbbwwAlgoSingleLepton::setMeasuredParticles(const std::vector<mem::MeasuredParticle>& measuredParticles)
{
  std::vector<MeasuredParticle> measuredParticles_rounded;
  for ( std::vector<MeasuredParticle>::const_iterator measuredParticle = measuredParticles.begin();
	measuredParticle != measuredParticles.end(); ++measuredParticle ) {
    MeasuredParticle measuredParticle_rounded(
      measuredParticle->type(), 
      roundToNdigits(measuredParticle->pt()), 
      roundToNdigits(measuredParticle->eta()), 
      roundToNdigits(measuredParticle->phi()), 
      roundToNdigits(measuredParticle->mass()), 
      measuredParticle->charge());
    measuredParticles_rounded.push_back(measuredParticle_rounded);
  }
  std::sort(measuredParticles_rounded.begin(), measuredParticles_rounded.end(), sortMeasuredParticles());
  measuredParticles_ = measuredParticles_rounded;
  measuredChargedLepton_ = nullptr;
  measuredHadWJets_.clear(); 
  measuredLeadingBJet_ = nullptr; 
  measuredSubleadingBJet_ = nullptr;
  for ( size_t idx = 0; idx < measuredParticles_.size(); ++idx ) {
    const MeasuredParticle& measuredParticle = measuredParticles_[idx];
    if ( verbosity_ >= 1 ) {
      std::cout << "measuredParticles #" << idx << " (type = " << measuredParticle.type() << "): Pt = " << measuredParticle.pt() << "," 
		<< " eta = " << measuredParticle.eta() << " (theta = " << measuredParticle.p3().theta() << ")" << ", phi = " << measuredParticle.phi() << "," 
		<< " mass = " << measuredParticle.mass() << ", charge = " << measuredParticle.charge() << std::endl;
    }
    if ( measuredParticle.type() == MeasuredParticle::kElectron || measuredParticle.type() == MeasuredParticle::kMuon ) {
      measuredChargedLepton_ = &measuredParticle;
    }
    if ( measuredParticle.type() == MeasuredParticle::kBJet ) {
      if      ( !measuredLeadingBJet_    ) measuredLeadingBJet_    = &measuredParticle;
      else if ( !measuredSubleadingBJet_ ) measuredSubleadingBJet_ = &measuredParticle;
    }
    if ( measuredParticle.type() == MeasuredParticle::kHadWJet ) {
      measuredHadWJets_.push_back(&measuredParticle);
    }
  }
  if ( !(measuredChargedLepton_ && measuredLeadingBJet_ && measuredSubleadingBJet_) ) {
    std::cerr << "<MEMbbwwAlgoSingleLepton::integrate>: Given measuredParticles are not of the expected type --> ABORTING !!\n";
    assert(0);
  }
  for ( std::vector<const mem::MeasuredParticle*>::const_iterator measuredHadWJet1 = measuredHadWJets_.begin();
	measuredHadWJet1 != measuredHadWJets_.end(); ++measuredHadWJet1 ) {
    for ( std::vector<const mem::MeasuredParticle*>::const_iterator measuredHadWJet2 = measuredHadWJets_.begin();
	measuredHadWJet2 != measuredHadWJets_.end(); ++measuredHadWJet2 ) {
      // CV: add both permutations of measuredHadWJet1 and measuredHadWJet2
      measuredHadWJetPairs_.push_back(MeasuredHadWJetPair(**measuredHadWJet1, **measuredHadWJet2));
    }
  }
  if ( sortHadJetPairOption_ == kSortHadWJetPairsByMass ) {
    std::sort(measuredHadWJetPairs_.begin(), measuredHadWJetPairs_.end(), sortMeasuredHadWJetPairsByMass());
  } else if ( sortHadJetPairOption_ == kSortHadWJetPairsByDeltaR ) {
    std::sort(measuredHadWJetPairs_.begin(), measuredHadWJetPairs_.end(), sortMeasuredHadWJetPairsByDeltaR());
  } else if ( sortHadJetPairOption_ == kSortHadWJetPairsByPt ) {
    std::sort(measuredHadWJetPairs_.begin(), measuredHadWJetPairs_.end(), sortMeasuredHadWJetPairsByPt());
  } else if ( sortHadJetPairOption_ == kSortHadWJetPairsByScalarPt ) {
    std::sort(measuredHadWJetPairs_.begin(), measuredHadWJetPairs_.end(), sortMeasuredHadWJetPairsByScalarPt());
  } else {
    std::cerr << "<MEMbbwwAlgoSingleLepton::integrate>: Invalid configuration parameter 'sortHadJetPairOption' = " << sortHadJetPairOption_ << " --> ABORTING !!\n";
    assert(0);
  }
  if ( !(measuredHadWJetPairs_.size() >= 1) ) {
    std::cerr << "<MEMbbwwAlgoSingleLepton::integrate>: Given measuredParticles do not contain at least two of type 'HadWJet' --> ABORTING !!\n";
    assert(0);
  }
}

void
MEMbbwwAlgoSingleLepton::integrate(const std::vector<MeasuredParticle>& measuredParticles, double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<MEMbbwwAlgoSingleLepton::integrate>:" << std::endl;
  }

  clock_->Reset();
  //const std::string label_total = "<MEMbbwwAlgoSingleLepton::integrate (total)>";
  const std::string label_total = "<MEMbbwwAlgoSingleLepton::integrate>";
  clock_->Start(label_total.data());

  setMeasuredParticles(measuredParticles);
  setMeasuredMEt_and_Cov(measuredMEtPx, measuredMEtPy, measuredMEtCov);

  size_t numHadWJetPairs = measuredHadWJetPairs_.size();
  if ( maxNumHadWJetPairs_ > 0 ) numHadWJetPairs = std::min(numHadWJetPairs, (size_t)maxNumHadWJetPairs_);
  
  //const std::string label_signal = "<MEMbbwwAlgoSingleLepton::integrate (signal hypothesis)>";
  //clock_->Start(label_signal.data());
  result_.prob_signal_ = 0.;
  result_.probErr_signal_ = 0.;
  result_.permutations_signal_.clear();
  integrand_signal_->resetNumMatrixElementEvaluations();
  unsigned numPermutations_signal = 0;
  for ( size_t idxHadWJetPair = 0; idxHadWJetPair < numHadWJetPairs; ++idxHadWJetPair ) {
    const MeasuredParticle* measuredHadWJet1 = measuredHadWJetPairs_[idxHadWJetPair].jet1();
    const MeasuredParticle* measuredHadWJet2 = measuredHadWJetPairs_[idxHadWJetPair].jet2();
    for ( unsigned idxPermutation = 0; idxPermutation < 4; ++idxPermutation ) {
      const MeasuredParticle* measuredBJet1 = nullptr;
      const MeasuredParticle* measuredBJet2 = nullptr;
      if ( idxPermutation == 0 || idxPermutation == 2 ) {
	measuredBJet1 = measuredLeadingBJet_;
	measuredBJet2 = measuredSubleadingBJet_;
      } else {
	measuredBJet1 = measuredSubleadingBJet_;
	measuredBJet2 = measuredLeadingBJet_;
      }
      integrand_signal_->setInputs(
        *measuredChargedLepton_, *measuredBJet1, *measuredBJet2, 
	*measuredHadWJet1, *measuredHadWJet2,
        measuredMEtPx_, measuredMEtPy_, measuredMEtCov_);
      int chargedLeptonPermutation = kPermutationUndefined1L;
      if ( idxPermutation == 0 || idxPermutation == 1 ) {
	chargedLeptonPermutation = kOnshellChargedLepton;
      } else {
	chargedLeptonPermutation = kOffshellChargedLepton;
      }
      integrand_signal_->setOnshellChargedLepton(chargedLeptonPermutation);
      MEMbbwwAlgoSingleLepton::gMEMIntegrand = integrand_signal_;
      initializeIntAlgo(maxObjFunctionCalls_signal_);
      double prob_permutation, probErr_permutation;
      MEMbbwwAlgoSingleLepton::runIntAlgo(integrand_signal_, prob_permutation, probErr_permutation);
      result_.prob_signal_ += prob_permutation;
      result_.probErr_signal_ += probErr_permutation;
      result_.permutations_signal_.push_back(MEMbbwwPermutationSingleLepton(
        prob_permutation, probErr_permutation, 
        *measuredChargedLepton_, *measuredBJet1, *measuredBJet2, *measuredHadWJet1, *measuredHadWJet2,
	chargedLeptonPermutation));
      ++numPermutations_signal;
      delete intAlgo_;
      intAlgo_ = nullptr;
    }
  }
  result_.prob_signal_ /= numPermutations_signal;
  result_.probErr_signal_ /= numPermutations_signal;
  numMatrixElementEvaluations_signal_ = integrand_signal_->getNumMatrixElementEvaluations();
  //clock_->Stop(label_signal.data().data());
  //if ( verbosity_ >= 0 ) {
  //  clock_->Show(label_signal.data().data());
  //}
 
  //const std::string label_background = "<MEMbbwwAlgoSingleLepton::integrate (background hypothesis)>";
  //clock_->Start(label_background.data());
  result_.prob_background_ = 0.;
  result_.probErr_background_ = 0.;
  result_.permutations_background_.clear();
  integrand_background_->resetNumMatrixElementEvaluations();
  unsigned numPermutations_background = 0;
  for ( size_t idxHadWJetPair = 0; idxHadWJetPair < numHadWJetPairs; ++idxHadWJetPair ) {
    const MeasuredParticle* measuredHadWJet1 = measuredHadWJetPairs_[idxHadWJetPair].jet1();
    const MeasuredParticle* measuredHadWJet2 = measuredHadWJetPairs_[idxHadWJetPair].jet2();
    for ( unsigned idxPermutation = 0; idxPermutation < 2; ++idxPermutation ) {
      const MeasuredParticle* measuredBJetFromTop = nullptr;
      const MeasuredParticle* measuredBJetFromAntiTop = nullptr;
      if ( idxPermutation == 0 ) {
	measuredBJetFromTop = measuredLeadingBJet_;
	measuredBJetFromAntiTop  = measuredSubleadingBJet_;
      } else {
	measuredBJetFromTop = measuredSubleadingBJet_;
	measuredBJetFromAntiTop  = measuredLeadingBJet_;
      }
      integrand_background_->setInputs(
        *measuredChargedLepton_, *measuredBJetFromTop, *measuredBJetFromAntiTop, 
	*measuredHadWJet1, *measuredHadWJet2,
        measuredMEtPx_, measuredMEtPy_, measuredMEtCov_);
      MEMbbwwAlgoSingleLepton::gMEMIntegrand = integrand_background_;
      initializeIntAlgo(maxObjFunctionCalls_background_);
      double prob_permutation, probErr_permutation;
      MEMbbwwAlgoSingleLepton::runIntAlgo(integrand_background_, prob_permutation, probErr_permutation);
      result_.prob_background_ += prob_permutation;
      result_.probErr_background_ += probErr_permutation;
      result_.permutations_background_.push_back(MEMbbwwPermutationSingleLepton(
	prob_permutation, probErr_permutation, 
        *measuredChargedLepton_, *measuredBJetFromTop, *measuredBJetFromAntiTop, *measuredHadWJet1, *measuredHadWJet2));
      ++numPermutations_background;								
      delete intAlgo_;
      intAlgo_ = nullptr;
    }
  }
  result_.prob_background_ /= numPermutations_background;
  result_.probErr_background_ /= numPermutations_background;
  numMatrixElementEvaluations_background_ = integrand_background_->getNumMatrixElementEvaluations();
  //clock_->Stop(label_background.data());
  //if ( verbosity_ >= 0 ) {
  //  clock_->Show(label_background.data());
  //}

  clock_->Stop(label_total.data());
  if ( verbosity_ >= 0 ) {
    clock_->Show(label_total.data());
  }
  numSeconds_cpu_ = clock_->GetCpuTime(label_total.data());
  numSeconds_real_ = clock_->GetRealTime(label_total.data());
}
