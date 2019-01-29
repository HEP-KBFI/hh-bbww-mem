
#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoDilepton.h"

#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorVEGAS.h"
#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorVAMP.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>

using namespace mem;

MEMbbwwAlgoDilepton::MEMbbwwAlgoDilepton(double sqrtS, 
					 const std::string& pdfName, 
					 const std::string& madgraphFileName_signal, const std::string& madgraphFileName_background, 
					 int verbosity) 
  : MEMbbwwAlgoBase(sqrtS, pdfName, madgraphFileName_signal, madgraphFileName_background, verbosity) 
  , integrand_signal_(nullptr)
  , integrand_background_(nullptr)
{ 
  integrand_signal_ = new MEMbbwwIntegrandDilepton_signal(sqrtS_, madgraphFileName_signal_, verbosity_);
  integrand_signal_->setPDF(pdf_);
  integrand_background_ = new MEMbbwwIntegrandDilepton_background(sqrtS_, madgraphFileName_background_, verbosity_);
  integrand_background_->setPDF(pdf_);
  
  result_.prob_signal_ = -1.;
  result_.probErr_signal_ = -1.;
  result_.prob_background_ = -1.;
  result_.probErr_background_ = -1.;
}

MEMbbwwAlgoDilepton::~MEMbbwwAlgoDilepton() 
{
  delete integrand_signal_;
  delete integrand_background_;
}

void 
MEMbbwwAlgoDilepton::applyOnshellWmassConstraint_signal(bool flag) 
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
}

void
MEMbbwwAlgoDilepton::setMeasuredParticles(const std::vector<MeasuredParticle>& measuredParticles)
{
  std::vector<MeasuredParticle> measuredParticles_rounded;
  for ( std::vector<MeasuredParticle>::const_iterator measuredParticle = measuredParticles.begin();
	measuredParticle != measuredParticles.end(); ++measuredParticle ) 
  {
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
  measuredChargedLeptonPlus_ = nullptr;
  measuredChargedLeptonMinus_ = nullptr;
  measuredLeadingBJet_ = nullptr; 
  measuredSubleadingBJet_ = nullptr;
  for ( size_t idx = 0; idx < measuredParticles_.size(); ++idx ) 
  {
    const MeasuredParticle& measuredParticle = measuredParticles_[idx];
    if ( verbosity_ >= 1 ) 
    {
      std::cout << "measuredParticles #" << idx << " (type = " << measuredParticle.type() << "): Pt = " << measuredParticle.pt() << "," 
		<< " eta = " << measuredParticle.eta() << " (theta = " << measuredParticle.p3().theta() << ")" << ", phi = " << measuredParticle.phi() << "," 
		<< " mass = " << measuredParticle.mass() << ", charge = " << measuredParticle.charge() << std::endl;
    }
    if ( measuredParticle.type() == MeasuredParticle::kElectron || measuredParticle.type() == MeasuredParticle::kMuon ) 
    {
      if ( measuredParticle.charge() > 0 ) measuredChargedLeptonPlus_  = &measuredParticle;
      if ( measuredParticle.charge() < 0 ) measuredChargedLeptonMinus_ = &measuredParticle;
    }
    if ( measuredParticle.type() == MeasuredParticle::kBJet ) 
    {
      if      ( !measuredLeadingBJet_    ) measuredLeadingBJet_    = &measuredParticle;
      else if ( !measuredSubleadingBJet_ ) measuredSubleadingBJet_ = &measuredParticle;
    }
  }
  if ( !(measuredChargedLeptonPlus_ && measuredChargedLeptonMinus_ && measuredLeadingBJet_) ) // CV: allow for one "missing" (non-reconstructed) b-jet 
  {
    std::cerr << "<MEMbbwwAlgoDilepton::integrate>: Given measuredParticles are not of the expected type --> ABORTING !!\n";
    assert(0);
  }
 
}

void
MEMbbwwAlgoDilepton::integrate(const std::vector<MeasuredParticle>& measuredParticles, double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ >= 1 ) 
  {
    std::cout << "<MEMbbwwAlgoDilepton::integrate>:" << std::endl;
  }

  clock_->Reset();
  //const std::string label_total = "<MEMbbwwAlgoDilepton::integrate (total)>";
  const std::string label_total = "<MEMbbwwAlgoDilepton::integrate>";
  clock_->Start(label_total.data());

  setMeasuredParticles(measuredParticles);
  setMeasuredMEt_and_Cov(measuredMEtPx, measuredMEtPy, measuredMEtCov);

  //const std::string label_signal = "<MEMbbwwAlgoDilepton::integrate (signal hypothesis)>";
  //clock_->Start(label_signal.data());
  result_.prob_signal_ = 0.;
  result_.probErr_signal_ = 0.;
  result_.permutations_signal_.clear();
  integrand_signal_->resetNumMatrixElementEvaluations();
  for ( unsigned idxPermutation = 0; idxPermutation < 2; ++idxPermutation ) 
  {
    //---------------------------------------------------------------------------
    // CV: MadGraph matrix element for HH->bbWW signal is symmetric under exchange b-jet1 vs b-jet2,
    //     so one association of 
    //       measuredBJet1 & measuredBJet1 to measuredLeadingBJet & measuredSubleadingBJet
    //     is sufficient
    //    (reduce number of permutations as much as possible, to save computing time !!)
    const MeasuredParticle* measuredBJet1 = measuredLeadingBJet_;
    const MeasuredParticle* measuredBJet2 = measuredSubleadingBJet_;
    //---------------------------------------------------------------------------
    integrand_signal_->setInputs(
      measuredChargedLeptonPlus_, measuredChargedLeptonMinus_, 
      measuredBJet1, measuredBJet2, 
      measuredMEtPx_, measuredMEtPy_, measuredMEtCov_);
    int chargedLeptonPermutation = mem::kPermutationUndefined2L;
    if ( idxPermutation == 0 ) 
    {
      chargedLeptonPermutation = mem::kOnshellChargedLeptonPlus;
    } 
    else 
    {
      chargedLeptonPermutation = mem::kOnshellChargedLeptonMinus;
    }
    integrand_signal_->setOnshellChargedLepton(chargedLeptonPermutation);
    MEMbbwwAlgoDilepton::gMEMIntegrand = integrand_signal_;
    initializeIntAlgo(maxObjFunctionCalls_signal_);
    double prob_permutation, probErr_permutation;
    MEMbbwwAlgoDilepton::runIntAlgo(integrand_signal_, prob_permutation, probErr_permutation);
    result_.prob_signal_ += prob_permutation;
    result_.probErr_signal_ += probErr_permutation;
    const MeasuredParticle* measuredChargedLeptonFromOnshellW = nullptr;
    const MeasuredParticle* measuredChargedLeptonFromOffshellW = nullptr;
    if ( chargedLeptonPermutation == mem::kOnshellChargedLeptonPlus ) 
    {
      measuredChargedLeptonFromOnshellW = measuredChargedLeptonPlus_;
      measuredChargedLeptonFromOffshellW = measuredChargedLeptonMinus_;
    } 
    else 
    {
      measuredChargedLeptonFromOnshellW = measuredChargedLeptonMinus_;
      measuredChargedLeptonFromOffshellW = measuredChargedLeptonPlus_;
    }
    result_.permutations_signal_.push_back(MEMbbwwPermutationDilepton(
      prob_permutation, probErr_permutation, 
      measuredBJet1, measuredBJet2, measuredChargedLeptonFromOnshellW, measuredChargedLeptonFromOffshellW,
      chargedLeptonPermutation));
    delete intAlgo_;
    intAlgo_ = nullptr;
  }
  result_.prob_signal_ /= 2.;
  result_.probErr_signal_ /= 2.;
  numMatrixElementEvaluations_signal_ = integrand_signal_->getNumMatrixElementEvaluations();
  //clock_->Stop(label_signal.data().data());
  //if ( verbosity_ >= 0 ) 
  //{
  //  clock_->Show(label_signal.data().data());
  //}

  //const std::string label_background = "<MEMbbwwAlgoDilepton::integrate (background hypothesis)>";
  //clock_->Start(label_background.data());
  result_.prob_background_ = 0.;
  result_.probErr_background_ = 0.;
  result_.permutations_background_.clear();
  integrand_background_->resetNumMatrixElementEvaluations();
  for ( unsigned idxPermutation = 0; idxPermutation < 2; ++idxPermutation ) 
  {
    const MeasuredParticle* measuredBJetFromTop = nullptr;
    const MeasuredParticle* measuredBJetFromAntiTop = nullptr;
    if ( idxPermutation == 0 ) 
    {
      measuredBJetFromTop = measuredLeadingBJet_;
      measuredBJetFromAntiTop  = measuredSubleadingBJet_;
    } 
    else 
    {
      measuredBJetFromTop = measuredSubleadingBJet_;
      measuredBJetFromAntiTop  = measuredLeadingBJet_;
    }
    integrand_background_->setInputs(
      measuredChargedLeptonPlus_, measuredChargedLeptonMinus_, 
      measuredBJetFromTop, measuredBJetFromAntiTop, 
      measuredMEtPx_, measuredMEtPy_, measuredMEtCov_);
    MEMbbwwAlgoDilepton::gMEMIntegrand = integrand_background_;
    initializeIntAlgo(maxObjFunctionCalls_background_);
    double prob_permutation, probErr_permutation;
    MEMbbwwAlgoDilepton::runIntAlgo(integrand_background_, prob_permutation, probErr_permutation);
    result_.prob_background_ += prob_permutation;
    result_.probErr_background_ += probErr_permutation;
    result_.permutations_background_.push_back(MEMbbwwPermutationDilepton(
      prob_permutation, probErr_permutation, 
      measuredBJetFromTop, measuredBJetFromAntiTop, measuredChargedLeptonPlus_, measuredChargedLeptonMinus_));
    delete intAlgo_;
    intAlgo_ = nullptr;
  }
  result_.prob_background_ /= 2.;
  result_.probErr_background_ /= 2.;
  numMatrixElementEvaluations_background_ = integrand_background_->getNumMatrixElementEvaluations();
  //clock_->Stop(label_background.data());
  //if ( verbosity_ >= 0 ) 
  //{
  //  clock_->Show(label_background.data());
  //}

  clock_->Stop(label_total.data());
  if ( verbosity_ >= 0 ) 
  {
    clock_->Show(label_total.data());
  }
  numSeconds_cpu_ = clock_->GetCpuTime(label_total.data());
  numSeconds_real_ = clock_->GetRealTime(label_total.data());
}
