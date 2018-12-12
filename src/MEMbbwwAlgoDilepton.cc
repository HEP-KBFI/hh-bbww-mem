#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoDilepton.h"

#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorVEGAS.h"
#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorVAMP.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>

using namespace mem;

const MEMbbwwIntegrandBase* MEMbbwwAlgoDilepton::gMEMIntegrand = nullptr;
LHAPDF::PDF * MEMbbwwAlgoDilepton::pdf_ = nullptr;

namespace 
{
  double g_C(double* x, size_t, void*)
  {    
    //std::cout << "<g_C>:" << std::endl;
    double retVal = MEMbbwwAlgoDilepton::gMEMIntegrand->Eval(x);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }

  double g_Fortran(double** x, size_t, void**)
  {    
    //std::cout << "<g_Fortran>:" << std::endl;
    double retVal = MEMbbwwAlgoDilepton::gMEMIntegrand->Eval(*x);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }
}

MEMbbwwAlgoDilepton::MEMbbwwAlgoDilepton(double sqrtS, 
					 const std::string& pdfName, 
					 const std::string& madgraphFileName_signal, const std::string& madgraphFileName_background, 
					 int verbosity) 
  : madgraphFileName_signal_(madgraphFileName_signal)
  , integrand_signal_(nullptr)
  , madgraphFileName_background_(madgraphFileName_background)
  , integrand_background_(nullptr)
  , sqrtS_(sqrtS)
  , intMode_(kVAMP)
  , intAlgo_(nullptr)
  , maxObjFunctionCalls_(20000)
  , precision_(1.e-3)
  , clock_(nullptr)
  , numSeconds_cpu_(-1.)
  , numSeconds_real_(-1.)
  , numMatrixElementEvaluations_signal_(0)
  , numMatrixElementEvaluations_background_(0)
  , verbosity_(verbosity)
{ 
  integrand_signal_ = new MEMbbwwIntegrandDilepton_signal(sqrtS_, madgraphFileName_signal_, verbosity_);
  integrand_background_ = new MEMbbwwIntegrandDilepton_background(sqrtS_, madgraphFileName_background_, verbosity_);

  if ( !pdf_ )
  {
    pdf_ = LHAPDF::mkPDF(pdfName);
  }
  integrand_signal_->setPDF(pdf_);
  integrand_background_->setPDF(pdf_);
  
  result_.prob_signal_ = -1.;
  result_.probErr_signal_ = -1.;
  result_.prob_background_ = -1.;
  result_.probErr_background_ = -1.;

  clock_ = new TBenchmark();
}

MEMbbwwAlgoDilepton::~MEMbbwwAlgoDilepton() 
{
  delete integrand_signal_;
  delete integrand_background_;

  if ( pdf_ )
  {
    delete pdf_;
    pdf_ = nullptr;
  }

  delete clock_;
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
MEMbbwwAlgoDilepton::integrate(const std::vector<MeasuredParticle>& measuredParticles, double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<MEMbbwwAlgoDilepton::integrate>:" << std::endl;
  }

  clock_->Reset();
  clock_->Start("<MEMbbwwAlgoDilepton::integrate>");
  
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
  measuredChargedLeptonPlus_ = nullptr;
  measuredChargedLeptonMinus_ = nullptr;
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
      if ( measuredParticle.charge() > 0 ) measuredChargedLeptonPlus_  = &measuredParticle;
      if ( measuredParticle.charge() < 0 ) measuredChargedLeptonMinus_ = &measuredParticle;
    }
    if ( measuredParticle.type() == MeasuredParticle::kBJet ) {
      if      ( !measuredLeadingBJet_    ) measuredLeadingBJet_    = &measuredParticle;
      else if ( !measuredSubleadingBJet_ ) measuredSubleadingBJet_ = &measuredParticle;
    }
  }
  if ( !(measuredChargedLeptonPlus_ && measuredChargedLeptonMinus_ && measuredLeadingBJet_ && measuredSubleadingBJet_) ) {
    std::cerr << "<MEMbbwwAlgoDilepton::integrate>: Given measuredParticles are not of the expected type --> ABORTING !!\n";
    assert(0);
  }
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

  result_.prob_signal_ = 0.;
  result_.probErr_signal_ = 0.;
  result_.permutations_signal_.clear();
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
      *measuredChargedLeptonPlus_, *measuredChargedLeptonMinus_, *measuredBJet1, *measuredBJet2, 
      measuredMEtPx_, measuredMEtPy_, measuredMEtCov_);
    int chargedLeptonPermutation = MEMbbwwIntegrandDilepton_signal::kPermutationUndefined;
    if ( idxPermutation == 0 || idxPermutation == 1 ) {
      chargedLeptonPermutation = MEMbbwwIntegrandDilepton_signal::kOnshellChargedLeptonPlus;
    } else {
      chargedLeptonPermutation = MEMbbwwIntegrandDilepton_signal::kOnshellChargedLeptonMinus;
    }
    integrand_signal_->setOnshellChargedLepton(chargedLeptonPermutation);
    MEMbbwwAlgoDilepton::gMEMIntegrand = integrand_signal_;
    initializeIntAlgo();
    double prob_permutation, probErr_permutation;
    MEMbbwwAlgoDilepton::runIntAlgo(integrand_signal_, prob_permutation, probErr_permutation);
    result_.prob_signal_ += prob_permutation;
    result_.probErr_signal_ += probErr_permutation;
    const MeasuredParticle* measuredChargedLeptonFromOnshellW = nullptr;
    const MeasuredParticle* measuredChargedLeptonFromOffshellW = nullptr;
    if ( chargedLeptonPermutation == MEMbbwwIntegrandDilepton_signal::kOnshellChargedLeptonPlus ) {
      measuredChargedLeptonFromOnshellW = measuredChargedLeptonPlus_;
      measuredChargedLeptonFromOffshellW = measuredChargedLeptonMinus_;
    } else {
      measuredChargedLeptonFromOnshellW = measuredChargedLeptonMinus_;
      measuredChargedLeptonFromOffshellW = measuredChargedLeptonPlus_;
    }
    result_.permutations_signal_.push_back(MEMbbwwPermutationDilepton_sig(
      prob_permutation, probErr_permutation, 
      *measuredBJet1, *measuredBJet2, *measuredChargedLeptonFromOnshellW, *measuredChargedLeptonFromOffshellW));
    delete intAlgo_;
    intAlgo_ = nullptr;
  }
  result_.prob_signal_ /= 4.;
  result_.probErr_signal_ /= 4.;
  numMatrixElementEvaluations_signal_ = integrand_signal_->getNumMatrixElementEvaluations();
 
  result_.prob_background_ = 0.;
  result_.probErr_background_ = 0.;
  result_.permutations_background_.clear();
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
      *measuredChargedLeptonPlus_, *measuredChargedLeptonMinus_, *measuredBJetFromTop, *measuredBJetFromAntiTop, 
      measuredMEtPx_, measuredMEtPy_, measuredMEtCov_);
    MEMbbwwAlgoDilepton::gMEMIntegrand = integrand_background_;
    initializeIntAlgo();
    double prob_permutation, probErr_permutation;
    MEMbbwwAlgoDilepton::runIntAlgo(integrand_background_, prob_permutation, probErr_permutation);
    result_.prob_background_ += prob_permutation;
    result_.probErr_background_ += probErr_permutation;
    result_.permutations_background_.push_back(MEMbbwwPermutationDilepton_bkg(
      prob_permutation, probErr_permutation, 
      *measuredBJetFromTop, *measuredChargedLeptonPlus_, *measuredBJetFromAntiTop, *measuredChargedLeptonMinus_));
    delete intAlgo_;
    intAlgo_ = nullptr;
  }
  result_.prob_background_ /= 2.;
  result_.probErr_background_ /= 2.;
  numMatrixElementEvaluations_background_ = integrand_background_->getNumMatrixElementEvaluations();

  clock_->Stop("<MEMbbwwAlgoDilepton::integrate>");
  if ( verbosity_ >= 0 ) {
    clock_->Show("<MEMbbwwAlgoDilepton::integrate>");
  }
  numSeconds_cpu_ = clock_->GetCpuTime("<MEMbbwwAlgoDilepton::integrate>");
  numSeconds_real_ = clock_->GetRealTime("<MEMbbwwAlgoDilepton::integrate>");
}

void MEMbbwwAlgoDilepton::initializeIntAlgo()
{
  if ( intMode_ == kVEGAS ) {
    unsigned numCallsGridOpt = TMath::Nint(0.20*maxObjFunctionCalls_);
    unsigned numCallsIntEval = TMath::Nint(0.80*maxObjFunctionCalls_);
    intAlgo_ = new MEMIntegratorVEGAS(
      numCallsGridOpt, numCallsIntEval, 
      2., 1);
  } else if ( intMode_ == kVAMP ) {
    unsigned numCallsGridOpt = TMath::Nint(0.20*maxObjFunctionCalls_);
    unsigned numCallsIntEval = TMath::Nint(0.80*maxObjFunctionCalls_);
    intAlgo_ = new MEMIntegratorVAMP(
      numCallsGridOpt, numCallsIntEval);
  } else {
    std::cerr << "<MEMbbwwAlgoDilepton::initializeIntAlgo>: Invalid configuration parameter 'intMode' = " << intMode_ << " --> ABORTING !!\n";
    assert(0);
  }
}

void MEMbbwwAlgoDilepton::runIntAlgo(mem::MEMbbwwIntegrandBase* integrand, double& prob, double& probErr)
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

