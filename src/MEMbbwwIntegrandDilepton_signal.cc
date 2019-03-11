#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandDilepton_signal.h"

#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TMath.h>
#include <Math/VectorUtil.h> // ROOT::Math::VectorUtil::boost

using namespace mem;

MEMbbwwIntegrandDilepton_signal::MEMbbwwIntegrandDilepton_signal(double sqrtS, const std::string& madgraphFileName, int verbosity)
  : MEMbbwwIntegrandDilepton(sqrtS, madgraphFileName, verbosity)
  , applyOnshellWmassConstraint_(true)
  , offsetBJet1Theta_(-1)
  , offsetBJet1Phi_(-1)  
  , offsetBJet2Theta_(-1)
  , offsetBJet2Phi_(-1)
  , chargedLeptonPermutation_(kPermutationUndefined2L)
{
  if ( verbosity_ ) 
  {
    std::cout << "<MEMbbwwIntegrandDilepton_signal::MEMbbwwIntegrandDilepton_signal>:" << std::endl;
  }

  // initialize MadGraph
  if ( madgraphFileName != "" ) 
  {
    std::cout << "initializing MadGraph ME for HH signal using " << madgraphFileName << " file." << std::endl;
    me_madgraph_.initProc(madgraphFileName);
    madgraphIsInitialized_ = true;
  } 
  else 
  {
    std::cerr << "Error in <MEMbbwwIntegrandDilepton_signal>: No param.dat file for MadGraph given !!" << std::endl;
    assert(0);
  }

  // define ordering of four-vectors passed to MadGraph when evaluating the matrix element:
  //   1) first incoming gluon
  //   2) second incoming gluon
  //   3) first outgoing b quark (b-jet)
  //   4) second outgoing b quark (b-jet)
  //   5) outgoing lepton of positive charge
  //   6) outgoing neutrino
  //   7) outgoing lepton of negative charge
  //   8) outgoing anti-neutrino
  //
  // Note: the ordering needs to match the labels 1-8 as they are defined in the Feynman diagram doc/mg5/hh.ps
  madgraphMomenta_.push_back(madgraphGluon1P4_);
  madgraphMomenta_.push_back(madgraphGluon2P4_);
  madgraphMomenta_.push_back(madgraphBJet1P4_);
  madgraphMomenta_.push_back(madgraphBJet2P4_);
  madgraphMomenta_.push_back(madgraphChargedLeptonPlusP4_);
  madgraphMomenta_.push_back(madgraphNeutrinoP4_);
  madgraphMomenta_.push_back(madgraphChargedLeptonMinusP4_);
  madgraphMomenta_.push_back(madgraphAntiNeutrinoP4_);
}

MEMbbwwIntegrandDilepton_signal::~MEMbbwwIntegrandDilepton_signal()
{}

void 
MEMbbwwIntegrandDilepton_signal::initializeIntVars()
{
  intNumDimensions_ = 5; 
  if ( !applyOnshellWmassConstraint_ ) {
    intNumDimensions_ += 1;
  }
  if ( !measuredBJet1_ ) 
  {
    intNumDimensions_ += 2; 
  }
  if ( !measuredBJet2_ ) {
    intNumDimensions_ += 2; 
  }
  delete [] intBounds_lower_;
  intBounds_lower_ = new double[intNumDimensions_];
  delete [] intBounds_upper_;
  intBounds_upper_ = new double[intNumDimensions_];
  intVarNames_.clear();
  if ( measuredBJet2_ && !measuredBJet1_ ) 
  {
    intVarNames_.push_back("BJet2En"); 
  }
  else 
  {
    intVarNames_.push_back("BJet1En"); 
  }
  intBounds_lower_[0] = -1.; // to be set as function of measured b-jet energy and expected b-jet energy resolution
  intBounds_upper_[0] = -1.;
  intVarNames_.push_back("Nu1Theta"); 
  intBounds_lower_[1] = 0.;
  intBounds_upper_[1] = TMath::Pi();
  intVarNames_.push_back("Nu1Phi"); 
  intBounds_lower_[2] = -TMath::Pi();
  intBounds_upper_[2] = +TMath::Pi();
  intVarNames_.push_back("Nu2Theta"); 
  intBounds_lower_[3] = 0.;
  intBounds_upper_[3] = TMath::Pi();
  intVarNames_.push_back("Nu2Phi"); 
  intBounds_lower_[4] = -TMath::Pi();
  intBounds_upper_[4] = +TMath::Pi();  
  int offset = 5;     
  if ( !applyOnshellWmassConstraint_ ) 
  {
    intVarNames_.push_back("q2W"); 
    intBounds_lower_[5] = 0.;
    intBounds_upper_[5] = square(wBosonMass + 3.*wBosonWidth);
    offset += 1;
  }  
  if ( !measuredBJet1_ ) 
  {
    offsetBJet1Theta_ = offset; 
    intVarNames_.push_back("BJet1Theta"); 
    intBounds_lower_[offsetBJet1Theta_] = 0.;
    intBounds_upper_[offsetBJet1Theta_] = TMath::Pi();
    offsetBJet1Phi_ = offset + 1; 
    intVarNames_.push_back("BJet1Phi"); 
    intBounds_lower_[offsetBJet1Phi_] = -TMath::Pi();
    intBounds_upper_[offsetBJet1Phi_] = +TMath::Pi();
    offset += 2; 
  }
  if ( !measuredBJet2_ ) 
  {
    offsetBJet2Theta_ = offset; 
    intVarNames_.push_back("BJet2Theta"); 
    intBounds_lower_[offsetBJet2Theta_] = 0.;
    intBounds_upper_[offsetBJet2Theta_] = TMath::Pi();
    offsetBJet2Phi_ = offset + 1; 
    intVarNames_.push_back("BJet2Phi"); 
    intBounds_lower_[offsetBJet2Phi_] = -TMath::Pi();
    intBounds_upper_[offsetBJet2Phi_] = +TMath::Pi();
    offset += 2; 
  }
}

void 
MEMbbwwIntegrandDilepton_signal::applyOnshellWmassConstraint(bool flag) 
{ 
  applyOnshellWmassConstraint_ = flag; 
}

void 
MEMbbwwIntegrandDilepton_signal::setInputs(const MeasuredParticle* measuredChargedLeptonPlus, const MeasuredParticle* measuredChargedLeptonMinus,
					   const MeasuredParticle* measuredBJet1, const MeasuredParticle* measuredBJet2,
					   double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ ) {
    std::cout << "<MEMbbwwIntegrandDilepton_signal::setInputs>:" << std::endl;
  }

  MEMbbwwIntegrandDilepton::setInputs(
    measuredChargedLeptonPlus, measuredChargedLeptonMinus,
    measuredBJet1, measuredBJet2,
    measuredMEtPx, measuredMEtPy, measuredMEtCov);
  
  /// define integration variables
  initializeIntVars();

  // set integration boundary for energy of first b-jet
  // to measured jet energy +/- 3 times expected energy resolution (rough estimate),
  // in order to increase efficiency of numeric integration
  double measuredBJetEn, measuredBJetEnRes;
  if ( measuredBJet2_ && !measuredBJet1_ ) 
  {
    assert(measuredBJet2_);
    measuredBJetEn = measuredBJet2_->energy();
    measuredBJetEnRes = 0.50*TMath::Sqrt(measuredBJet2_->energy());
  }
  else 
  {
    assert(measuredBJet1_);
    measuredBJetEn = measuredBJet1_->energy();
    measuredBJetEnRes = 0.50*TMath::Sqrt(measuredBJet1_->energy());
  }
  intBounds_lower_[0] = TMath::Max(10., measuredBJetEn - 3.*measuredBJetEnRes);
  intBounds_upper_[0] = measuredBJetEn + 3.*measuredBJetEnRes;

  // Standard Model (SM) cross section for (non-resonant) HH production @ 13 TeV center-of-mass energy
  // time branching fraction for the decay HH->bbWW->bblnulnu (excluding electrons and muons from tau decays)
  //
  // Note: SM cross section is taken from next-to-next-to-leading order (NNLO) computation, 
  //       while MadGraph matrix element is leading order (LO). 
  //       We expect this inconsistency to have little practical effect.
  const double crossSection_signal = 33.53e-3*2.*0.577*0.215*square(0.216); // [pb]

  // compute product of constant terms in the integrand (terms that do not depend on the integration variables),
  // in order to reduce computing time required for numeric integration
  double numerator = square(higgsBosonMass*higgsBosonWidth)*wBosonMass*wBosonWidth;
  double denominator = 1.;
  if ( applyOnshellWmassConstraint_ ) {
    denominator *= (TMath::Pi()*wBosonMass*wBosonWidth);
  }
  denominator *= (TMath::Power(2., 22)*TMath::Power(TMath::Pi(), 15));
  denominator *= measuredChargedLeptonPlus_->energy();
  denominator *= measuredChargedLeptonMinus_->energy();
  if ( measuredBJet1_ ) 
  {
    denominator *= (measuredBJet1_->p()*measuredBJet1_->energy());
  }
  if ( measuredBJet2_ ) 
  {
    denominator *= (measuredBJet2_->p()*measuredBJet2_->energy());
  }
  denominator *= (crossSection_signal*square(sqrtS_));
  assert(denominator > 0.);
  normFactor_ = numerator/denominator;
}

void 
MEMbbwwIntegrandDilepton_signal::setOnshellChargedLepton(int chargedLeptonPermutation)
{
  chargedLeptonPermutation_ = chargedLeptonPermutation;
}

double MEMbbwwIntegrandDilepton_signal::Eval(const double* x) const
{
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "<MEMbbwwIntegrandDilepton_signal::Eval>:" << std::endl;
    printVDouble(" x", x, intNumDimensions_);
  }

  assert(chargedLeptonPermutation_ == kOnshellChargedLeptonPlus || chargedLeptonPermutation_ == kOnshellChargedLeptonMinus);

  double trueBJet1Theta, trueBJet1Phi;
  if ( measuredBJet1_ ) 
  {
    trueBJet1Theta = measuredBJet1_->theta();
    trueBJet1Phi = measuredBJet1_->phi();
  } 
  else 
  {
    trueBJet1Theta = x[offsetBJet1Theta_];
    trueBJet1Phi = x[offsetBJet1Phi_];
  }  
  double trueBJet2Theta, trueBJet2Phi;
  if ( measuredBJet2_ ) 
  {
    trueBJet2Theta = measuredBJet2_->theta();
    trueBJet2Phi = measuredBJet2_->phi();
  } 
  else 
  {
    trueBJet2Theta = x[offsetBJet2Theta_];
    trueBJet2Phi = x[offsetBJet2Phi_];
  }  
  LorentzVector trueBJet1P4, trueBJet2P4;
  if ( measuredBJet2_ && !measuredBJet1_ ) 
  {
    double trueBJet2En = x[0];
    trueBJet2P4 = buildLorentzVector(trueBJet2En, trueBJet2Theta, trueBJet2Phi, bottomQuarkMass);
    std::vector<double> trueBJet1En_solutions = compBJetEn_Hbb(trueBJet2P4, trueBJet1Theta, trueBJet1Phi);
    bool trueBJet1En_foundSolution = false;
    trueBJet1P4 = findBJetEn_solution_Hbb(trueBJet1En_solutions, trueBJet1Theta, trueBJet1Phi, trueBJet2P4, trueBJet1En_foundSolution);
    if ( !trueBJet1En_foundSolution ) return 0.;
  }
  else 
  {
    double trueBJet1En = x[0];
    trueBJet1P4 = buildLorentzVector(trueBJet1En, trueBJet1Theta, trueBJet1Phi, bottomQuarkMass);
    std::vector<double> trueBJet2En_solutions = compBJetEn_Hbb(trueBJet1P4, trueBJet2Theta, trueBJet2Phi);
    bool trueBJet2En_foundSolution = false;
    trueBJet2P4 = findBJetEn_solution_Hbb(trueBJet2En_solutions, trueBJet2Theta, trueBJet2Phi, trueBJet1P4, trueBJet2En_foundSolution);
    if ( !trueBJet2En_foundSolution ) return 0.;
  }

  const LorentzVector& trueChargedLeptonPlusP4 = measuredChargedLeptonPlus_->p4();
  const LorentzVector& trueChargedLeptonMinusP4 = measuredChargedLeptonMinus_->p4();

  double trueNuTheta = x[1];
  double trueNuPhi = x[2];
  LorentzVector trueNuP4;
  double trueAntiNuTheta = x[3];
  double trueAntiNuPhi = x[4];
  LorentzVector trueAntiNuP4;
  if ( chargedLeptonPermutation_ == kOnshellChargedLeptonPlus ) 
  { 
    // lepton of positive charge (and hence neutrino) originates from on-shell W boson,
    // while lepton of negative charge (and hence anti-neutrino) originates from off-shell W boson
    //
    // Note: the energy (and four-vector) of the neutrino originating from the on-shell W boson has to be computed first!
    double trueNuEn;
    if ( applyOnshellWmassConstraint_ ) 
    {
      trueNuEn = compNuEn_Wlnu(trueChargedLeptonPlusP4, trueNuTheta, trueNuPhi);
    } 
    else 
    {
      double q2W = x[5];
      trueNuEn = compNuEn_Wlnu(trueChargedLeptonPlusP4, trueNuTheta, trueNuPhi, q2W);
    }
    if ( !(trueNuEn > 0.) ) return 0.;
    trueNuP4 = buildLorentzVector(trueNuEn, trueNuTheta, trueNuPhi);
    LorentzVector trueChargedLeptonPlus_Nu_ChargedLeptonMinusP4 = trueChargedLeptonPlusP4 + trueNuP4 + trueChargedLeptonMinusP4;
    if ( trueChargedLeptonPlus_Nu_ChargedLeptonMinusP4.mass() >= higgsBosonMass ) return 0.;
    double trueAntiNuEn = compNuStarEn_Hww(trueChargedLeptonPlus_Nu_ChargedLeptonMinusP4, trueAntiNuTheta, trueAntiNuPhi);
    if ( !(trueAntiNuEn > 0.) ) return 0.;
    trueAntiNuP4 = buildLorentzVector(trueAntiNuEn, trueAntiNuTheta, trueAntiNuPhi);
  } 
  else 
  {
    // lepton of negative charge (and hence anti-neutrino) originates from on-shell W boson,
    // while lepton of positive charge (and hence neutrino) originates from off-shell W boson
    //
    // Note: the energy (and four-vector) of the anti-neutrino originating from the on-shell W boson has to be computed first!
    double trueAntiNuEn;
    if ( applyOnshellWmassConstraint_ ) 
    {
      trueAntiNuEn = compNuEn_Wlnu(trueChargedLeptonMinusP4, trueAntiNuTheta, trueAntiNuPhi);
    } 
    else 
    {
      double q2W = x[5];
      trueAntiNuEn = compNuEn_Wlnu(trueChargedLeptonMinusP4, trueAntiNuTheta, trueAntiNuPhi, q2W);
    }
    if ( !(trueAntiNuEn > 0.) ) return 0.;
    trueAntiNuP4 = buildLorentzVector(trueAntiNuEn, trueAntiNuTheta, trueAntiNuPhi);
    LorentzVector trueChargedLeptonMinus_AntiNu_ChargedLeptonPlusP4 = trueChargedLeptonMinusP4 + trueAntiNuP4 + trueChargedLeptonPlusP4;
    if ( trueChargedLeptonMinus_AntiNu_ChargedLeptonPlusP4.mass() >= higgsBosonMass ) return 0.;
    double trueNuEn = compNuStarEn_Hww(trueChargedLeptonMinus_AntiNu_ChargedLeptonPlusP4, trueNuTheta, trueNuPhi);
    if ( !(trueNuEn > 0.) ) return 0.;
    trueNuP4 = buildLorentzVector(trueNuEn, trueNuTheta, trueNuPhi);
  }

  LorentzVector trueSumP4 = trueBJet1P4 + trueBJet2P4 + trueChargedLeptonPlusP4 + trueNuP4 + trueChargedLeptonMinusP4 + trueAntiNuP4;

  // perform boost into zero-transverse-momentum (ZTM) frame 
  double ztmFramePx = trueSumP4.px();
  double ztmFramePy = trueSumP4.py();
  double ztmFrameEn = trueSumP4.energy();
  Vector boost(-ztmFramePx/ztmFrameEn, -ztmFramePy/ztmFrameEn, 0.);
  LorentzVector trueBJet1P4_ztm = ROOT::Math::VectorUtil::boost(trueBJet1P4, boost); 
  LorentzVector trueBJet2P4_ztm = ROOT::Math::VectorUtil::boost(trueBJet2P4, boost); 
  LorentzVector trueChargedLeptonPlusP4_ztm = ROOT::Math::VectorUtil::boost(trueChargedLeptonPlusP4, boost);
  LorentzVector trueNuP4_ztm = ROOT::Math::VectorUtil::boost(trueNuP4, boost);
  LorentzVector trueChargedLeptonMinusP4_ztm = ROOT::Math::VectorUtil::boost(trueChargedLeptonMinusP4, boost);
  LorentzVector trueAntiNuP4_ztm = ROOT::Math::VectorUtil::boost(trueAntiNuP4, boost);
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "laboratory frame:" << std::endl;
    if ( measuredBJet1_ ) 
    {
      printLorentzVector("b-jet1", trueBJet1P4, measuredBJet1_->p4());
    } 
    else 
    {
      printLorentzVector_NA("b-jet1", trueBJet1P4);
    }
    if ( measuredBJet2_ ) 
    {
      printLorentzVector("b-jet2", trueBJet2P4, measuredBJet2_->p4());
    } 
    else 
    {
      printLorentzVector_NA("b-jet2", trueBJet2P4);
    }
    printLorentzVector("lepton+", trueChargedLeptonPlusP4, measuredChargedLeptonPlus_->p4());
    printLorentzVector("neutrino", trueNuP4);
    printLorentzVector("lepton-", trueChargedLeptonMinusP4, measuredChargedLeptonMinus_->p4());
    printLorentzVector("anti-neutrino", trueAntiNuP4);
    std::cout << "m(b bbar) = " << (trueBJet1P4 + trueBJet2P4).mass() << std::endl;
    std::cout << "m(lep+ nu) = " << (trueChargedLeptonPlusP4 + trueNuP4).mass() << std::endl;
    std::cout << "m(lep- nu) = " << (trueChargedLeptonMinusP4 + trueAntiNuP4).mass() << std::endl;
    std::cout << "m(lep+ nu lep- nu) = " << (trueChargedLeptonPlusP4 + trueNuP4 + trueChargedLeptonMinusP4 + trueAntiNuP4).mass() << std::endl;
    printLorentzVector("sum", trueSumP4);
    std::cout << "zero-transverse-momentum frame:" << std::endl;
    printLorentzVector("b-jet1", trueBJet1P4_ztm);
    printLorentzVector("b-jet2", trueBJet2P4_ztm);
    printLorentzVector("lepton+", trueChargedLeptonPlusP4_ztm);
    printLorentzVector("neutrino", trueNuP4_ztm);
    printLorentzVector("lepton-", trueChargedLeptonMinusP4_ztm);
    printLorentzVector("anti-neutrino", trueAntiNuP4_ztm);
    std::cout << "m(b bbar) = " << (trueBJet1P4_ztm + trueBJet2P4_ztm).mass() << std::endl;
    std::cout << "m(lep+ nu) = " << (trueChargedLeptonPlusP4_ztm + trueNuP4_ztm).mass() << std::endl;
    std::cout << "m(lep- nu) = " << (trueChargedLeptonMinusP4_ztm + trueAntiNuP4_ztm).mass() << std::endl;
    std::cout << "m(lep+ nu lep- nu) = " << (trueChargedLeptonPlusP4_ztm + trueNuP4_ztm + trueChargedLeptonMinusP4_ztm + trueAntiNuP4_ztm).mass() << std::endl;
    LorentzVector trueSumP4_ztm = trueBJet1P4_ztm + trueBJet2P4_ztm + trueChargedLeptonPlusP4_ztm + trueNuP4_ztm + trueChargedLeptonMinusP4_ztm + trueAntiNuP4_ztm;
    printLorentzVector("sum", trueSumP4_ztm);
  }

  // compute Bjorken-x of incoming protons and evaluate PDF factor
  double trueSumPz = trueSumP4.pz();
  double trueSumEn = trueSumP4.E();
  double xa = (trueSumEn + trueSumPz)/sqrtS_;
  double xb = (trueSumEn - trueSumPz)/sqrtS_;
  if ( xa <= 0. || xa >= 1. ) return 0.;
  if ( xb <= 0. || xb >= 1. ) return 0.;
  double Q = 0.5*trueSumP4.mass();
  assert(pdf_);
  double fa = pdf_->xfxQ(21, xa, Q)/xa; // gluon distribution
  double fb = pdf_->xfxQ(21, xb, Q)/xb;
  double prob_PDF = (fa*fb);

  // evaluate flux factor 
  //
  // Note: the factor 1/s is aready included in the "normFactor" data-member
  double prob_flux = (1./(xa*xb)); 

  // evaluate LO matrix element, generated using MadGraph
  madgraphGluon1P4_[0] =  0.5*xa*sqrtS_; 
  madgraphGluon1P4_[3] = +0.5*xa*sqrtS_;
  madgraphGluon2P4_[0] =  0.5*xb*sqrtS_;
  madgraphGluon2P4_[3] = -0.5*xb*sqrtS_;
  madgraphChargedLeptonPlusP4_[0] = trueChargedLeptonPlusP4_ztm.energy();
  madgraphChargedLeptonPlusP4_[1] = trueChargedLeptonPlusP4_ztm.px();
  madgraphChargedLeptonPlusP4_[2] = trueChargedLeptonPlusP4_ztm.py();
  madgraphChargedLeptonPlusP4_[3] = trueChargedLeptonPlusP4_ztm.pz();
  madgraphNeutrinoP4_[0] = trueNuP4_ztm.energy();
  madgraphNeutrinoP4_[1] = trueNuP4_ztm.px();
  madgraphNeutrinoP4_[2] = trueNuP4_ztm.py();
  madgraphNeutrinoP4_[3] = trueNuP4_ztm.pz();
  madgraphChargedLeptonMinusP4_[0] = trueChargedLeptonMinusP4_ztm.energy();
  madgraphChargedLeptonMinusP4_[1] = trueChargedLeptonMinusP4_ztm.px();
  madgraphChargedLeptonMinusP4_[2] = trueChargedLeptonMinusP4_ztm.py();
  madgraphChargedLeptonMinusP4_[3] = trueChargedLeptonMinusP4_ztm.pz();
  madgraphAntiNeutrinoP4_[0] = trueAntiNuP4_ztm.energy();
  madgraphAntiNeutrinoP4_[1] = trueAntiNuP4_ztm.px();
  madgraphAntiNeutrinoP4_[2] = trueAntiNuP4_ztm.py();
  madgraphAntiNeutrinoP4_[3] = trueAntiNuP4_ztm.pz();
  madgraphBJet1P4_[0] = trueBJet1P4_ztm.energy();
  madgraphBJet1P4_[1] = trueBJet1P4_ztm.px();
  madgraphBJet1P4_[2] = trueBJet1P4_ztm.py();
  madgraphBJet1P4_[3] = trueBJet1P4_ztm.pz();
  madgraphBJet2P4_[0] = trueBJet2P4_ztm.energy();
  madgraphBJet2P4_[1] = trueBJet2P4_ztm.px();
  madgraphBJet2P4_[2] = trueBJet2P4_ztm.py();
  madgraphBJet2P4_[3] = trueBJet2P4_ztm.pz();
  if ( verbosity_ >= 2 )
  {
    printMadGraphMomenta(madgraphMomenta_);
  }
  double prob_ME = -1.;
  if ( madgraphIsInitialized_ ) {
    me_madgraph_.setHiggsWidth(2.); // CV: enlarge Higgs boson width to make ME evaluation robust against rounding errors
    me_madgraph_.setMomenta(madgraphMomenta_);
    me_madgraph_.sigmaKin();
    prob_ME = me_madgraph_.getMatrixElements()[0];
    ++numMatrixElementEvaluations_;
    if ( TMath::IsNaN(prob_ME) ) 
    {
      std::cerr << "Warning: MadGraph returned NaN --> skipping event !!" << std::endl;
      printLorentzVector("b-jet1", trueBJet1P4_ztm);
      printLorentzVector("b-jet2", trueBJet2P4_ztm);
      printLorentzVector("lepton+", trueChargedLeptonPlusP4_ztm);
      printLorentzVector("neutrino", trueNuP4_ztm);
      printLorentzVector("lepton-", trueChargedLeptonMinusP4_ztm);
      printLorentzVector("anti-neutrino", trueAntiNuP4_ztm);
      return 0.;
    }
  }
  assert(prob_ME >= 0.);
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "prob_ME = " << prob_ME << std::endl;
  }

  // CV: do not include "missing", i.e. non-reconstructed, b-jets in trueHadRecoil
  //    (for consistency with computation of measuredHadRecoil in MEMbbwwIntegrandDilepton::setInputs function)
  double trueHadRecoilPx = -(trueChargedLeptonPlusP4.px() + trueChargedLeptonMinusP4.px() + trueNuP4.px() + trueAntiNuP4.px());
  double trueHadRecoilPy = -(trueChargedLeptonPlusP4.py() + trueChargedLeptonMinusP4.py() + trueNuP4.py() + trueAntiNuP4.py());
  if ( measuredBJet1_ ) 
  {
    trueHadRecoilPx -= trueBJet1P4.px();
    trueHadRecoilPy -= trueBJet1P4.py();
  }
  if ( measuredBJet2_ ) 
  {
    trueHadRecoilPx -= trueBJet2P4.px();
    trueHadRecoilPy -= trueBJet2P4.py();
  }
  if ( verbosity_ >= 2 )
  {
    std::cout << "hadRecoil:" << std::endl;
    std::cout << " true Px = " << trueHadRecoilPx << ", Py = " << trueHadRecoilPy << std::endl;
    std::cout << " rec. Px = " << measuredHadRecoilPx_ << ", Py = " << measuredHadRecoilPy_ << std::endl;
    if ( !measuredBJet1_ || !measuredBJet2_ ) 
    {
      std::cout << "Note:";
      if ( !measuredBJet1_ ) std::cout << " b-jet1 is 'missing', i.e. not reconstructed, and not included in hadRecoil.";
      if ( !measuredBJet2_ ) std::cout << " b-jet2 is 'missing', i.e. not reconstructed, and not included in hadRecoil.";
      std::cout << std::endl;
    }
  }

  double prob_TF = 1.;
  if ( measuredBJet1_ ) 
  {
    prob_TF *= bjet1TF_->Eval(trueBJet1P4.energy());
  }
  if ( measuredBJet2_ ) 
  {
    prob_TF *= bjet2TF_->Eval(trueBJet2P4.energy());
  }
  prob_TF *= hadRecoilTF_->Eval(trueHadRecoilPx, trueHadRecoilPy);
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "prob_TF = " << prob_TF << std::endl;
  }

  double jacobiFactor = compJacobiFactor_Hbb(trueBJet1P4, trueBJet2P4);
  if ( chargedLeptonPermutation_ == kOnshellChargedLeptonPlus ) { 
    //---------------------------------------------------------------------------
    // CV: Jacobi factor is the same for case that W mass constraint is applied to on-shell W boson and for the case that it is not applied
    jacobiFactor *= compJacobiFactor_Wlnu(trueChargedLeptonPlusP4, trueNuP4);
    //---------------------------------------------------------------------------
    LorentzVector trueChargedLeptonPlus_Nu_ChargedLeptonMinusP4 = trueChargedLeptonPlusP4 + trueNuP4 + trueChargedLeptonMinusP4;
    jacobiFactor *= compJacobiFactor_Hww(trueChargedLeptonPlus_Nu_ChargedLeptonMinusP4, trueAntiNuP4);
  } else {
    //---------------------------------------------------------------------------
    // CV: Jacobi factor is the same for case that W mass constraint is applied to on-shell W boson and for the case that it is not applied
    jacobiFactor *= compJacobiFactor_Wlnu(trueChargedLeptonMinusP4, trueAntiNuP4);
    //---------------------------------------------------------------------------
    LorentzVector trueChargedLeptonMinus_AntiNu_ChargedLeptonPlusP4 = trueChargedLeptonMinusP4 + trueAntiNuP4 + trueChargedLeptonPlusP4;
    jacobiFactor *= compJacobiFactor_Hww(trueChargedLeptonMinus_AntiNu_ChargedLeptonPlusP4, trueNuP4);
  }  
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "jacobiFactor = " << jacobiFactor << std::endl;
  }

  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  double integrandValue = conversionFactor*normFactor_;
  const double fudgeFactor = 0.05*1.6e+29;
  const double fudgeFactor_missingBJet = 0.02*2.9e+24;
  if ( measuredBJet1_ && measuredBJet2_ ) integrandValue *= fudgeFactor;
  else integrandValue *= fudgeFactor_missingBJet;
  integrandValue *= (trueNuP4.pt()*trueAntiNuP4.pt());
  integrandValue *= prob_PDF;
  integrandValue *= prob_flux;
  integrandValue *= prob_ME;
  if ( measuredBJet1_ ) 
  {
    integrandValue *= trueBJet1P4.P();
  }
  else 
  {
    integrandValue *= (1./trueBJet1P4.energy());    
  }
  if ( measuredBJet2_ ) 
  {
    integrandValue *= trueBJet2P4.P();
  }
  else 
  {
    integrandValue *= (1./trueBJet2P4.energy());    
  }
  integrandValue *= prob_TF;
  integrandValue *= jacobiFactor;
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "--> integrandValue = " << integrandValue << std::endl;
  }

  return integrandValue;
}

