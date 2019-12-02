#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandSingleLepton_signal.h"

#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TMath.h>
#include <Math/VectorUtil.h> // ROOT::Math::VectorUtil::boost

using namespace mem;

MEMbbwwIntegrandSingleLepton_signal::MEMbbwwIntegrandSingleLepton_signal(double sqrtS, const std::string& madgraphFileName, int verbosity)
  : MEMbbwwIntegrandSingleLepton(sqrtS, madgraphFileName, verbosity)
  , applyOnshellWmassConstraint_(true)
  , offsetHadWJet1Theta_(-1)
  , offsetHadWJet1Phi_(-1)  
  , offsetHadWJet2Theta_(-1)
  , offsetHadWJet2Phi_(-1)	
  , offsetBJet1Theta_(-1)
  , offsetBJet1Phi_(-1)  
  , offsetBJet2Theta_(-1)
  , offsetBJet2Phi_(-1)
  , chargedLeptonPermutation_(kPermutationUndefined1L)
{
  if ( verbosity_ ) 
  {
    std::cout << "<MEMbbwwIntegrandSingleLepton_signal::MEMbbwwIntegrandSingleLepton_signal>:" << std::endl;
  }

  /// define integration variables
  initializeIntVars();

  // initialize MadGraph
  if ( madgraphFileName != "" ) 
  {
    std::cout << "initializing MadGraph ME for HH signal using " << madgraphFileName << " file." << std::endl;
    me_madgraph_chargedLeptonPlus_.initProc(madgraphFileName);
    me_madgraph_chargedLeptonMinus_.initProc(madgraphFileName);
    madgraphIsInitialized_ = true;
  } 
  else 
  {
    std::cerr << "Error in <MEMbbwwIntegrandSingleLepton_signal>: No param.dat file for MadGraph given !!" << std::endl;
    assert(0);
  }

  // define ordering of four-vectors passed to MadGraph when evaluating the matrix element
  // for events containing leptons of positive charge:
  //   1) first incoming gluon
  //   2) second incoming gluon
  //   3) first outgoing b quark (b-jet)
  //   4) second outgoing b quark (b-jet)
  //   5) outgoing charged lepton
  //   6) outgoing neutrino
  //   7) first outgoing jet from W->jj decay
  //   8) second outgoing jet from W->jj decay
  //
  // Note: the ordering needs to match the labels 1-8 as they are defined in the Feynman diagram doc/mg5/hh_sl_pos.ps
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphGluon1P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphGluon2P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphBJet1P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphBJet2P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphChargedLeptonP4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphNeutrinoP4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphHadWJet1P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphHadWJet2P4_);

  // define ordering of four-vectors passed to MadGraph when evaluating the matrix element
  // for events containing leptons of negative charge:
  //   1) first incoming gluon
  //   2) second incoming gluon
  //   3) first outgoing b quark (b-jet)
  //   4) second outgoing b quark (b-jet)
  //   5) first outgoing jet from W->jj decay
  //   6) second outgoing jet from W->jj decay
  //   7) outgoing charged lepton
  //   8) outgoing neutrino
  //
  // Note: the ordering needs to match the labels 1-8 as they are defined in the Feynman diagram doc/mg5/hh_sl_neg.ps
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphGluon1P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphGluon2P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphBJet1P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphBJet2P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphHadWJet1P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphHadWJet2P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphChargedLeptonP4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphNeutrinoP4_);
}

MEMbbwwIntegrandSingleLepton_signal::~MEMbbwwIntegrandSingleLepton_signal()
{}

void 
MEMbbwwIntegrandSingleLepton_signal::initializeIntVars()
{
  intNumDimensions_ = 4;
  if ( !applyOnshellWmassConstraint_ ) 
  {
    intNumDimensions_ += 1;
  }
  if ( !measuredHadWJet1_ ) 
  {
    intNumDimensions_ += 2; 
  }
  if ( !measuredHadWJet2_ ) 
  {
    intNumDimensions_ += 2; 
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
  intVarNames_.push_back("NuTheta"); 
  intBounds_lower_[1] = 0.;
  intBounds_upper_[1] = TMath::Pi();
  intVarNames_.push_back("NuPhi"); 
  intBounds_lower_[2] = -TMath::Pi();
  intBounds_upper_[2] = +TMath::Pi();
  if ( measuredHadWJet2_ && !measuredHadWJet1_ ) 
  {
    intVarNames_.push_back("HadWJet2En"); 
  }
  else 
  {	
    intVarNames_.push_back("HadWJet1En"); 
  }
  intBounds_lower_[3] = -1.; // to be set as function of measured jet energy and expected jet energy resolution
  intBounds_upper_[3] = -1.;
  int offset = 4;    
  if ( !applyOnshellWmassConstraint_ ) {
    intVarNames_.push_back("q2W"); 
    intBounds_lower_[4] = 0.;
    intBounds_upper_[4] = square(wBosonMass + 3.*wBosonWidth);
    offset += 1;
  }
  if ( !measuredHadWJet1_ ) 
  {
    offsetHadWJet1Theta_ = offset; 
    intVarNames_.push_back("HadWJet1Theta"); 
    intBounds_lower_[offsetHadWJet1Theta_] = 0.;
    intBounds_upper_[offsetHadWJet1Theta_] = TMath::Pi();
    offsetHadWJet1Phi_ = offset + 1; 
    intVarNames_.push_back("HadWJet1Phi"); 
    intBounds_lower_[offsetHadWJet1Phi_] = -TMath::Pi();
    intBounds_upper_[offsetHadWJet1Phi_] = +TMath::Pi();
    offset += 2; 
  }
  if ( !measuredHadWJet2_ ) 
  {
    offsetHadWJet2Theta_ = offset; 
    intVarNames_.push_back("HadWJet2Theta"); 
    intBounds_lower_[offsetHadWJet2Theta_] = 0.;
    intBounds_upper_[offsetHadWJet2Theta_] = TMath::Pi();
    offsetBJet2Phi_ = offset + 1; 
    intVarNames_.push_back("HadWJet2Phi"); 
    intBounds_lower_[offsetHadWJet2Phi_] = -TMath::Pi();
    intBounds_upper_[offsetHadWJet2Phi_] = +TMath::Pi();
    offset += 2; 
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
MEMbbwwIntegrandSingleLepton_signal::applyOnshellWmassConstraint(bool flag) 
{ 
  applyOnshellWmassConstraint_ = flag; 
}

void 
MEMbbwwIntegrandSingleLepton_signal::setInputs(const MeasuredParticle* measuredChargedLepton, 
					       const MeasuredParticle* measuredHadWJet1, const MeasuredParticle* measuredHadWJet2,
					       const MeasuredParticle* measuredBJet1, const MeasuredParticle* measuredBJet2,
					       double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ ) 
  {
    std::cout << "<MEMbbwwIntegrandSingleLepton_signal::setInputs>:" << std::endl;
  }

  MEMbbwwIntegrandSingleLepton::setInputs(
    measuredChargedLepton, 
    measuredHadWJet1, measuredHadWJet2,
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

  // set integration boundary for energy of first jet from W->jj decay
  // to measured jet energy +/- 3 times expected energy resolution (rough estimate),
  // in order to increase efficiency of numeric integration
  double measuredHadWJetEn, measuredHadWJetEnRes;
  if ( measuredHadWJet2_ && !measuredHadWJet1_ ) 
  {
    assert(measuredHadWJet2_);
    measuredHadWJetEn = measuredHadWJet2_->energy();
    measuredHadWJetEnRes = 0.50*TMath::Sqrt(measuredHadWJet2_->energy());
  }
  else 
  {
    assert(measuredHadWJet1_);
    measuredHadWJetEn = measuredHadWJet1_->energy();
    measuredHadWJetEnRes = 0.50*TMath::Sqrt(measuredHadWJet1_->energy());
  }	
  intBounds_lower_[3] = TMath::Max(10., measuredHadWJetEn - 3.*measuredHadWJetEnRes);
  intBounds_upper_[3] = measuredHadWJetEn + 3.*measuredHadWJetEnRes;
  
  // Standard Model (SM) cross section for (non-resonant) HH production @ 13 TeV center-of-mass energy
  // time branching fraction for the decay HH->bbWW->bblnulnu (excluding electrons and muons from tau decays)
  //
  // Note: SM cross section is taken from next-to-next-to-leading order (NNLO) computation, 
  //       while MadGraph matrix element is leading order (LO). 
  //       We expect this inconsistency to have little practical effect.
  const double crossSection_signal = 33.53e-3*2*0.577*0.215*0.216*0.6742; // [pb]

  // compute product of constant terms in the integrand (terms that do not depend on the integration variables),
  // in order to reduce computing time required for numeric integration
  double numerator = square(higgsBosonMass*higgsBosonWidth)*wBosonMass*wBosonWidth;
  double denominator = 1.;
  if ( applyOnshellWmassConstraint_ ) 
  {
    denominator *= (TMath::Pi()*wBosonMass*wBosonWidth);
  }
  denominator *= (TMath::Power(2., 25)*TMath::Power(TMath::Pi(), 18));
  denominator *= measuredChargedLepton_->energy();
  if ( measuredHadWJet1_ ) 
  {
    denominator *= (measuredHadWJet1_->p()*measuredHadWJet1_->energy());
  }
  if ( measuredHadWJet2_ ) 
  {
    denominator *= (measuredHadWJet2_->p()*measuredHadWJet2_->energy());
  }
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
MEMbbwwIntegrandSingleLepton_signal::setOnshellChargedLepton(int chargedLeptonPermutation)
{
  chargedLeptonPermutation_ = chargedLeptonPermutation;
}

double MEMbbwwIntegrandSingleLepton_signal::Eval(const double* x) const
{
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "<MEMbbwwIntegrandSingleLepton_signal::Eval>:" << std::endl;
    printVDouble(" x", x, intNumDimensions_);
  }

  assert(chargedLeptonPermutation_ == kOnshellChargedLepton || chargedLeptonPermutation_ == kOffshellChargedLepton);

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

  const LorentzVector& trueChargedLeptonP4 = measuredChargedLepton_->p4();

  double trueNuTheta = x[1];
  double trueNuPhi = x[2];
  LorentzVector trueNuP4;
  double trueHadWJet1Theta, trueHadWJet1Phi;
  if ( measuredHadWJet1_ ) 
  {
    trueHadWJet1Theta = measuredHadWJet1_->theta();
    trueHadWJet1Phi = measuredHadWJet1_->phi();
  } 
  else 
  {
    trueHadWJet1Theta = x[offsetHadWJet1Theta_];
    trueHadWJet1Phi = x[offsetHadWJet1Phi_];
  }  	
  double trueHadWJet2Theta, trueHadWJet2Phi;
  if ( measuredHadWJet2_ ) 
  {
    trueHadWJet2Theta = measuredHadWJet2_->theta();
    trueHadWJet2Phi = measuredHadWJet2_->phi();
  } 
  else 
  {
    trueHadWJet2Theta = x[offsetHadWJet2Theta_];
    trueHadWJet2Phi = x[offsetHadWJet2Phi_];
  }	
  LorentzVector trueHadWJet1P4, trueHadWJet2P4;
  if ( chargedLeptonPermutation_ == kOnshellChargedLepton ) 
  { 
    // charged lepton and neutrino originate from on-shell W boson,
    // while the two jets from W->jj decay originate from off-shell W boson
    //
    // Note: the energy (and four-vector) of the neutrino originating from the on-shell W boson has to be computed first!
    double trueNuEn;
    if ( applyOnshellWmassConstraint_ ) 
    {
      trueNuEn = compNuEn_Wlnu(trueChargedLeptonP4, trueNuTheta, trueNuPhi);
    } 
    else 
    {
      double q2W = x[4];
      trueNuEn = compNuEn_Wlnu(trueChargedLeptonP4, trueNuTheta, trueNuPhi, q2W);
    }
    if ( !(trueNuEn > 0.) ) return 0.;
    trueNuP4 = buildLorentzVector(trueNuEn, trueNuTheta, trueNuPhi);
    if ( measuredBJet2_ && !measuredBJet1_ ) 
    {
      double trueHadWJet2En = x[3]; 
      trueHadWJet2P4 = buildLorentzVector(trueHadWJet2En, trueHadWJet2Theta, trueHadWJet2Phi);
      LorentzVector trueChargedLepton_Nu_HadWJet2P4 = trueChargedLeptonP4 + trueNuP4 + trueHadWJet2P4;
      if ( trueChargedLepton_Nu_HadWJet2P4.mass() >= higgsBosonMass ) return 0.;
      double trueHadWJet1En = compHadWJet1En_Hww(trueChargedLepton_Nu_HadWJet2P4, trueHadWJet1Theta, trueHadWJet1Phi);
      if ( !(trueHadWJet1En > 0.) ) return 0.;
      trueHadWJet1P4 = buildLorentzVector(trueHadWJet1En, trueHadWJet1Theta, trueHadWJet1Phi);
    } 
    else
    {
      double trueHadWJet1En = x[3];
      trueHadWJet1P4 = buildLorentzVector(trueHadWJet1En, trueHadWJet1Theta, trueHadWJet1Phi);
      LorentzVector trueChargedLepton_Nu_HadWJet1P4 = trueChargedLeptonP4 + trueNuP4 + trueHadWJet1P4;
      if ( trueChargedLepton_Nu_HadWJet1P4.mass() >= higgsBosonMass ) return 0.;
      double trueHadWJet2En = compHadWJet2En_Hww(trueChargedLepton_Nu_HadWJet1P4, trueHadWJet2Theta, trueHadWJet2Phi);
      if ( !(trueHadWJet2En > 0.) ) return 0.;
      trueHadWJet2P4 = buildLorentzVector(trueHadWJet2En, trueHadWJet2Theta, trueHadWJet2Phi);
    }
  } 
  else 
  {
    // the two jets from W->jj decay originate from on-shell W boson,
    // while charged lepton and neutrino originate from off-shell W boson
    //
    // Note: the energy (and four-vector) of the two jets from W->jj decay originating from the on-shell W boson has to be computed first!
    if ( measuredBJet2_ && !measuredBJet1_ ) 
    {	
      double trueHadWJet2En = x[3];
      trueHadWJet2P4 = buildLorentzVector(trueHadWJet2En, trueHadWJet2Theta, trueHadWJet2Phi);	
      double trueHadWJet1En;
      if ( applyOnshellWmassConstraint_ ) 
      {
        trueHadWJet1En = compHadWJet1En_Wjj(trueHadWJet2P4, trueHadWJet1Theta, trueHadWJet1Phi);
      } 
      else 
      {
        double q2W = x[4];
        trueHadWJet1En = compHadWJet1En_Wjj(trueHadWJet2P4, trueHadWJet1Theta, trueHadWJet1Phi, q2W);
      }
      if ( !(trueHadWJet1En > 0.) ) return 0.;
      trueHadWJet1P4 = buildLorentzVector(trueHadWJet1En, trueHadWJet1Theta, trueHadWJet1Phi);
    }
    else 
    { 
      double trueHadWJet1En = x[3];
      trueHadWJet1P4 = buildLorentzVector(trueHadWJet1En, trueHadWJet1Theta, trueHadWJet1Phi);	
      double trueHadWJet2En;
      if ( applyOnshellWmassConstraint_ ) 
      {
        trueHadWJet2En = compHadWJet2En_Wjj(trueHadWJet1P4, trueHadWJet2Theta, trueHadWJet2Phi);
      } 
      else 
      {
        double q2W = x[4];
        trueHadWJet2En = compHadWJet2En_Wjj(trueHadWJet1P4, trueHadWJet2Theta, trueHadWJet2Phi, q2W);
      }
      if ( !(trueHadWJet2En > 0.) ) return 0.;
      trueHadWJet2P4 = buildLorentzVector(trueHadWJet2En, trueHadWJet2Theta, trueHadWJet2Phi);
    }
    LorentzVector trueHadWJet1_HadWJet2_ChargedLeptonP4 = trueHadWJet1P4 + trueHadWJet2P4 + trueChargedLeptonP4;
    if ( trueHadWJet1_HadWJet2_ChargedLeptonP4.mass() >= higgsBosonMass ) return 0.;
    double trueNuEn = compNuStarEn_Hww(trueHadWJet1_HadWJet2_ChargedLeptonP4, trueNuTheta, trueNuPhi);
    if ( !(trueNuEn > 0.) ) return 0.;
    trueNuP4 = buildLorentzVector(trueNuEn, trueNuTheta, trueNuPhi);
  }

  LorentzVector trueSumP4 = trueBJet1P4 + trueBJet2P4 + trueChargedLeptonP4 + trueNuP4 + trueHadWJet1P4 + trueHadWJet2P4;

  // perform boost into zero-transverse-momentum (ZTM) frame 
  double ztmFramePx = trueSumP4.px();
  double ztmFramePy = trueSumP4.py();
  double ztmFrameEn = trueSumP4.energy();
  Vector boost(-ztmFramePx/ztmFrameEn, -ztmFramePy/ztmFrameEn, 0.);
  LorentzVector trueBJet1P4_ztm = ROOT::Math::VectorUtil::boost(trueBJet1P4, boost); 
  LorentzVector trueBJet2P4_ztm = ROOT::Math::VectorUtil::boost(trueBJet2P4, boost); 
  LorentzVector trueChargedLeptonP4_ztm = ROOT::Math::VectorUtil::boost(trueChargedLeptonP4, boost);
  LorentzVector trueNuP4_ztm = ROOT::Math::VectorUtil::boost(trueNuP4, boost);
  LorentzVector trueHadWJet1P4_ztm = ROOT::Math::VectorUtil::boost(trueHadWJet1P4, boost);
  LorentzVector trueHadWJet2P4_ztm = ROOT::Math::VectorUtil::boost(trueHadWJet2P4, boost);
  if ( verbosity_ >= 2 ) 
  {
    std::string chargedLepton_shortLabel;
    std::string chargedLepton_longLabel;
    if ( measuredChargedLepton_->charge() > 0 ) 
    {
      chargedLepton_shortLabel = "lep+";
      chargedLepton_longLabel = "lepton+";
    } 
    else if ( measuredChargedLepton_->charge() < 0 ) 
    {
      chargedLepton_shortLabel = "lep-";
      chargedLepton_longLabel = "lepton-";
    } 
    else assert(0);
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
    printLorentzVector(Form("%s", chargedLepton_longLabel.data()), trueChargedLeptonP4, measuredChargedLepton_->p4());
    printLorentzVector("neutrino", trueNuP4);
    if ( measuredHadWJet1_ ) 
    {
      printLorentzVector("jet1 from W->jj", trueHadWJet1P4, measuredHadWJet1_->p4());
    } 
    else 
    {
      printLorentzVector_NA("jet1 from W->jj", trueHadWJet1P4);
    }
    if ( measuredHadWJet2_ ) 
    {
      printLorentzVector("jet2 from W->jj", trueHadWJet2P4, measuredHadWJet2_->p4());
    } 
    else 
    {
      printLorentzVector_NA("jet2 from W->jj", trueHadWJet2P4);
    }
    std::cout << "m(b bbar) = " << (trueBJet1P4 + trueBJet2P4).mass() << std::endl;
    std::cout << Form("m(%s nu) = ", chargedLepton_shortLabel.data()) << (trueChargedLeptonP4 + trueNuP4).mass() << std::endl;
    std::cout << "m(W->jj) = " << (trueHadWJet1P4 + trueHadWJet2P4).mass() << std::endl;
    std::cout << Form("m(%s nu W->jj) = ", chargedLepton_shortLabel.data()) << (trueChargedLeptonP4 + trueNuP4 + trueHadWJet1P4 + trueHadWJet2P4).mass() << std::endl;
    printLorentzVector("sum", trueSumP4);
    std::cout << "zero-transverse-momentum frame:" << std::endl;
    printLorentzVector("b-jet1", trueBJet1P4_ztm);
    printLorentzVector("b-jet2", trueBJet2P4_ztm);
    printLorentzVector(Form("%s", chargedLepton_longLabel.data()), trueChargedLeptonP4_ztm);
    printLorentzVector("neutrino", trueNuP4_ztm);
    printLorentzVector("jet1 from W->jj", trueHadWJet1P4_ztm);
    printLorentzVector("jet2 from W->jj", trueHadWJet2P4_ztm);
    std::cout << "m(b bbar) = " << (trueBJet1P4_ztm + trueBJet2P4_ztm).mass() << std::endl;
    std::cout << Form("m(%s nu) = ", chargedLepton_shortLabel.data()) << (trueChargedLeptonP4_ztm + trueNuP4_ztm).mass() << std::endl;
    std::cout << "m(W->jj) = " << (trueHadWJet1P4_ztm + trueHadWJet2P4_ztm).mass() << std::endl;
    std::cout << Form("m(%s nu W->jj) = ", chargedLepton_shortLabel.data()) << (trueChargedLeptonP4_ztm + trueNuP4_ztm + trueHadWJet1P4_ztm + trueHadWJet2P4_ztm).mass() << std::endl;
    LorentzVector trueSumP4_ztm = trueBJet1P4_ztm + trueBJet2P4_ztm + trueChargedLeptonP4_ztm + trueNuP4_ztm + trueHadWJet1P4_ztm + trueHadWJet2P4_ztm;
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
  madgraphChargedLeptonP4_[0] = trueChargedLeptonP4_ztm.energy();
  madgraphChargedLeptonP4_[1] = trueChargedLeptonP4_ztm.px();
  madgraphChargedLeptonP4_[2] = trueChargedLeptonP4_ztm.py();
  madgraphChargedLeptonP4_[3] = trueChargedLeptonP4_ztm.pz();
  madgraphNeutrinoP4_[0] = trueNuP4_ztm.energy();
  madgraphNeutrinoP4_[1] = trueNuP4_ztm.px();
  madgraphNeutrinoP4_[2] = trueNuP4_ztm.py();
  madgraphNeutrinoP4_[3] = trueNuP4_ztm.pz();
  madgraphHadWJet1P4_[0] = trueHadWJet1P4_ztm.energy();
  madgraphHadWJet1P4_[1] = trueHadWJet1P4_ztm.px();
  madgraphHadWJet1P4_[2] = trueHadWJet1P4_ztm.py();
  madgraphHadWJet1P4_[3] = trueHadWJet1P4_ztm.pz();
  madgraphHadWJet2P4_[0] = trueHadWJet2P4_ztm.energy();
  madgraphHadWJet2P4_[1] = trueHadWJet2P4_ztm.px();
  madgraphHadWJet2P4_[2] = trueHadWJet2P4_ztm.py();
  madgraphHadWJet2P4_[3] = trueHadWJet2P4_ztm.pz();
  madgraphBJet1P4_[0] = trueBJet1P4_ztm.energy();
  madgraphBJet1P4_[1] = trueBJet1P4_ztm.px();
  madgraphBJet1P4_[2] = trueBJet1P4_ztm.py();
  madgraphBJet1P4_[3] = trueBJet1P4_ztm.pz();
  madgraphBJet2P4_[0] = trueBJet2P4_ztm.energy();
  madgraphBJet2P4_[1] = trueBJet2P4_ztm.px();
  madgraphBJet2P4_[2] = trueBJet2P4_ztm.py();
  madgraphBJet2P4_[3] = trueBJet2P4_ztm.pz();
  double prob_ME = -1.;
  if ( madgraphIsInitialized_ ) 
  {
    if ( measuredChargedLepton_->charge() > 0 ) 
    {
      if ( verbosity_ >= 2 )
      {
        printMadGraphMomenta(madgraphMomenta_chargedLeptonPlus_);
      }
      me_madgraph_chargedLeptonPlus_.setHiggsWidth(2.); // CV: enlarge Higgs boson width to make ME evaluation robust against rounding errors
      me_madgraph_chargedLeptonPlus_.setMomenta(madgraphMomenta_chargedLeptonPlus_);
      me_madgraph_chargedLeptonPlus_.sigmaKin();
      prob_ME = me_madgraph_chargedLeptonPlus_.getMatrixElements()[0];
    } 
    else if ( measuredChargedLepton_->charge() < 0 ) 
    {
      if ( verbosity_ >= 2 )
      {
        printMadGraphMomenta(madgraphMomenta_chargedLeptonMinus_);
      }
      me_madgraph_chargedLeptonMinus_.setHiggsWidth(2.); // CV: enlarge Higgs boson width to make ME evaluation robust against rounding errors
      me_madgraph_chargedLeptonMinus_.setMomenta(madgraphMomenta_chargedLeptonMinus_);
      me_madgraph_chargedLeptonMinus_.sigmaKin();
      prob_ME = me_madgraph_chargedLeptonMinus_.getMatrixElements()[0];
    } 
    else assert(0);
    ++numMatrixElementEvaluations_;
    if ( TMath::IsNaN(prob_ME) ) 
    {
      std::cerr << "Warning: MadGraph returned NaN --> skipping event !!" << std::endl;
      std::string chargedLepton_label;
      if ( measuredChargedLepton_->charge() > 0 ) 
      {
	chargedLepton_label = "lepton+";
      } 
      else if ( measuredChargedLepton_->charge() < 0 ) 
      {
	chargedLepton_label = "lepton-";
      } 
      else assert(0);
      printLorentzVector("b-jet1", trueBJet1P4_ztm);
      printLorentzVector("b-jet2", trueBJet2P4_ztm);
      printLorentzVector(Form("%s", chargedLepton_label.data()), trueChargedLeptonP4_ztm);
      printLorentzVector("neutrino", trueNuP4_ztm);
      printLorentzVector("jet1 from W->jj", trueHadWJet1P4_ztm);
      printLorentzVector("jet2 from W->jj", trueHadWJet2P4_ztm);
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
  double trueHadRecoilPx = -(trueChargedLeptonP4.px() + trueNuP4.px());
  double trueHadRecoilPy = -(trueChargedLeptonP4.py() + trueNuP4.py());
  if ( measuredHadWJet1_ ) 
  {
    trueHadRecoilPx -= trueHadWJet1P4.px();
    trueHadRecoilPy -= trueHadWJet1P4.py();
  }
  if ( measuredHadWJet2_ ) 
  {
    trueHadRecoilPx -= trueHadWJet2P4.px();
    trueHadRecoilPy -= trueHadWJet2P4.py();
  }
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
    if ( !measuredHadWJet1_ || !measuredHadWJet2_ ) 
    {
      std::cout << "Note:";
      if ( !measuredHadWJet1_ ) std::cout << " jet1 from W->jj is 'missing', i.e. not reconstructed, and not included in hadRecoil.";
      if ( !measuredHadWJet2_ ) std::cout << " jet2 from W->jj is 'missing', i.e. not reconstructed, and not included in hadRecoil.";
      std::cout << std::endl;
    }
    if ( !measuredBJet1_ || !measuredBJet2_ ) 
    {
      std::cout << "Note:";
      if ( !measuredBJet1_ ) std::cout << " b-jet1 is 'missing', i.e. not reconstructed, and not included in hadRecoil.";
      if ( !measuredBJet2_ ) std::cout << " b-jet2 is 'missing', i.e. not reconstructed, and not included in hadRecoil.";
      std::cout << std::endl;
    }
  }

  double prob_TF = 1.;
  if ( measuredHadWJet1_ ) 
  {
    prob_TF *= hadWJet1TF_->Eval(trueHadWJet1P4.energy());
  }
  if ( measuredHadWJet2_ ) 
  {
    prob_TF *= hadWJet2TF_->Eval(trueHadWJet2P4.energy());
  }
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
  if ( chargedLeptonPermutation_ == kOnshellChargedLepton ) 
  { 
    //---------------------------------------------------------------------------
    // CV: Jacobi factor is the same for case that W mass constraint is applied to on-shell W boson and for the case that it is not applied
    jacobiFactor *= compJacobiFactor_Wlnu(trueChargedLeptonP4, trueNuP4);
    //---------------------------------------------------------------------------
    LorentzVector trueChargedLepton_Nu_HadWJet1P4 = trueChargedLeptonP4 + trueNuP4 + trueHadWJet1P4;
    jacobiFactor *= compJacobiFactor_Hww(trueChargedLepton_Nu_HadWJet1P4, trueHadWJet2P4);
  } 
  else 
  {
    //---------------------------------------------------------------------------
    // CV: Jacobi factor is the same for case that W mass constraint is applied to on-shell W boson and for the case that it is not applied
    jacobiFactor *= compJacobiFactor_Wjj(trueHadWJet1P4, trueHadWJet2P4);
    //---------------------------------------------------------------------------
    LorentzVector trueHadWJet1_HadWJet2_ChargedLeptonP4 = trueHadWJet1P4 + trueHadWJet2P4 + trueChargedLeptonP4;
    jacobiFactor *= compJacobiFactor_Hww(trueHadWJet1_HadWJet2_ChargedLeptonP4, trueNuP4);
  }  
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "jacobiFactor = " << jacobiFactor << std::endl;
  }

  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  double integrandValue = conversionFactor*normFactor_;
  const double fudgeFactor = 1.e+15;
  integrandValue *= fudgeFactor;
  integrandValue *= trueNuP4.pt();
  integrandValue *= prob_PDF;
  integrandValue *= prob_flux;
  integrandValue *= prob_ME;
  if ( measuredHadWJet1_ ) 
  {
    integrandValue *= trueHadWJet1P4.P();
  }
  else 
  {
    integrandValue *= (1./trueHadWJet1P4.energy());    
  }
  if ( measuredHadWJet2_ ) 
  {
    integrandValue *= trueHadWJet2P4.P();
  }
  else 
  {
    integrandValue *= (1./trueHadWJet2P4.energy());    
  }
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

