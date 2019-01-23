#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandSingleLepton_background.h"

#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TMath.h>
#include <Math/VectorUtil.h> // ROOT::Math::VectorUtil::boost
#include <TString.h> // Form

using namespace mem;

MEMbbwwIntegrandSingleLepton_background::MEMbbwwIntegrandSingleLepton_background(double sqrtS, const std::string& madgraphFileName, int verbosity)
  : MEMbbwwIntegrandSingleLepton(sqrtS, madgraphFileName, verbosity)
{
  if ( verbosity_ ) 
  {
    std::cout << "<MEMbbwwIntegrandSingleLepton_background::MEMbbwwIntegrandSingleLepton_background>:" << std::endl;
  }

  /// define integration variables
  intNumDimensions_ = 3;  
  intBounds_lower_ = new double[intNumDimensions_];
  intBounds_upper_ = new double[intNumDimensions_];
  intVarNames_.clear();
  intVarNames_.push_back("NuTheta"); 
  intBounds_lower_[0] = 0.;
  intBounds_upper_[0] = TMath::Pi();
  intVarNames_.push_back("NuPhi"); 
  intBounds_lower_[1] = -TMath::Pi();
  intBounds_upper_[1] = +TMath::Pi();
  intVarNames_.push_back("HadWJet1En"); 
  intBounds_lower_[2] = -1.; // to be set as function of measured jet energy and expected jet energy resolution
  intBounds_upper_[2] = -1.;

  // initialize MadGraph
  if ( madgraphFileName != "" ) 
  {
    std::cout << "initializing MadGraph ME for ttbar background using " << madgraphFileName << " file." << std::endl;
    me_madgraph_chargedLeptonPlus_.initProc(madgraphFileName);
    me_madgraph_chargedLeptonMinus_.initProc(madgraphFileName);
    madgraphIsInitialized_ = true;
  } 
  else 
  {
    std::cerr << "Error in <MEMbbwwIntegrandSingleLepton_background>: No param.dat file for MadGraph given !!" << std::endl;
    assert(0);
  }

  // define ordering of four-vectors passed to MadGraph when evaluating the matrix element
  // for events containing leptons of positive charge:
  //   1) first incoming gluon
  //   2) second incoming gluon
  //   3) outgoing charged lepton
  //   4) outgoing neutrino
  //   5) outgoing b quark (b-jet)
  //   6) first outgoing jet from W->jj decay
  //   7) second outgoing jet from W->jj decay
  //   8) outgoing anti-b quark (b-jet)
  //
  // Note: the ordering needs to match the labels 1-8 as they are defined in the Feynman diagram doc/mg5/ttbar_sl_pos.ps
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphGluon1P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphGluon2P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphChargedLeptonP4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphNeutrinoP4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphBJet1P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphHadWJet1P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphHadWJet2P4_);
  madgraphMomenta_chargedLeptonPlus_.push_back(madgraphBJet2P4_);

  // define ordering of four-vectors passed to MadGraph when evaluating the matrix element
  // for events containing leptons of negative charge:
  //   1) first incoming gluon
  //   2) second incoming gluon
  //   3) first outgoing jet from W->jj decay
  //   4) second outgoing jet from W->jj decay
  //   5) outgoing b quark (b-jet)
  //   6) outgoing charged lepton
  //   7) outgoing neutrino
  //   8) outgoing anti-b quark (b-jet)
  //
  // Note: the ordering needs to match the labels 1-8 as they are defined in the Feynman diagram doc/mg5/ttbar_sl_neg.ps
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphGluon1P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphGluon2P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphHadWJet1P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphHadWJet2P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphBJet1P4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphChargedLeptonP4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphNeutrinoP4_);
  madgraphMomenta_chargedLeptonMinus_.push_back(madgraphBJet2P4_);
}

MEMbbwwIntegrandSingleLepton_background::~MEMbbwwIntegrandSingleLepton_background()
{}

void 
MEMbbwwIntegrandSingleLepton_background::setInputs(const MeasuredParticle* measuredChargedLepton, 
						   const MeasuredParticle* measuredHadWJet1, const MeasuredParticle* measuredHadWJet2,
						   const MeasuredParticle* measuredBJet1, const MeasuredParticle* measuredBJet2,
						   double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
{
  if ( verbosity_ ) 
  {
    std::cout << "<MEMbbwwIntegrandSingleLepton_background::setInputs>:" << std::endl;
  }

  MEMbbwwIntegrandSingleLepton::setInputs(
    measuredChargedLepton, 
    measuredHadWJet1, measuredHadWJet2,
    measuredBJet1, measuredBJet2,
    measuredMEtPx, measuredMEtPy, measuredMEtCov);

  // set integration boundary for energy of first jet from W->jj decay
  // to measured jet energy +/- 3 times expected energy resolution (rough estimate),
  // in order to increase efficiency of numeric integration
  double measuredHadWJet1En = measuredHadWJet1_->energy();
  double measuredHadWJet1EnRes = 0.50*TMath::Sqrt(measuredHadWJet1_->energy());
  intBounds_lower_[2] = TMath::Max(10., measuredHadWJet1En - 3.*measuredHadWJet1EnRes);
  intBounds_upper_[2] = measuredHadWJet1En + 3.*measuredHadWJet1EnRes;

  // Cross section for Standard Model (SM) ttbar production @ 13 TeV center-of-mass energy
  // time branching fraction for the decay ttbar->bbWW->bblnulnu (excluding electrons and muons from tau decays)
  //
  // Note: SM cross section is taken from next-to-next-to-leading order (NNLO) computation, 
  //       while MadGraph matrix element is leading order (LO). 
  //       We expect this inconsistency to have little practical effect.
  const double crossSection_background = 831.76*2.*0.216*0.6742; // [pb]

  // compute product of constant terms in the integrand (terms that do not depend on the integration variables),
  // in order to reduce computing time required for numeric integration
  double numerator = square(topQuarkMass*topQuarkWidth)*square(wBosonMass*wBosonWidth);
  double denominator = 1.;
  denominator *= (TMath::Power(2., 26)*TMath::Power(TMath::Pi(), 17));
  denominator *= measuredChargedLepton_->energy();
  denominator *= (measuredHadWJet1_->p()*measuredHadWJet1_->energy());
  denominator *= (measuredHadWJet2_->p()*measuredHadWJet2_->energy());
  denominator *= (measuredBJet1_->p()*measuredBJet1_->energy());
  denominator *= (measuredBJet2_->p()*measuredBJet2_->energy());
  denominator *= (crossSection_background*square(sqrtS_));
  assert(denominator > 0.);
  normFactor_ = numerator/denominator;
}

double MEMbbwwIntegrandSingleLepton_background::Eval(const double* x) const
{
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "<MEMbbwwIntegrandSingleLepton_background::Eval>:" << std::endl;
  }

  const LorentzVector& trueChargedLeptonP4 = measuredChargedLepton_->p4();

  double trueNuTheta = x[0];
  double trueNuPhi = x[1];
  double trueNuEn = compNuEn_Wlnu(trueChargedLeptonP4, trueNuTheta, trueNuPhi);
  if ( !(trueNuEn > 0.) ) return 0.;
  LorentzVector trueNuP4 = buildLorentzVector(trueNuEn, trueNuTheta, trueNuPhi);

  std::vector<double> trueBJet1En_solutions = compBJetEn_top(trueChargedLeptonP4 + trueNuP4, measuredBJet1_->theta(), measuredBJet1_->phi());
  bool trueBJet1En_foundSolution = false;
  LorentzVector trueBJet1P4 = findBJetEn_solution_top(trueBJet1En_solutions, measuredBJet1_->theta(), measuredBJet1_->phi(), trueChargedLeptonP4 + trueNuP4, trueBJet1En_foundSolution);
  if ( !trueBJet1En_foundSolution ) return 0.;

  double trueHadWJet1En = x[2];
  LorentzVector trueHadWJet1P4 = buildLorentzVector(trueHadWJet1En, measuredHadWJet1_->theta(), measuredHadWJet1_->phi());

  double trueHadWJet2En = compHadWJet2En_Wjj(trueHadWJet1P4, measuredHadWJet2_->theta(), measuredHadWJet2_->phi());
  if ( !(trueHadWJet2En > 0.) ) return 0.;
  LorentzVector trueHadWJet2P4 = buildLorentzVector(trueHadWJet2En, measuredHadWJet2_->theta(), measuredHadWJet2_->phi());

  std::vector<double> trueBJet2En_solutions = compBJetEn_top(trueHadWJet1P4 + trueHadWJet2P4, measuredBJet2_->theta(), measuredBJet2_->phi());
  bool trueBJet2En_foundSolution = false;
  LorentzVector trueBJet2P4 = findBJetEn_solution_top(trueBJet2En_solutions, measuredBJet2_->theta(), measuredBJet2_->phi(), trueHadWJet1P4 + trueHadWJet2P4, trueBJet2En_foundSolution);
  if ( !trueBJet2En_foundSolution ) return 0.;

  LorentzVector trueSumP4 = trueBJet1P4 + trueChargedLeptonP4 + trueNuP4 + trueBJet2P4 + trueHadWJet1P4 + trueHadWJet2P4;

  // perform boost into zero-transverse-momentum (ZTM) frame 
  double ztmFramePx = trueSumP4.px();
  double ztmFramePy = trueSumP4.py();
  double ztmFrameEn = trueSumP4.energy();
  Vector boost(-ztmFramePx/ztmFrameEn, -ztmFramePy/ztmFrameEn, 0.);
  LorentzVector trueChargedLeptonP4_ztm = ROOT::Math::VectorUtil::boost(trueChargedLeptonP4, boost);
  LorentzVector trueNuP4_ztm = ROOT::Math::VectorUtil::boost(trueNuP4, boost);
  LorentzVector trueBJet1P4_ztm = ROOT::Math::VectorUtil::boost(trueBJet1P4, boost); 
  LorentzVector trueHadWJet1P4_ztm = ROOT::Math::VectorUtil::boost(trueHadWJet1P4, boost);
  LorentzVector trueHadWJet2P4_ztm = ROOT::Math::VectorUtil::boost(trueHadWJet2P4, boost);
  LorentzVector trueBJet2P4_ztm = ROOT::Math::VectorUtil::boost(trueBJet2P4, boost); 
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "laboratory frame:" << std::endl;
    std::string chargedLepton_shortLabel;
    std::string chargedLepton_longLabel;
    std::string bjet1_label;
    std::string bjet2_label;
    if ( measuredChargedLepton_->charge() > 0 ) 
    {
      chargedLepton_shortLabel = "lep+";
      chargedLepton_longLabel = "lepton+";
      bjet1_label = "b";
      bjet2_label = "bbar";
    } 
    else if ( measuredChargedLepton_->charge() < 0 ) 
    {
      chargedLepton_shortLabel = "lep-";
      chargedLepton_longLabel = "lepton-";
      bjet1_label = "bbar";
      bjet2_label = "b";
    } else assert(0);
    printLorentzVector(Form("%s", chargedLepton_longLabel.data()), trueChargedLeptonP4, measuredChargedLepton_->p4());
    printLorentzVector("neutrino", trueNuP4);
    printLorentzVector(Form("%s", bjet1_label.data()), trueBJet1P4, measuredBJet1_->p4());
    printLorentzVector("jet1 from W->jj", trueHadWJet1P4, measuredHadWJet1_->p4());
    printLorentzVector("jet2 from W->jj", trueHadWJet2P4, measuredHadWJet2_->p4());
    printLorentzVector(Form("%s", bjet2_label.data()), trueBJet2P4, measuredBJet2_->p4());
    std::cout << Form("m(%s nu) = ", chargedLepton_shortLabel.data()) << (trueChargedLeptonP4 + trueNuP4).mass() << std::endl;
    std::cout << "m(W->jj) = " << (trueHadWJet1P4 + trueHadWJet2P4).mass() << std::endl;
    std::cout << Form("m(%s %s nu) = ", bjet1_label.data(), chargedLepton_shortLabel.data()) << (trueBJet1P4 + trueChargedLeptonP4 + trueNuP4).mass() << std::endl;
    std::cout << Form("m(%s W->jj) = ", bjet2_label.data()) << (trueBJet2P4 + trueHadWJet1P4 + trueHadWJet2P4).mass() << std::endl;
    printLorentzVector("sum", trueSumP4);
    std::cout << "zero-transverse-momentum frame:" << std::endl;
    printLorentzVector(Form("%s", chargedLepton_longLabel.data()), trueChargedLeptonP4_ztm);
    printLorentzVector("neutrino", trueNuP4_ztm);
    printLorentzVector(Form("%s", bjet1_label.data()), trueBJet1P4_ztm);
    printLorentzVector("jet1 from W->jj", trueHadWJet1P4_ztm);
    printLorentzVector("jet2 from W->jj", trueHadWJet2P4_ztm);
    printLorentzVector(Form("%s", bjet2_label.data()), trueBJet2P4_ztm);
    std::cout << Form("m(%s nu) = ", chargedLepton_shortLabel.data()) << (trueChargedLeptonP4_ztm + trueNuP4_ztm).mass() << std::endl;
    std::cout << "m(W->jj) = " << (trueHadWJet1P4_ztm + trueHadWJet2P4_ztm).mass() << std::endl;    
    std::cout << Form("m(%s %s nu) = ", bjet1_label.data(), chargedLepton_shortLabel.data()) << (trueBJet1P4_ztm + trueChargedLeptonP4_ztm + trueNuP4_ztm).mass() << std::endl;
    std::cout << Form("m(%s W->jj) = ", bjet2_label.data()) << (trueBJet2P4_ztm + trueHadWJet1P4_ztm + trueHadWJet2P4_ztm).mass() << std::endl;
    LorentzVector trueSumP4_ztm = trueBJet1P4_ztm + trueChargedLeptonP4_ztm + trueNuP4_ztm + trueBJet2P4_ztm + trueHadWJet1P4_ztm + trueHadWJet2P4_ztm;
    printLorentzVector("sum", trueSumP4_ztm);
  }

  // compute Bjorken-x of incoming protons and evaluate PDF factor
  double trueSumPz = trueSumP4.pz();
  double trueSumEn = trueSumP4.E();
  double xa = (trueSumEn + trueSumPz)/sqrtS_;
  double xb = (trueSumEn - trueSumPz)/sqrtS_;
  if ( xa <= 0. || xa >= 1. ) return 0.;
  if ( xb <= 0. || xb >= 1. ) return 0.;
  double Q = 2.*topQuarkMass;
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
  madgraphBJet1P4_[0] = trueBJet1P4_ztm.energy();
  madgraphBJet1P4_[1] = trueBJet1P4_ztm.px();
  madgraphBJet1P4_[2] = trueBJet1P4_ztm.py();
  madgraphBJet1P4_[3] = trueBJet1P4_ztm.pz();
  madgraphHadWJet1P4_[0] = trueHadWJet1P4_ztm.energy();
  madgraphHadWJet1P4_[1] = trueHadWJet1P4_ztm.px();
  madgraphHadWJet1P4_[2] = trueHadWJet1P4_ztm.py();
  madgraphHadWJet1P4_[3] = trueHadWJet1P4_ztm.pz();
  madgraphHadWJet2P4_[0] = trueHadWJet2P4_ztm.energy();
  madgraphHadWJet2P4_[1] = trueHadWJet2P4_ztm.px();
  madgraphHadWJet2P4_[2] = trueHadWJet2P4_ztm.py();
  madgraphHadWJet2P4_[3] = trueHadWJet2P4_ztm.pz();
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
      me_madgraph_chargedLeptonMinus_.setMomenta(madgraphMomenta_chargedLeptonMinus_);
      me_madgraph_chargedLeptonMinus_.sigmaKin();
      prob_ME = me_madgraph_chargedLeptonMinus_.getMatrixElements()[0];
    } else assert(0);
    ++numMatrixElementEvaluations_;
    if ( TMath::IsNaN(prob_ME) ) 
    {
      std::cerr << "Warning: MadGraph returned NaN --> skipping event !!" << std::endl;
      std::string chargedLepton_label;
      std::string bjet1_label;
      std::string bjet2_label;
      if ( measuredChargedLepton_->charge() > 0 ) 
      {
	chargedLepton_label = "lepton+";
	bjet1_label = "b";
	bjet2_label = "bbar";
      } 
      else if ( measuredChargedLepton_->charge() < 0 ) 
      {
	chargedLepton_label = "lepton-";
	bjet1_label = "bbar";
	bjet2_label = "b";
      } 
      else assert(0);
      printLorentzVector(Form("%s", chargedLepton_label.data()), trueChargedLeptonP4_ztm);
      printLorentzVector("neutrino", trueNuP4_ztm);
      printLorentzVector(Form("%s", bjet1_label.data()), trueBJet1P4_ztm);
      printLorentzVector("jet1 from W->jj", trueHadWJet1P4_ztm);
      printLorentzVector("jet2 from W->jj", trueHadWJet2P4_ztm);
      printLorentzVector(Form("%s", bjet2_label.data()), trueBJet2P4_ztm);
      return 0.;
    }
  }
  assert(prob_ME >= 0.);
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "prob_ME = " << prob_ME << std::endl;
  }

  double trueHadRecoilPx = -(trueChargedLeptonP4.px() + trueNuP4.px() + trueBJet1P4.px() + trueHadWJet1P4.px() + trueHadWJet2P4.px() + trueBJet2P4.px());
  double trueHadRecoilPy = -(trueChargedLeptonP4.py() + trueNuP4.py() + trueBJet1P4.py() + trueHadWJet1P4.py() + trueHadWJet2P4.py() + trueBJet2P4.py());
  if ( verbosity_ >= 2 )
  {
    std::cout << "hadRecoil:" << std::endl;
    std::cout << " true Px = " << trueHadRecoilPx << ", Py = " << trueHadRecoilPy << std::endl;
    std::cout << " rec. Px = " << measuredHadRecoilPx_ << ", Py = " << measuredHadRecoilPy_ << std::endl;
  }

  double prob_TF = bjet1TF_->Eval(trueBJet1P4.energy())*bjet2TF_->Eval(trueBJet2P4.energy());
  prob_TF *= (hadWJet1TF_->Eval(trueHadWJet1P4.energy())*hadWJet2TF_->Eval(trueHadWJet2P4.energy()));
  prob_TF *= hadRecoilTF_->Eval(trueHadRecoilPx, trueHadRecoilPy);
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "prob_TF = " << prob_TF << std::endl;
  }

  double jacobiFactor = 1.;
  jacobiFactor *= (compJacobiFactor_Wlnu(trueChargedLeptonP4, trueNuP4)*compJacobiFactor_top(trueChargedLeptonP4 + trueNuP4, trueBJet1P4));
  jacobiFactor *= (compJacobiFactor_Wjj(trueHadWJet1P4, trueHadWJet2P4)*compJacobiFactor_top(trueHadWJet1P4 + trueHadWJet2P4, trueBJet2P4));
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "jacobiFactor = " << jacobiFactor << std::endl;
  }

  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  double integrandValue = conversionFactor*normFactor_;
  integrandValue *= trueNuP4.pt();
  integrandValue *= prob_PDF;
  integrandValue *= prob_flux;
  integrandValue *= prob_ME;
  integrandValue *= (trueBJet1P4.P()*trueBJet2P4.P()*trueHadWJet1P4.P()*trueHadWJet2P4.P()*prob_TF);
  integrandValue *= jacobiFactor;
  if ( verbosity_ >= 2 ) 
  {
    std::cout << "--> integrandValue = " << integrandValue << std::endl;
  }

  return integrandValue;
}

