/**
   \class testMEMbbwwDilepton testMEMbbwwDilepton.cc "hhAnalysis/bbwwMEM/bin/testMEMbbwwDilepton.cc"
   \brief Basic example for the use of the matrix element method (MEM) in the HH->bbWW dilepton channel
*/

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoDilepton.h"
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <FWCore/ParameterSet/interface/FileInPath.h>

namespace
{
  std::string
  findFile(const std::string & fileName)
  {
    const edm::FileInPath inputFile(fileName);
    if(inputFile.fullPath().empty())
    {
      std::cerr << "Error: Cannot find file = " << fileName;
      assert(0);
    }
    return inputFile.fullPath();
  }

  void computeLikelihoodRatio(double prob_S, double probErr_S, double prob_B, double probErr_B, double& likelihoodRatio, double& likelihoodRatioErr)
  {
    double prob_SplusB = prob_S + prob_B;    
    if ( prob_SplusB > 0. ) {
      likelihoodRatio = prob_S/prob_SplusB;
      double prob2_SplusB = mem::square(prob_SplusB);
      likelihoodRatioErr = mem::square((prob_S/prob2_SplusB)*probErr_S) + mem::square((prob_B/prob2_SplusB)*probErr_B);
    } else {
      likelihoodRatio = 0.;
      likelihoodRatioErr = 0.;
    }
  }
}

int
main(int argc __attribute__((unused)), char ** argv __attribute__((unused)))
{

  /* This is a single HH->bbWW->bb lnulnu signal event for testing purposes */
  using namespace mem;

  // define measured momenta of b-jets and charged leptons
  // event taken from /store/mc/RunIIFall17MiniAODv2/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph_correctedcfg/MINIAODSIM
  //                  /PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/EC1B2B42-F5AF-E811-84F1-ECB1D79E5C40.root
  // 1:13:12012
  //
  // at generator level:
  // pt =   7.891 eta = -2.125 phi = -0.064 mass = 0.000 pdgId = +13 status = 1
  // pt = 191.000 eta = -0.941 phi = -0.857 mass = 0.000 pdgId = -14 status = 1
  // pt = 192.500 eta = -0.750 phi = -0.855 mass = 0.000 pdgId = -11 status = 1
  // pt = 157.500 eta = -0.654 phi = -0.869 mass = 0.000 pdgId = +12 status = 1
  // pt = 294.000 eta = +0.011 phi = +2.617 mass = 0.000 pdgId =  -5 status = 23
  // pt =  93.500 eta = +0.016 phi = -2.898 mass = 0.000 pdgId =  +5 status = 23
  const std::vector<MeasuredParticle> measuredParticles_signal = {
    { MeasuredParticle::kElectron, 192.500, -0.750, -0.855,    electronMass, +1 },
    { MeasuredParticle::kMuon,       7.891, -2.125, -0.064,        muonMass, -1 },
    { MeasuredParticle::kBJet,     294.000, +0.011, +2.617, bottomQuarkMass     },
    { MeasuredParticle::kBJet,      93.500, +0.016, -2.898, bottomQuarkMass     },
  };

//std::cout << "trueLep+: Pt = " << measuredParticles_signal[0].pt() << ", eta = " << measuredParticles_signal[0].eta() << ", phi = " << measuredParticles_signal[0].phi() << std::endl;
//MeasuredParticle trueNu(MeasuredParticle::kElectron, 157.500, -0.654, -0.869, 0.000);
//std::cout << "trueNu: En = " << trueNu.energy() << ", Theta = " << trueNu.theta() << std::endl;
//std::cout << "m(lep+ nu) = " << (measuredParticles_signal[0].p4() + trueNu.p4()).mass() << std::endl;
//std::cout << "trueLep-: Pt = " << measuredParticles_signal[1].pt() << ", eta = " << measuredParticles_signal[1].eta() << ", phi = " << measuredParticles_signal[1].phi() << std::endl;
//MeasuredParticle trueAntiNu(MeasuredParticle::kElectron, 191.000, -0.941, -0.857, 0.000);
//std::cout << "trueAntiNu: En = " << trueAntiNu.energy() << ", Theta = " << trueAntiNu.theta() << std::endl;
//std::cout << "m(lep- nu) = " << (measuredParticles_signal[1].p4() + trueAntiNu.p4()).mass() << std::endl;
//std::cout << "m(lep+ nu lep- nu) = " << (measuredParticles_signal[0].p4() + trueNu.p4() + measuredParticles_signal[1].p4() + trueAntiNu.p4()).mass() << std::endl;
  // gen jets:
  // pt = 205.977 eta = -0.011 phi = +2.634 mass = 23.656 partonFlavour = -5
  // pt =  87.232 eta = +0.049 phi = -2.891 mass =  6.625 partonFlavour = +5
  //
  // at reco level:
  //const std::vector<MeasuredParticle> measuredParticles_signal = {
  //  { MeasuredParticle::kElectron, 190.399, -0.750, -0.855,    electronMass, +1 },
  //  { MeasuredParticle::kMuon,       7.945, -2.128, -0.064,        muonMass, -1 },
  //  { MeasuredParticle::kBJet,     185.875, -0.006, +2.630, bottomQuarkMass     },
  //  { MeasuredParticle::kBJet,      94.812, +0.037, -2.917, bottomQuarkMass     },
  //};

  // define measured missing transverse momentum (MET)
  // at generator level:
//std::cout << "trueNuPx = " << 157.500 * std::cos(-0.869) << ", trueAntiNuPx = " << 191.000 * std::cos(-0.857) << std::endl;
  const double measuredMEtPx_signal = 157.500 * std::cos(-0.869) + 191.000 * std::cos(-0.857);
  const double measuredMEtPy_signal = 157.500 * std::sin(-0.869) + 191.000 * std::sin(-0.857);
  // at reco level:
  //const double measuredMEtPt_signal = 214.285;
  //const double measuredMEtPhi_signal = -0.806;
  //const double measuredMEtPx_signal = measuredMEtPt_signal * std::cos(measuredMEtPhi_signal);
  //const double measuredMEtPy_signal = measuredMEtPt_signal * std::sin(measuredMEtPhi_signal);

  // define MET uncertainty matrix
  TMatrixD measuredMEtCov_signal(2,2);
  measuredMEtCov_signal[0][0] = 1364.000;
  measuredMEtCov_signal[1][0] =  -43.125;
  measuredMEtCov_signal[0][1] =  -43.125;
  measuredMEtCov_signal[1][1] = 1006.000;
 
  /* This is a single ttbar->bW bW->blnu blnu background event for testing purposes */

  // define measured momenta of b-jets and charged leptons
  // event taken from /store/mc/RunIIFall17MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM
  //                  /94X_mc2017_realistic_v10-v2/60000/1688A03F-E9EF-E711-9A58-008CFAF356FA.root,
  // 1:3293:3188303
  //
  // at generator level:
  // pt =  51.000 eta = +0.621 phi = -1.824 mass =   0.000 pdgId = +13 status = 1
  // pt =  35.125 eta = +0.889 phi = +2.391 mass =   0.000 pdgId = -14 status = 1
  // pt =  33.625 eta = -1.328 phi = +2.852 mass =   0.000 pdgId = -13 status = 1
  // pt =  53.000 eta = -2.430 phi = -1.762 mass =   0.000 pdgId = +14 status = 1
  // pt = 121.500 eta = +0.807 phi = +1.703 mass =   0.000 pdgId =  -5 status = 23
  // pt =  68.000 eta = -1.410 phi = +0.281 mass =   0.000 pdgId =  +5 status = 23
  // pt =  60.000 eta = -2.484 phi = -2.359 mass =  79.500 pdgId = +24 status = 22
  // pt =  46.125 eta = +1.195 phi = -2.555 mass =  73.500 pdgId = -24 status = 22
  // pt =  32.875 eta = -3.398 phi = -0.803 mass = 171.500 pdgId =  +6 status = 62
  // pt = 109.500 eta = +1.262 phi = +2.094 mass = 177.500 pdgId =  -6 status = 62
  const std::vector<MeasuredParticle> measuredParticles_background = {
    { MeasuredParticle::kMuon,  51.000, +0.621, -1.824,        muonMass, -1 },
    { MeasuredParticle::kMuon,  33.625, -1.328, +2.852,        muonMass, +1 },
    { MeasuredParticle::kBJet, 121.500, +0.807, +1.703, bottomQuarkMass     },
    { MeasuredParticle::kBJet,  68.000, -1.410, +0.281, bottomQuarkMass     },
  };

  // gen jets:
  // pt = 121.405 eta = +0.809 phi = +1.705 mass = 7.441 partonFlavour = -5
  // pt =  67.967 eta = -1.407 phi = +0.277 mass = 6.715 partonFlavour = +5
  //
  // at parton level:
  // pt =  54.490 eta = +0.621 phi = -1.811 mass = 0.106 E = 65.329 pdgId = +13
  // pt =  33.320 eta = +0.963 phi = +2.455 mass = 0.001 E = 49.990 pdgId = -14
  // pt =  33.643 eta = -1.356 phi = +2.856 mass = 0.106 E = 69.639 pdgId = -13
  // pt =  53.635 eta = -2.450 phi = -1.766 mass = 0.002 E = 313.106 pdgId = +14
  // pt = 112.500 eta = +0.895 phi = +1.718 mass = 4.800 E = 160.707 pdgId = -5
  // pt =  67.773 eta = -1.441 phi = +0.278 mass = 4.800 E = 151.271 pdgId = +5
  //
  // at reco level:
  //const std::vector<MeasuredParticle> measuredParticles_background = {
  //  { MeasuredParticle::kMuon,  54.490, +0.621, -1.811,        muonMass, -1 },
  //  { MeasuredParticle::kMuon,  33.643, -1.356, +2.852,        muonMass, +1 },
  //  { MeasuredParticle::kBJet, 140.125, +0.800, +1.708, bottomQuarkMass     },
  //  { MeasuredParticle::kBJet,  52.969, -1.397, +0.261, bottomQuarkMass     },
  //};

  // define measured missing transverse momentum (MET)
  // at generator level:
  const double measuredMEtPx_background = 35.125 * std::cos(+2.391) + 53.000 * std::cos(-1.762);
  const double measuredMEtPy_background = 35.125 * std::sin(+2.391) + 53.000 * std::sin(-1.762);
  // at reco level:
  //const double measuredMEtPt_background = 84.232;
  //const double measuredMEtPhi_background = -2.952;
  //const double measuredMEtPx_background = measuredMEtPt_background * std::cos(measuredMEtPhi_background);
  //const double measuredMEtPy_background = measuredMEtPt_background * std::sin(measuredMEtPhi_background);
  
  // define MET uncertainty matrix
  TMatrixD measuredMEtCov_background(2,2);
  measuredMEtCov_background[0][0] = 688.000;
  measuredMEtCov_background[1][0] = -25.125;
  measuredMEtCov_background[0][1] = -25.125;
  measuredMEtCov_background[1][1] = 846.000;

  // set center-of-mass energy to 13 TeV (LHC run 2)
  const double sqrtS = 13.e+3;
  const std::string pdfName = "MSTW2008lo68cl";

  const std::string madgraphFileName_signal     = "hhAnalysis/bbwwMEM/data/param_hh.dat";
  const std::string madgraphFileName_background = "hhAnalysis/bbwwMEM/data/param_ttbar.dat";

  /// fix (flag=true) mass of charged lepton plus neutrino originating from the decay of the "on-shell" W boson to mW,
  /// or allow the mass to vary during the integration (flag=false)
  ///
  /// Note: flag has an effect on the likelihood of the HH->bbWW signal hypothesis only (not on the likelihood of the ttbar background hypothesis)
  //const bool applyOnshellWmassConstraint_signal = true;
  const bool applyOnshellWmassConstraint_signal = false;

  const int verbosity = 1;
  MEMbbwwAlgoDilepton memAlgo(sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), verbosity);
  memAlgo.applyOnshellWmassConstraint_signal(applyOnshellWmassConstraint_signal);
  memAlgo.setIntMode(MEMbbwwAlgoDilepton::kVAMP);
  memAlgo.setMaxObjFunctionCalls(20000);

  std::cout << "processing signal event:\n";
  std::cout << " m(bb) = " << (measuredParticles_signal[2].p4() + measuredParticles_signal[3].p4()).mass() << std::endl;
  std::cout << " m(ll) = " << (measuredParticles_signal[0].p4() + measuredParticles_signal[1].p4()).mass() << std::endl;
  std::cout << " m(bl+):" 
	    << " 1st permutation = " << (measuredParticles_signal[0].p4() + measuredParticles_signal[2].p4()).mass() << "," 
	    << " 2nd permutation = " << (measuredParticles_signal[0].p4() + measuredParticles_signal[3].p4()).mass() << std::endl;
  std::cout << " m(bl-):" 
	    << " 1st permutation = " << (measuredParticles_signal[1].p4() + measuredParticles_signal[2].p4()).mass() << "," 
	    << " 2nd permutation = " << (measuredParticles_signal[1].p4() + measuredParticles_signal[3].p4()).mass() << std::endl;
  memAlgo.integrate(measuredParticles_signal, measuredMEtPx_signal, measuredMEtPy_signal, measuredMEtCov_signal);
  std::cout << "numMatrixElementEvaluations:" 
	    << " signal = " << memAlgo.getNumMatrixElementEvaluations_signal() << "," 
	    << " background = " << memAlgo.getNumMatrixElementEvaluations_background() << std::endl;
  MEMbbwwAlgoDilepton::resultType result_sig = memAlgo.getResult();
  double likelihoodRatio_sig, likelihoodRatioErr_sig;
  computeLikelihoodRatio(
    result_sig.prob_signal_, result_sig.probErr_signal_, result_sig.prob_background_, result_sig.probErr_background_, 
    likelihoodRatio_sig, likelihoodRatioErr_sig);
  std::cout << " probability for signal hypothesis = "     << result_sig.prob_signal_     << " +/- " << result_sig.probErr_signal_     << "\n"
	    << " probability for background hypothesis = " << result_sig.prob_background_ << " +/- " << result_sig.probErr_background_ << "\n"
	    << "--> likelihood ratio = "                   << likelihoodRatio_sig         << " +/- " << likelihoodRatioErr_sig         << "\n"
  ;
 
  std::cout << "processing background event:\n";
  std::cout << " m(bb) = " << (measuredParticles_background[2].p4() + measuredParticles_background[3].p4()).mass() << std::endl;
  std::cout << " m(ll) = " << (measuredParticles_background[0].p4() + measuredParticles_background[1].p4()).mass() << std::endl;
  std::cout << " m(bl+):" 
	    << " 1st permutation = " << (measuredParticles_background[1].p4() + measuredParticles_background[2].p4()).mass() << "," 
	    << " 2nd permutation = " << (measuredParticles_background[1].p4() + measuredParticles_background[3].p4()).mass() << std::endl;
  std::cout << " m(bl-):" 
	    << " 1st permutation = " << (measuredParticles_background[0].p4() + measuredParticles_background[2].p4()).mass() << "," 
	    << " 2nd permutation = " << (measuredParticles_background[0].p4() + measuredParticles_background[3].p4()).mass() << std::endl;
  memAlgo.integrate(measuredParticles_background, measuredMEtPx_background, measuredMEtPy_background, measuredMEtCov_background);
  std::cout << "numMatrixElementEvaluations:" 
	    << " signal = " << memAlgo.getNumMatrixElementEvaluations_signal() << "," 
	    << " background = " << memAlgo.getNumMatrixElementEvaluations_background() << std::endl;
  MEMbbwwAlgoDilepton::resultType result_bkg = memAlgo.getResult();
  double likelihoodRatio_bkg, likelihoodRatioErr_bkg;
  computeLikelihoodRatio(
    result_bkg.prob_signal_, result_bkg.probErr_signal_, result_bkg.prob_background_, result_bkg.probErr_background_, 
    likelihoodRatio_bkg, likelihoodRatioErr_bkg);
  std::cout << " probability for signal hypothesis = "     << result_bkg.prob_signal_     << " +/- " << result_bkg.probErr_signal_     << "\n"
	    << " probability for background hypothesis = " << result_bkg.prob_background_ << " +/- " << result_bkg.probErr_background_ << "\n"
	    << "--> likelihood ratio = "                   << likelihoodRatio_bkg         << " +/- " << likelihoodRatioErr_bkg         << "\n"
  ;

  return EXIT_SUCCESS;
}
