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
}

int
main(int argc __attribute__((unused)),
     char ** argv __attribute__((unused)))
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
  // pt = 294.000 eta = +0.011 phi = +2.617 mass = 0.000 pdgId = -5 status = 23
  // pt =  93.500 eta = +0.016 phi = -2.898 mass = 0.000 pdgId = +5 status = 23
  //
  // gen jets:
  // pt = 205.977 eta = -0.011 phi = +2.634 mass = 23.656 partonFlavour = -5
  // pt = 87.232 eta = +0.049 phi = -2.891 mass = 6.625 partonFlavour = +5
  //
  // at reco level:
  const std::vector<MeasuredParticle> measuredParticles_signal = {
    { MeasuredParticle::kElectron, 190.399, -0.750, -0.855, -0.047, +1 },
    { MeasuredParticle::kMuon,       7.945, -2.128, -0.064,  0.106, -1 },
    { MeasuredParticle::kBJet,     185.875, -0.006, +2.630, 21.625     },
    { MeasuredParticle::kBJet,      94.812, +0.037, -2.917, 11.852     },
  };

  // define measured missing transverse momentum (MET)
  const double measuredMEtPt_signal = 214.285;
  const double measuredMEtPhi_signal = -0.806;
  const double measuredMEtPx_signal = measuredMEtPt_signal * std::cos(measuredMEtPhi_signal);
  const double measuredMEtPy_signal = measuredMEtPt_signal * std::sin(measuredMEtPhi_signal);

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
  // pt = 51.000 eta = +0.621 phi = -1.824 mass = 0.000 pdgId = +13 status = 1
  // pt = 35.125 eta = +0.889 phi = +2.391 mass = 0.000 pdgId = -14 status = 1
  // pt = 33.625 eta = -1.328 phi = +2.852 mass = 0.000 pdgId = -13 status = 1
  // pt = 53.000 eta = -2.430 phi = -1.762 mass = 0.000 pdgId = +14 status = 1
  // pt = 121.500 eta = +0.807 phi = +1.703 mass = 0.000 pdgId = -5 status = 23
  // pt = 68.000 eta = -1.410 phi = +0.281 mass = 0.000 pdgId = +5 status = 23
  // pt = 60.000 eta = -2.484 phi = -2.359 mass = 79.500 pdgId = +24 status = 22
  // pt = 46.125 eta = +1.195 phi = -2.555 mass = 73.500 pdgId = -24 status = 22
  // pt = 32.875 eta = -3.398 phi = -0.803 mass = 171.500 pdgId = +6 status = 62
  // pt = 109.500 eta = +1.262 phi = +2.094 mass = 177.500 pdgId = -6 status = 62
  //
  // gen jets:
  // pt = 121.405 eta = +0.809 phi = +1.705 mass = 7.441 partonFlavour = -5
  // pt = 67.967 eta = -1.407 phi = +0.277 mass = 6.715 partonFlavour = +5
  //
  // at parton level:
  // pt = 54.490 eta = +0.621 phi = -1.811 mass = 0.106 E = 65.329 pdgId = +13
  // pt = 33.320 eta = +0.963 phi = +2.455 mass = 0.001 E = 49.990 pdgId = -14
  // pt = 33.643 eta = -1.356 phi = +2.856 mass = 0.106 E = 69.639 pdgId = -13
  // pt = 53.635 eta = -2.450 phi = -1.766 mass = 0.002 E = 313.106 pdgId = +14
  // pt = 112.500 eta = +0.895 phi = +1.718 mass = 4.800 E = 160.707 pdgId = -5
  // pt = 67.773 eta = -1.441 phi = +0.278 mass = 4.800 E = 151.271 pdgId = +5
  //
  // at reco level:
  const std::vector<MeasuredParticle> measuredParticles_background = {
    { MeasuredParticle::kMuon,  51.867, +0.621, -1.825, 0.106, -1 },
    { MeasuredParticle::kMuon,  33.838, -1.327, +2.852, 0.106, +1 },
    { MeasuredParticle::kBJet, 140.125, +0.800, +1.708, 9.930     },
    { MeasuredParticle::kBJet,  52.969, -1.397, +0.261, 8.555     },
  };

  // define measured missing transverse momentum (MET)
  // GenMET pt = 45.500 phi = -2.477
  const double measuredMEtPt_background = 84.232;
  const double measuredMEtPhi_background = -2.952;
  const double measuredMEtPx_background = measuredMEtPt_background * std::cos(measuredMEtPhi_background);
  const double measuredMEtPy_background = measuredMEtPt_background * std::sin(measuredMEtPhi_background);
  
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

  const int verbosity = 0;
  MEMbbwwAlgoDilepton memAlgo(
    sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), verbosity
  );
  memAlgo.setIntMode(MEMbbwwAlgoDilepton::kVAMP);
  memAlgo.setMaxObjFunctionCalls(20000);

  std::cout << "processing signal event:\n";
  memAlgo.integrate(measuredParticles_signal, measuredMEtPx_signal, measuredMEtPy_signal, measuredMEtCov_signal);
  MEMbbwwAlgoDilepton::resultType result_sig = memAlgo.getResult();
  const double ratio_signal    = result_sig.prob_signal_ / result_sig.prob_background_;
  const double ratioErr_signal = ratio_signal * std::sqrt(
    square(result_sig.probErr_signal_ / result_sig.prob_signal_) + square(result_sig.probErr_background_ / result_sig.prob_background_)
  );
  std::cout <<
    " probability for signal hypothesis = "     << result_sig.prob_signal_     << " +/- " << result_sig.probErr_signal_     << "\n"
    " probability for background hypothesis = " << result_sig.prob_background_ << " +/- " << result_sig.probErr_background_ << "\n"
    "--> likelihood ratio = "                   << ratio_signal            << " +/- " << ratioErr_signal                    << '\n'
  ;

  std::cout << "processing background event:\n";
  memAlgo.integrate(measuredParticles_background, measuredMEtPx_background, measuredMEtPy_background, measuredMEtCov_background);
  MEMbbwwAlgoDilepton::resultType result_bkg = memAlgo.getResult();
  const double ratio_background    = result_bkg.prob_signal_ / result_bkg.prob_background_;
  const double ratioErr_background = ratio_background * std::sqrt(
    square(result_bkg.probErr_signal_ / result_bkg.prob_signal_) + square(result_bkg.probErr_background_ / result_bkg.prob_background_)
  );
  std::cout <<
    " probability for signal hypothesis = "     << result_bkg.prob_signal_     << " +/- " << result_bkg.probErr_signal_     << "\n"
    " probability for background hypothesis = " << result_bkg.prob_background_ << " +/- " << result_bkg.probErr_background_ << "\n"
    "--> likelihood ratio = "                   << ratio_background        << " +/- " << ratioErr_background                << '\n'
  ;

  return EXIT_SUCCESS;
}
