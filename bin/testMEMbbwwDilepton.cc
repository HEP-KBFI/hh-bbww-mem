
/**
   \class testMEMbbwwDilepton testMEMbbwwDilepton.cc "hhAnalysis/bbwwMEM/bin/testMEMbbwwDilepton.cc"
   \brief Basic example for the use of the matrix element method (MEM) in the HH->bbWW dilepton channel
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoDilepton.h"
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"

#include <TMath.h>

namespace
{
  std::string findFile(const std::string& fileName)
  {
    edm::FileInPath inputFile(fileName);
    if ( inputFile.fullPath() == "" ) {
      std::cerr << "Error: Cannot find file = " << fileName << " !!" << std::endl;
      assert(0);
    }
    return inputFile.fullPath().data();
  }
}

using namespace mem;

int main(int argc, char* argv[]) 
{
  /* 
     This is a single HH->bbWW->bb lnulnu signal event for testing purposes.
   */

  // define measured momenta of b-jets and charged leptons
  std::vector<MeasuredParticle> measuredParticles_signal;
  measuredParticles_signal.push_back(MeasuredParticle(MeasuredParticle::kElectron, 25.00, -1.00, +1.00, electronMass, +1)); // positron
  measuredParticles_signal.push_back(MeasuredParticle(MeasuredParticle::kMuon,     15.00, +2.00, -1.00, muonMass,     -1)); // muon
  measuredParticles_signal.push_back(MeasuredParticle(MeasuredParticle::kBJet,     35.00,  0.00, +2.00, bottomQuarkMass));  // first b-jet
  measuredParticles_signal.push_back(MeasuredParticle(MeasuredParticle::kBJet,     25.00, +1.00, -2.00, bottomQuarkMass));  // second b-jet
  
  // define measured missing transverse momentum (MET)
  double measuredMEtPx_signal = +18.24;
  double measuredMEtPy_signal = -23.07;
  
  // define MET uncertainty matrix
  TMatrixD measuredMEtCov_signal(2,2);
  measuredMEtCov_signal[0][0] = 100.00;
  measuredMEtCov_signal[1][0] =   0.00;
  measuredMEtCov_signal[0][1] =   0.00;
  measuredMEtCov_signal[1][1] = 100.00;

  /* 
     This is a single ttbar->bW bW->blnu blnu background event for testing purposes.
   */

  // define measured momenta of b-jets and charged leptons
  std::vector<MeasuredParticle> measuredParticles_background;
  measuredParticles_background.push_back(MeasuredParticle(MeasuredParticle::kElectron, 25.00, -1.00, +1.00, electronMass, +1)); // positron
  measuredParticles_background.push_back(MeasuredParticle(MeasuredParticle::kMuon,     15.00, +2.00, -1.00, muonMass,     -1)); // muon
  measuredParticles_background.push_back(MeasuredParticle(MeasuredParticle::kBJet,     35.00,  0.00, +2.00, bottomQuarkMass));  // first b-jet
  measuredParticles_background.push_back(MeasuredParticle(MeasuredParticle::kBJet,     25.00, +1.00, -2.00, bottomQuarkMass));  // second b-jet
  
  // define measured missing transverse momentum (MET)
  double measuredMEtPx_background = +18.24;
  double measuredMEtPy_background = -23.07;
  
  // define MET uncertainty matrix
  TMatrixD measuredMEtCov_background(2,2);
  measuredMEtCov_background[0][0] = 100.00;
  measuredMEtCov_background[1][0] =   0.00;
  measuredMEtCov_background[0][1] =   0.00;
  measuredMEtCov_background[1][1] = 100.00;

  // CV: set center-of-mass energy to 13 TeV (LHC run 2)
  double sqrtS = 13.e+3; 

  //std::string pdfName = "cteq66";
  std::string pdfName = "MSTW2008lo68cl";

  std::string madgraphFileName_signal = "hhAnalysis/bbwwMEM/data/param_hh.dat";
  std::string madgraphFileName_background = "hhAnalysis/bbwwMEM/data/param_ttbar.dat";

  int verbosity = 1;
  MEMbbwwAlgoDilepton memAlgo(sqrtS, pdfName.data(), findFile(madgraphFileName_signal), findFile(madgraphFileName_background), verbosity);
  memAlgo.setIntMode(MEMbbwwAlgoDilepton::kVAMP);
  memAlgo.setMaxObjFunctionCalls(20000);

  std::cout << "processing signal event:" << std::endl;
  memAlgo.integrate(measuredParticles_signal, measuredMEtPx_signal, measuredMEtPy_signal, measuredMEtCov_signal);
  MEMbbwwAlgoDilepton::resultType result = memAlgo.getResult();
  std::cout << " probability for signal hypothesis = " << result.prob_signal_ << " +/- " << result.probErr_signal_ << std::endl;
  std::cout << " probability for background hypothesis = " << result.prob_background_ << " +/- " << result.probErr_background_ << std::endl;
  double ratio_signal = result.prob_signal_/result.prob_background_;
  double ratioErr_signal = ratio_signal*TMath::Sqrt(square(result.probErr_signal_/result.prob_signal_) + square(result.probErr_background_/result.prob_background_));
  std::cout << "--> likelihood ratio = " << ratio_signal << " +/- " << ratioErr_signal << std::endl;

  std::cout << "processing background event:" << std::endl;
  memAlgo.integrate(measuredParticles_background, measuredMEtPx_background, measuredMEtPy_background, measuredMEtCov_background);
  result = memAlgo.getResult();
  std::cout << " probability for signal hypothesis = " << result.prob_signal_ << " +/- " << result.probErr_signal_ << std::endl;
  std::cout << " probability for background hypothesis = " << result.prob_background_ << " +/- " << result.probErr_background_ << std::endl;
  double ratio_background = result.prob_signal_/result.prob_background_;
  double ratioErr_background = ratio_signal*TMath::Sqrt(square(result.probErr_signal_/result.prob_signal_) + square(result.probErr_background_/result.prob_background_));
  std::cout << "--> likelihood ratio = " << ratio_background << " +/- " << ratioErr_background << std::endl;

  return 0;
}
