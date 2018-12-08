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
  const std::vector<MeasuredParticle> measuredParticles_signal = {
    { MeasuredParticle::kElectron, 25.00, -1.00, +1.00, electronMass,   +1 }, // positron
    { MeasuredParticle::kMuon,     15.00, +2.00, -1.00, muonMass,       -1 }, // muon
    { MeasuredParticle::kBJet,     35.00,  0.00, +2.00, bottomQuarkMass    }, // first b-jet
    { MeasuredParticle::kBJet,     25.00, +1.00, -2.00, bottomQuarkMass    }, // second b-jet
  };
  
  // define measured missing transverse momentum (MET)
  const double measuredMEtPx_signal = +18.24;
  const double measuredMEtPy_signal = -23.07;

  // define MET uncertainty matrix
  TMatrixD measuredMEtCov_signal(2,2);
  measuredMEtCov_signal[0][0] = 100.00;
  measuredMEtCov_signal[1][0] =   0.00;
  measuredMEtCov_signal[0][1] =   0.00;
  measuredMEtCov_signal[1][1] = 100.00;

  /* This is a single ttbar->bW bW->blnu blnu background event for testing purposes */

  // define measured momenta of b-jets and charged leptons
  const std::vector<MeasuredParticle> measuredParticles_background = {
    { MeasuredParticle::kElectron, 25.00, -1.00, +1.00, electronMass,   +1 }, // positron
    { MeasuredParticle::kMuon,     15.00, +2.00, -1.00, muonMass,       -1 }, // muon
    { MeasuredParticle::kBJet,     35.00,  0.00, +2.00, bottomQuarkMass    }, // first b-jet
    { MeasuredParticle::kBJet,     25.00, +1.00, -2.00, bottomQuarkMass    }, // second b-jet
  };
  
  // define measured missing transverse momentum (MET)
  const double measuredMEtPx_background = +18.24;
  const double measuredMEtPy_background = -23.07;
  
  // define MET uncertainty matrix
  TMatrixD measuredMEtCov_background(2,2);
  measuredMEtCov_background[0][0] = 100.00;
  measuredMEtCov_background[1][0] =   0.00;
  measuredMEtCov_background[0][1] =   0.00;
  measuredMEtCov_background[1][1] = 100.00;

  // set center-of-mass energy to 13 TeV (LHC run 2)
  const double sqrtS = 13.e+3;
  const std::string pdfName = "MSTW2008lo68cl";

  const std::string madgraphFileName_signal     = "hhAnalysis/bbwwMEM/data/param_hh.dat";
  const std::string madgraphFileName_background = "hhAnalysis/bbwwMEM/data/param_ttbar.dat";

  const int verbosity = 2;
  MEMbbwwAlgoDilepton memAlgo(
    sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), verbosity
  );
  memAlgo.setIntMode(MEMbbwwAlgoDilepton::kVAMP);
  memAlgo.setMaxObjFunctionCalls(20000);

  std::cout << "processing signal event:\n";
  memAlgo.integrate(measuredParticles_signal, measuredMEtPx_signal, measuredMEtPy_signal, measuredMEtCov_signal);
  MEMbbwwAlgoDilepton::resultType result = memAlgo.getResult();
  const double ratio_signal    = result.prob_signal_ / result.prob_background_;
  const double ratioErr_signal = ratio_signal * std::sqrt(
    square(result.probErr_signal_ / result.prob_signal_) + square(result.probErr_background_ / result.prob_background_)
  );
  std::cout <<
    " probability for signal hypothesis = "     << result.prob_signal_     << " +/- " << result.probErr_signal_     << "\n"
    " probability for background hypothesis = " << result.prob_background_ << " +/- " << result.probErr_background_ << "\n"
    "--> likelihood ratio = "                   << ratio_signal            << " +/- " << ratioErr_signal            << '\n'
  ;

  std::cout << "processing background event:\n";
  memAlgo.integrate(measuredParticles_background, measuredMEtPx_background, measuredMEtPy_background, measuredMEtCov_background);
  result = memAlgo.getResult();
  const double ratio_background    = result.prob_signal_ / result.prob_background_;
  const double ratioErr_background = ratio_signal * std::sqrt(
    square(result.probErr_signal_ / result.prob_signal_) + square(result.probErr_background_ / result.prob_background_)
  );
  std::cout <<
    " probability for signal hypothesis = "     << result.prob_signal_     << " +/- " << result.probErr_signal_     << "\n"
    " probability for background hypothesis = " << result.prob_background_ << " +/- " << result.probErr_background_ << "\n"
    "--> likelihood ratio = "                   << ratio_background        << " +/- " << ratioErr_background        << '\n'
  ;

  return EXIT_SUCCESS;
}
