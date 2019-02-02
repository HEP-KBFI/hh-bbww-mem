#include "hhAnalysis/bbwwMEM/interface/BJetTF.h"

#include "tthAnalysis/tthMEM/interface/JetTransferFunction.h" // bJetParams, gaussianPDF

using namespace mem;

BJetTF::BJetTF(int verbosity)
  : verbosity_(verbosity)
{}

BJetTF::~BJetTF()
{}

void BJetTF::setInputs(const mem::LorentzVector& measuredP4)
{
  measuredEn_ = measuredP4.energy();
  measuredEta_ = measuredP4.eta();
}

double BJetTF::Eval(double trueEn) const
{
  double fb = tthMEM::constants::fb;
  tthMEM::structs::bJetParams b(trueEn, measuredEta_); // trueEta = measuredEta
  double gauss1 =       fb  * tthMEM::functions::gaussianPDF(measuredEn_, b.first.mu, b.first.sigma);
  double gauss2 = (1. - fb) * tthMEM::functions::gaussianPDF(measuredEn_, b.second.mu, b.second.sigma);
  double prob = gauss1 + gauss2;
  return prob;
}
