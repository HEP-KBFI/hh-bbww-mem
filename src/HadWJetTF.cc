#include "hhAnalysis/bbwwMEM/interface/HadWJetTF.h"

#include "tthAnalysis/tthMEM/interface/JetTransferFunction.h" // qJetParams

using namespace mem;

HadWJetTF::HadWJetTF(int verbosity)
  : verbosity_(verbosity)
{}

HadWJetTF::~HadWJetTF()
{}

void HadWJetTF::setInputs(const mem::LorentzVector& measuredP4)
{
  measuredEn_ = measuredP4.energy();
  measuredEta_ = measuredP4.eta();
}

double HadWJetTF::Eval(double trueEn) const
{
  tthMEM::structs::qJetParams q(trueEn, measuredEta_); // trueEta = measuredEta
  double prob = tthMEM::functions::gaussianPDF(measuredEn_, q.mu, q.sigma);
  return prob;
}
