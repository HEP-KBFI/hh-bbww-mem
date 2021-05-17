#ifndef hhAnalysis_bbwwMEM_MEMbbwwAlgoDilepton_h
#define hhAnalysis_bbwwMEM_MEMbbwwAlgoDilepton_h

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoBase.h"
#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandDilepton_signal.h"
#include "hhAnalysis/bbwwMEM/interface/MEMbbwwIntegrandDilepton_background.h"
#include "hhAnalysis/bbwwMEM/interface/MEMIntegratorBase.h"
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"
#include "hhAnalysis/bbwwMEM/interface/MEMResult.h"

#include <TBenchmark.h>
#include <TMatrixD.h>
#include <TMath.h>

#include <vector> 

class MEMbbwwAlgoDilepton : public MEMbbwwAlgoBase
{
 public:
  MEMbbwwAlgoDilepton(double, const std::string&, const std::string&, const std::string&, int = 0);
  ~MEMbbwwAlgoDilepton();

  /// set transfer functions for b-jets and MET
  void setBJet1TF(mem::BJetTF*);
  void setBJet2TF(mem::BJetTF*);
  void setHadRecoilTF(mem::HadRecoilTF*);

  /// fix (flag=true) mass of charged lepton plus neutrino originating from the decay of the "on-shell" W boson to mW,
  /// or allow the mass to vary during the integration (flag=false)
  ///
  /// Note: flag has an effect on the likelihood of the HH->bbWW signal hypothesis only (not on the likelihood of the ttbar background hypothesis)
  void applyOnshellWmassConstraint_signal(bool);
  
  /// run integration 
  void integrate(const std::vector<mem::MeasuredParticle>&, double, double, const TMatrixD&);

  /// return probabilities for signal and background hypotheses
  MEMbbwwResultDilepton getResult() const { return result_; }

 protected:
  /// set measured momenta of charged leptons, b-jets, and (light quark) jets originating from hadronic decay of W boson
  void setMeasuredParticles(const std::vector<mem::MeasuredParticle>&);

  /// pointers to integration classes for signal and background hypotheses
  mem::MEMbbwwIntegrandDilepton_signal* integrand_signal_;
  bool integrand_signal_applyOnshellWmassConstraint_;	
  mem::MEMbbwwIntegrandDilepton_background* integrand_background_;

  /// measured momenta of charged leptons and b-jets
  std::vector<mem::MeasuredParticle> measuredParticles_;
  const mem::MeasuredParticle* measuredChargedLeptonPlus_;
  const mem::MeasuredParticle* measuredChargedLeptonMinus_;
  const mem::MeasuredParticle* measuredLeadingBJet_; 
  const mem::MeasuredParticle* measuredSubleadingBJet_;

  /// result of integration (probabilities for signal and background hypotheses)
  MEMbbwwResultDilepton result_;
};

#endif
