#ifndef hhAnalysis_bbwwMEM_HadWJetTF_h
#define hhAnalysis_bbwwMEM_HadWJetTF_h

#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h" // mem::LorentzVector

namespace mem
{
  
class HadWJetTF
{
 public:
  HadWJetTF(int = 0);
  ~HadWJetTF();
  
  /// set measured missing transverse momentum (MET)
  /// and MET uncertainty matrix
  void setInputs(const mem::LorentzVector&);

  /// evaluate transfer function (TF)
  double Eval(double) const;

 protected:  
  /// measured (light quark) jet energy and pseudo-rapidity
  double measuredEn_;
  double measuredEta_;

  /// verbosity level
  int verbosity_;
};

}

#endif
