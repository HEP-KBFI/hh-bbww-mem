#ifndef hhAnalysis_bbwwMEM_BJetTF_h
#define hhAnalysis_bbwwMEM_BJetTF_h

#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h" // mem::LorentzVector

namespace mem
{
  
class BJetTF
{
 public:
  BJetTF(int = 0);
  ~BJetTF();
  
  /// set measured b-jet energy and pseudo-rapidity
  void setInputs(const mem::LorentzVector&);

  /// evaluate transfer function (TF)
  double Eval(double) const;

 protected:  
  /// measured b-jet energy and pseudo-rapidity
  double measuredEn_;
  double measuredEta_;

  /// verbosity level
  int verbosity_;
};

}

#endif
