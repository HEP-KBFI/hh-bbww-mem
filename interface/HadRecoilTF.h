#ifndef hhAnalysis_bbwwMEM_HadRecoilTF_h
#define hhAnalysis_bbwwMEM_HadRecoilTF_h

#include <TMatrixD.h>

namespace mem
{
  
class HadRecoilTF
{
 public:
  /// error codes that can be read out by SVfitMEM class
  enum ErrorCodes {
    None            = 0x00000000,
    MatrixInversion = 0x00000001
  };

  HadRecoilTF(int = 0);
  ~HadRecoilTF();
  
  /// set measured missing transverse momentum (MET)
  /// and MET uncertainty matrix
  void setInputs(double, double, const TMatrixD&);

  int getErrorCode() const { return errorCode_; }

  /// evaluate transfer function (TF)
  double Eval(double, double) const;

 protected:  
  /// measured momentum of hadronic recoil (rho) and uncertainty matrix 
  /// that quantifies the expected experimental resolution on rho
  double measuredPx_;
  double measuredPy_;
  TMatrixD cov_;
  TMatrixD covInverse_;
  double covInverse_xx_;
  double covInverse_xy_;
  double covInverse_yx_; 
  double covInverse_yy_;
  double normFactor_;

  /// error code to indicate that inversion of the uncertainty matrix failed
  int errorCode_;

  /// verbosity level
  int verbosity_;
};

}

#endif
