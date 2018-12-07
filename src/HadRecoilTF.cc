#include "hhAnalysis/bbwwMEM/interface/HadRecoilTF.h"

#include <TMath.h>

#include <iostream>

using namespace mem;

HadRecoilTF::HadRecoilTF(int verbosity)
  : errorCode_(0)
  , verbosity_(verbosity)
{}

HadRecoilTF::~HadRecoilTF()
{}

void HadRecoilTF::setInputs(double measuredPx, double measuredPy, const TMatrixD& cov)
{
  measuredPx_ = measuredPx;
  measuredPy_ = measuredPy;
  cov_.ResizeTo(2,2);
  cov_ = cov;

  // reset 'MatrixInversion' error code
  errorCode_ &= (errorCode_ ^ MatrixInversion);

  // invert uncertainty matrix
  covInverse_.ResizeTo(2,2);
  covInverse_ = cov_;
  double covDet = covInverse_.Determinant();
  normFactor_ = 0.;
  if ( covDet != 0. ) { 
    covInverse_.Invert(); 
    covInverse_xx_ = covInverse_(0,0);
    covInverse_xy_ = covInverse_(0,1);
    covInverse_yx_ = covInverse_(1,0);
    covInverse_yy_ = covInverse_(1,1);
    normFactor_ = 1./(2.*TMath::Pi()*TMath::Sqrt(covDet));
  } else {
    std::cerr << "Error: Failed to invert uncertainty matrix for hadronic recoil (det=0) !!" << std::endl;
    errorCode_ |= MatrixInversion;
  }
}

double HadRecoilTF::Eval(double truePx, double truePy) const
{
  if ( errorCode_ ) {
    return 0.;
  }
  double residualPx = measuredPx_ - truePx;
  double residualPy = measuredPy_ - truePy;
  double pull2 = residualPx*(covInverse_xx_*residualPx + covInverse_xy_*residualPy) + residualPy*(covInverse_yx_*residualPx + covInverse_yy_*residualPy);
  double prob = normFactor_*TMath::Exp(-0.5*pull2);
  return prob;
}
 
