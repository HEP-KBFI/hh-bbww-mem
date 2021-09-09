#ifndef hhAnalysis_bbwwMEM_measuredParticleAuxFunctions_h
#define hhAnalysis_bbwwMEM_measuredParticleAuxFunctions_h

#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h" // MeasuredParticle

namespace mem
{

template <typename T>
const T * findGenMatch(const MeasuredParticle * measuredParticle, const std::vector<T *> & genParticles)
{
  const T * bestMatch = nullptr;
  double dRmatch = 1.e+3;
  for ( typename std::vector<T *>::const_iterator genParticle = genParticles.begin();
        genParticle != genParticles.end(); ++genParticle ) {
    double dR = deltaR(measuredParticle->eta(), measuredParticle->phi(), (*genParticle)->eta(), (*genParticle)->phi());
    if ( dR < 0.3 && dR < dRmatch ) 
    {
      bestMatch = *genParticle;
      dRmatch = dR;
    }
  }
  return bestMatch;
}

bool
isHigherPt(const MeasuredParticle * particle1,
           const MeasuredParticle * particle2);

}

#endif
