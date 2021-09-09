#include "hhAnalysis/bbwwMEM/interface/measuredParticleAuxFunctions.h"

bool
mem::isHigherPt(const MeasuredParticle * particle1,
                const MeasuredParticle * particle2)
{
  return particle1->pt() > particle2->pt();
}
