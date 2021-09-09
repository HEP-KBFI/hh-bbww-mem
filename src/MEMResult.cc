#include "hhAnalysis/bbwwMEM/interface/MEMResult.h"

std::ostream&
operator<<(std::ostream& stream, const MEMResultBase& result)
{
  stream << " probability for signal hypothesis = " << result.getProb_signal() << " +/- " << result.getProbErr_signal() << std::endl;
  stream << " probability for background hypothesis = " << result.getProb_background() << " +/- " << result.getProbErr_background() << std::endl;
  stream << "--> likelihood ratio = " << result.getLikelihoodRatio() << " +/- " << result.getLikelihoodRatioErr() << std::endl;
  return stream;
}

std::ostream&
operator<<(std::ostream& stream, const MEMbbwwResultDilepton& result)
{
  stream << static_cast<const MEMResultBase &>(result);
  return stream;
}

std::ostream&
operator<<(std::ostream& stream, const MEMbbwwResultSingleLepton& result)
{
  stream << static_cast<const MEMResultBase &>(result);
  return stream;
}
