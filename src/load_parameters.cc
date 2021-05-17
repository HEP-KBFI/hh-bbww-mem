#include "tthAnalysis/tthMEM/interface/mg5/read_slha.h"
void read_BSM_couplings(SLHAReader& slha, double & kl, double & kt, double & c2, double & cg, double & c2g)
{
  kl  = slha.get_block_entry("bsm", 188, 1.000000e+00);
  kt  = slha.get_block_entry("bsm", 189, 1.000000e+00);
  c2g = slha.get_block_entry("bsm", 32, 1.000000e+00);
  cg  = slha.get_block_entry("bsm", 31, 1.000000e+00);
  c2  = slha.get_block_entry("bsm", 30, -1.000000e+00);
}
