# hh-bbww-mem
Matrix element method (MEM) for HH->bbWW analysis

## Generating the MG files

1. download & unpack Madgraph version 2.3.3: https://launchpad.net/mg5amcnlo/2.0/2.3.x/+download/MG5_aMC_v2.3.3.tar.gz
2. download UFO files: https://cms-project-generators.web.cern.ch/cms-project-generators/BSM_gg_hh.tar
3. unpack the UFO files to `MG5_aMC_v2_3_3/models/BSM_gg_hh`
4. in `./MG5_aMC_v2_3_3/bin/mg5_aMC` run:
```
generate g g > t t~, ( t > w+ b, w+ > l+ vl ), ( t~ > w- b~, w- > l- vl~ )
output standalone ttbar
output standalone_cpp ttbar
display diagrams

generate g g > t t~, ( t > w+ b, w+ > l+ vl ), ( t~ > w- b~, w- > j j )
output standalone ttbar_sl_pos
output standalone_cpp ttbar_sl_pos
display diagrams

generate g g > t t~, ( t > w+ b, w+ > j j ), ( t~ > w- b~, w- > l- vl~ )
output standalone ttbar_sl_neg
output standalone_cpp ttbar_sl_neg
display diagrams

import BSM_gg_hh
generate g g > h h, ( h > b b~ ), ( h > w+ w-, w+ > l+ vl, w- > l- vl~ )
output standalone hh_bbww
output standalone_cpp hh_bbww
display diagrams

import BSM_gg_hh
generate g g > h h, ( h > b b~ ), ( h > w+ w-, w+ > l+ vl, w- > j j )
output standalone hh_bbww_sl_pos
output standalone_cpp hh_bbww_sl_pos
display diagrams

import BSM_gg_hh
generate g g > h h, ( h > b b~ ), ( h > w+ w-, w+ > j j, w- > l- vl~ )
output standalone hh_bbww_sl_neg
output standalone_cpp hh_bbww_sl_neg
display diagrams
```

The fortran files are taken from `MG5_aMC_v2_3_3/models/BSM_gg_hh/Fortran`.

## Building instructions

```bash
git clone https://github.com/HEP-KBFI/hh-bbww-mem.git $CMSSW_BASE/src/hhAnalysis/bbwwMEM
```

For dependencies, see: https://github.com/HEP-KBFI/tth-analyze-mem#building-instructions

