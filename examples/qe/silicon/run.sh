#!/bin/bash -l

PW=pw.x
PH=ph.x
P2Y=p2y
YAMBO=yambo
LELPH=../../../../src/lelphc

#run scf
cd scf
$PW < scf.in | tee scf.out
cd ..
# run phonons
cd ph
cp -r ../scf/si.* .
$PH < ph.in | tee ph.out
cd ..
# run nscf
cd nscf
cp -r ../scf/si.* .
$PW < nscf.in | tee nscf.out
cd si.save
## create save dir
$P2Y
$YAMBO
cd ../..
## now create ph_save
cd ph
$LELPH -pp --code=qe -F ph.in
cd ..
cd elph
$LELPH -F elph.in
cd ..
