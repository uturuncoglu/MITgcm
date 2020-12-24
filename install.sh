#!/bin/bash

HOME=$PWD

rm -rf $HOME/verification/MED12_ALADIN_75_COUPLED/build
mkdir $HOME/verification/MED12_ALADIN_75_COUPLED/build
cd $HOME/verification/MED12_ALADIN_75_COUPLED/build

$HOME/tools/genmake2 \
	-rootdir=$HOME \
	-of=$HOME/opt/ictp_argo.cfg \
	-mods=$HOME/verification/MED12_ALADIN_75_COUPLED/code \
	-make=gmake -mpi

make depend
make
