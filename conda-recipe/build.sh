#!/bin/bash
 
# $PYTHON -m pip install . -vv
 
mkdir -p $PREFIX/bin/
 
cp $RECIPE_DIR/../CAMPneu.nf $RECIPE_DIR/../23SsnpAnalysis.py  $RECIPE_DIR/../environment.yml  $RECIPE_DIR/../nextflow.config $PREFIX/bin/

chmod +x $PREFIX/bin/CAMPneu.nf $PREFIX/bin/23SsnpAnalysis.py $PREFIX/bin/environment.yml $PREFIX/bin/nextflow.config
