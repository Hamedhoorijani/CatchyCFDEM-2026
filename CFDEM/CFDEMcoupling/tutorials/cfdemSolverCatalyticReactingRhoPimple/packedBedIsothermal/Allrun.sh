#!/bin/bash

export PBS_NP=4 
#! manually check decomposeParDict and in.liggghts_run to be consistent with number of processors !

# run liggghts init
cd DEM
mpirun -np $PBS_NP $CFDEM_LIGGGHTS_BIN_DIR/liggghts -in in.liggghts_init > ../log.liggghts_init
cd ..

# run CFDDEM
cd CFD
cp -r 0.orig 0
blockMesh > log.blockMesh
if [ -x "$(command -v canteraToFoam)" ]; then
    canteraToFoam constant/mechanism/ocm_polimi31_srla2o3.cti gas constant/mechanism/transportProperties constant/mechanism/chem constant/mechanism/thermo -surface surface1 -transport > log.canteraToFoam
else
    echo "canteraToFoam not available - copying input files from backup"
    cp $CATCHY_ETC/backupDicts/mechanisms/ocm_polimi31_srla2o3/* constant/mechanism/
fi
decomposePar -force > log.decomposePar
mpirun -np $PBS_NP cfdemSolverCatalyticReactingRhoPimple -parallel > log.cfdemSolverCatalyticReactingRhoPimple
reconstructPar -noLagrangian -latestTime > log.reconstructPar
rm -r processor*

# post-processing
postProcess -func "components(U)" -latestTime > log.postProcess
postProcess -func sample -latestTime >> log.postProcess
calcMixingCup -fields "(voidfraction T p Ux Uy Uz CH4 O2 C2H4 C2H6 CO CO2)" -nPoints 15 -latestTime >> log.postProcess
cd ..

# compare cantera
python surf_pfr.py
python plot_comparison.py
