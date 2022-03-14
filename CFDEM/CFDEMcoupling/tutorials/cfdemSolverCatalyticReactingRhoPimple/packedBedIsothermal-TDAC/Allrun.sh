#!/bin/bash

export PBS_NP=2
#! manually check decomposeParDict and in.liggghts_run to be consistent with number of processors !

# run liggghts init
cd DEM
mpirun -np $PBS_NP $CFDEM_LIGGGHTS_BIN_DIR/liggghts -in in.liggghts_init > ../log.liggghts_init
cd ..

# run CFDDEM
cd DEM
cp in.liggghts_runfirst in.liggghts_run

cd ../CFD
cp -r 0.orig 0
sed -i 's/endTime         0.15;/endTime         0.10;/g' system/controlDict
blockMesh > log.blockMesh
if [ -x "$(command -v canteraToFoam)" ]; then
    canteraToFoam constant/mechanism/ocm_polimi31_srla2o3.cti gas constant/mechanism/transportProperties constant/mechanism/chem constant/mechanism/thermo -surface surface1 -transport > log.canteraToFoam
else
    echo "canteraToFoam not available - copying input files from backup"
    cp $CATCHY_ETC/backupDicts/mechanisms/ocm_polimi31_srla2o3/* constant/mechanism/
fi
decomposePar -force > log.decomposePar
#valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes mpirun -np $PBS_NP cfdemSolverCatalyticReactingRhoPimple -parallel > log.cfdemSolverCatalyticReactingRhoPimple 2>&1
mpirun -np $PBS_NP cfdemSolverCatalyticReactingRhoPimple -parallel > log.cfdemSolverCatalyticReactingRhoPimple 2>&1
reconstructPar -noLagrangian -latestTime > log.reconstructPar
#rm -r processor*

# post-processing
postProcess -func "components(U)" -latestTime > log.postProcess
postProcess -func sample -latestTime >> log.postProcess
calcMixingCup -fields "(voidfraction T p Ux Uy Uz CH4 O2 C2H4 C2H6 CO CO2)" -nPoints 20 -latestTime >> log.postProcess
cd ..

# compare cantera
python3 validate.py
