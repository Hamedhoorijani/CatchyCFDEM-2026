#of8
#catchy8
#cfdem8


cd DEM
rm post/dump_liggghts_*
rm post/restart/*

liggghts -in in.liggghts_init


cd ../CFD
#rm -r processor*
#decomposePar
#mpirun -np 2 cfdemSolverRhoPimple -parallel
#paraview
