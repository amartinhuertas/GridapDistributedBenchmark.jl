* HPC machine:    NCI GADI
* Julia version:  1.4.2 (self-installed)
* command example: 
```bash
mpirun -np 768 /home/565/am6349/julia-1.4.2/bin/julia \
-J /home/565/am6349/GridapDistributedBenchmark.jl/GridapDistributedBenchmark.so \ 
--project=/home/565/am6349/GridapDistributedBenchmark.jl \ 
/home/565/am6349/GridapDistributedBenchmark.jl/src/test.jl \ 
-n 10 -p 12288 16384 -s 24 32
```
* GridapDistributed.jl           commit ID: 0828c7924eaf77b4427a0d6e92e243c50af0ce7f branch: mpi_petsc_communicator
* PETSc.jl                       commit ID: d757fe04a8c1d64ad12b15cc19715418cee0fe9b branch: uptodate1.4
* GridapDistributedBenchmark.jl  commit ID: c495b18688d8d2da8619fda9db716cd989a3cbb5 branch: master
