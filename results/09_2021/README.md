HPC machine:    NCI GADI
Julia version:  1.6.2 (self-installed)
command example: 
mpirun --report-bindings -np 3888 /home/565/am6349/julia-1.6.2/bin/julia -J /home/565/am6349/GridapDistributedBenchmark.jl/GridapDistributedBenchmark.so --project=/home/565/am6349/GridapDistributedBenchmark.jl /home/565/am6349/GridapDistributedBenchmark.jl/src/test.jl -n 10 -p 378 504 -s 54 72 -k gamg -m p4est -r 7

GridapDistributed.jl           v0.1.0 (registered) 
GridapDistributedBenchmark.jl  afc1becfee1aa97ec93778d22ba022976a1fd3ba branch: main 
