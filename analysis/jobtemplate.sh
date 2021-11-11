#!/bin/bash
#PBS -P bt62
#PBS -q {{q}} 
#PBS -l walltime={{walltime}}
#PBS -l ncpus={{ncpus}}
#PBS -l mem={{mem}}
#PBS -N {{{name}}}
#PBS -l wd
#PBS -o {{{o}}}
#PBS -e {{{e}}} 
#PBS -l software=GridapDistributed.jl

source {{{modules}}}

$HOME/.julia/bin/mpiexecjl --project={{{projectdir}}} -n {{n}}\
    julia -J {{{sysimage}}} -O3 --check-bounds=no -e\
      'using GridapDistributedBenchmark; GridapDistributedBenchmark.main(mesh=:{{mesh}},solver=:{{solver}},nc={{nc}},np={{np}},numrefs={{numrefs}},nr={{nr}},title="{{{title}}}")'

