* HPC machine: NCI GADI
* Julia version: 1.6.2 (self-installed)
* OpenMPI 4.1.0 
* PETSc 3.15.4
* command example: 
```bash
$HOME/.julia/bin/mpiexecjl --project=/scratch/bt62/am6349/GridapDistributedBenchmark.jl/ -n 48\
    julia -O3 --check-bounds=no -e\
      'using GridapDistributedBenchmark; 
       GridapDistributedBenchmark.main(mesh=:cartesian,
                                       solver=:gamg,nc=(32, 32, 24),
                                       np=(4, 4, 3),numrefs=-1,nr=5,
                                       title="data/d_3_mesh_cartesian_nc_24576_np_48_nr_5_solver_gamg")'
```

GridapDistributedBenchmark version

* `GridapDistributedBenchmark.jl`  1bcfc6ecc2845f980c387d6bec8191b7f1467f61 branch: main 

Rest of packages details

```
(GridapDistributedBenchmark) pkg> status
     Project GridapDistributedBenchmark v0.1.0
      Status `/scratch/bt62/am6349/GridapDistributedBenchmark.jl/Project.toml`
  [c7e460c6] ArgParse v1.1.4
  [5789e2e9] FileIO v1.11.2
  [56d4f2e9] Gridap v0.17.5
  [f9701e48] GridapDistributed v0.2.0 `https://github.com/gridap/GridapDistributed.jl#master`
  [c2c8e14b] GridapP4est v0.1.0 `https://github.com/gridap/GridapP4est.jl#main`
  [bcdc36c2] GridapPETSc v0.3.0 `https://github.com/gridap/GridapPETSc.jl#exploring_petsc_garbage_collection`
  [da04e1cc] MPI v0.19.1
  [9b87118b] PackageCompiler v2.0.0
  [5a9dfac6] PartitionedArrays v0.2.7 `https://github.com/fverdugo/PartitionedArrays.jl#master`
```

```
[[GridapDistributed]]
deps = ["FillArrays", "Gridap", "LinearAlgebra", "MPI", "PartitionedArrays", "SparseArrays", "WriteVTK"]
git-tree-sha1 = "c67ac74777366e68ce0d12f9220327c3dd6b4291"
repo-rev = "master"
repo-url = "https://github.com/gridap/GridapDistributed.jl"
uuid = "f9701e48-63b3-45aa-9a63-9bc6c271f355"
version = "0.2.0"

[[GridapP4est]]
deps = ["ArgParse", "FillArrays", "Gridap", "GridapDistributed", "MPI", "P4est_wrapper", "PartitionedArrays", "Test"]
git-tree-sha1 = "4ded9f0dfe1031fba8730f652b5e1d3a3b3ec124"
repo-rev = "main"
repo-url = "https://github.com/gridap/GridapP4est.jl"
uuid = "c2c8e14b-f5fd-423d-9666-1dd9ad120af9"
version = "0.1.0"

[[GridapPETSc]]
deps = ["Gridap", "GridapDistributed", "Libdl", "LinearAlgebra", "MPI", "PETSc_jll", "PartitionedArrays", "SparseArrays", "SparseMatricesCSR"]
git-tree-sha1 = "dd671660a9e7d5f75a48cea11d3d89f9ad47e85a"
repo-rev = "exploring_petsc_garbage_collection"
repo-url = "https://github.com/gridap/GridapPETSc.jl"
uuid = "bcdc36c2-0c3e-11ea-095a-c9dadae499f1"
version = "0.3.0"

[[PartitionedArrays]]
deps = ["Distances", "IterativeSolvers", "LinearAlgebra", "MPI", "Printf", "SparseArrays", "SparseMatricesCSR"]
git-tree-sha1 = "d67ba9b3b6c4ff0024c6e5d9aa500f0e3a40cba2"
repo-rev = "master"
repo-url = "https://github.com/fverdugo/PartitionedArrays.jl"
uuid = "5a9dfac6-5c52-46f7-8278-5e2210713be9"
version = "0.2.7"
```
