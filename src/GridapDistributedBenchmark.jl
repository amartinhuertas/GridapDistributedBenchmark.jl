module GridapDistributedBenchmark

using Gridap
using GridapDistributed
using GridapPETSc
using GridapP4est
using PartitionedArrays
using FileIO
using MPI
const PArrays=PartitionedArrays


function petsc_gamg_options()
  """
    -ksp_type cg -ksp_rtol 1.0e-06 -ksp_atol 0.0
    -ksp_monitor -pc_type gamg -pc_gamg_type agg -pc_gamg_est_ksp_type cg
    -mg_levels_esteig_ksp_type cg -mg_coarse_sub_pc_type cholesky
    -mg_coarse_sub_pc_factor_mat_ordering_type nd -pc_gamg_process_eq_limit 50
    -pc_gamg_square_graph 9 pc_gamg_agg_nsmooths 1
  """
end

function mytic!(t,comm)
  MPI.Barrier(comm)
  PArrays.tic!(t)
end

function _from_setup_fe_space_to_the_end(t,model,order=1)
  # Manufactured solution
  u(x) = x[1] + x[2]
  f(x) = -Δ(u,x)

  comm = model.models.comm

  # FESpaces
  mytic!(t,comm)
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(model,reffe,dirichlet_tags="boundary")
  U = TrialFESpace(V, u)
  PArrays.toc!(t,"FESpaces")

  trian=Triangulation(model)
  dΩ=Measure(trian,2*(order+1))
  function a(u,v)
    ∫(∇(v)⋅∇(u))dΩ
  end
  function l(v)
    ∫(v*f)dΩ
  end

  # FE Affine Operator
  mytic!(t,comm)
  op = AffineFEOperator(a,l,U,V)
  PArrays.toc!(t,"AffineFEOperator")

  # Linear Solver
  mytic!(t,comm)
  ls = PETScLinearSolver()
  fels = LinearFESolver(ls)
  uh = solve(fels, op)
  PArrays.toc!(t,"Solve")

  # Error norms and print solution
  mytic!(t,comm)
  dΩ=Measure(trian,2*order)
  e = u-uh
  e_l2 = sum(∫(e*e)dΩ)
  PArrays.toc!(t,"L2Norm")

  map_main(get_part_ids(model.models)) do part
    println("$(e_l2)\n")
  end

  ngdofs  = length(V.gids)
  ngcells = length(model.gids)
  ngdofs, ngcells, e_l2
end


function generate_model_cartesian(parts,subdomains,partition)
  d = length(subdomains)
  domain = Vector{Float64}(undef, 2*d)
  for i = 1:2:2 * d
    domain[i]=0
    domain[i+1]=1
  end
  domain = Tuple(domain)
  model = CartesianDiscreteModel(parts, domain, partition)
end

function main_cartesian(parts,subdomains,partition,title,ir,order=1)
  t = PArrays.PTimer(parts,verbose=true)
  PArrays.tic!(t)
  model=generate_model_cartesian(parts,subdomains,partition)
  PArrays.toc!(t,"Model")
  ngcells, ngdofs, enorm = _from_setup_fe_space_to_the_end(t,model,order)
  display(t)
  nparts = length(parts)
  map_main(t.data) do data
    out = Dict{String,Any}()
    merge!(out,data)
    out["d"] = length(subdomains)
    out["mesh"] = "cartesian"
    out["enorm"] = enorm
    out["nparts"] = nparts
    out["ngdofs"] = ngdofs
    out["ngcells"] = ngcells
    out["nc"] = partition
    out["np"] = subdomains
    out["ir"] = ir
    save("$title.bson",out)
  end
end

function generate_model_p4est(parts,subdomains,numrefs)
  d = length(subdomains)
  domain = Vector{Float64}(undef, 2*d)
  for i = 1:2:2 * d
    domain[i]=0
    domain[i+1]=1
  end
  domain = Tuple(domain)
  coarse_model = CartesianDiscreteModel(domain,subdomains)
  model = UniformlyRefinedForestOfOctreesDiscreteModel(parts,coarse_model,numrefs)
end

function main_p4est(parts,subdomains,numrefs,title,ir,order=1)
  t = PArrays.PTimer(parts,verbose=true)
  PArrays.tic!(t)
  model=generate_model_p4est(parts,subdomains,numrefs)
  PArrays.toc!(t,"Model")
  ngcells, ngdofs, enorm = _from_setup_fe_space_to_the_end(t,model,order)
  display(t)
  nparts = length(parts)
  map_main(t.data) do data
    out = Dict{String,Any}()
    merge!(out,data)
    out["d"] = length(subdomains)
    out["mesh"] = "p4est"
    out["enorm"] = enorm
    out["nparts"] = nparts
    out["ngdofs"] = ngdofs
    out["ngcells"] = ngcells
    out["nc"] = numrefs
    out["np"] = subdomains
    out["ir"] = ir
    save("$title.bson",out)
  end
end

########

function main(;
  mesh::Symbol,
  solver::Symbol,
  np::Tuple,
  nr::Integer,
  title::AbstractString,
  nc::Tuple=(-1,-1),
  numrefs::Integer=-1,
  k::Integer=1,
  verbose::Bool=true)

  mesh   in (:cartesian,:p4est) || throw(ArgumentError("mesh should be :cartesian or :p4est"))
  solver in (:gamg,:mumps)      || throw(ArgumentError("solver should be :gamg or :mumps"))

  # Process parameters of mesh
  if mesh == :p4est
    numrefs>=1 || throw(ArgumentError("numrefs should be larger or equal than 1"))
  else
    length(np) == length(nc) || throw(ArgumentError("np and nc must be of same length"))
    all(nc .> 0) || throw(ArgumentError("all values in nc should be larger than 0"))
  end

  if solver == :gamg
    options=petsc_gamg_options()
  else
    #options=petsc_mumps_options()
  end

  if mesh == :p4est
    prun(mpi,prod(np)) do parts
      for ir in 1:nr
        GridapPETSc.with(args=split(options)) do
           str_r   = lpad(ir,ceil(Int,log10(nr)),'0')
           title_r = "$(title)_ir$(str_r)"
           main_p4est(parts,np,numrefs,title_r,ir,1)
           GridapPETSc.gridap_petsc_gc()
        end
      end
    end
  else
    prun(mpi,np) do parts
      for ir in 1:nr
        GridapPETSc.with(args=split(options)) do
           str_r   = lpad(ir,ceil(Int,log10(nr)),'0')
           title_r = "$(title)_ir$(str_r)"
           main_cartesian(parts,np,nc,title_r,ir,1)
           GridapPETSc.gridap_petsc_gc()
        end
      end
    end
  end
end

#   GridapDistributedPETScWrappers.C.PCFactorGetMatrix(pc[],mumpsmat)
#   MatMumpsSetIcntl(mumpsmat[],4 ,2)     # level of printing (0 to 4)
#   MatMumpsSetIcntl(mumpsmat[],28,2)     # use 1 for sequential analysis and ictnl(7) ordering,
#                                       # or 2 for parallel analysis and ictnl(29) ordering
#   MatMumpsSetIcntl(mumpsmat[],29,2)     # parallel ordering 1 = ptscotch, 2 = parmetis
#   MatMumpsSetCntl(mumpsmat[] ,3,1.0e-6)  # threshhold for row pivot detection
#   GridapDistributedPETScWrappers.C.KSPSetUp(ksp[])
# end
end
