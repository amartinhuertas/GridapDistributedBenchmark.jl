module GridapDistributedBenchmark

using Gridap
using GridapDistributed
using GridapDistributedPETScWrappers

function PCFactorSetMatSolverType(arg1::GridapDistributedPETScWrappers.C.PC{Float64},
                                  arg2::Union{String,Cstring,Symbol,Array{UInt8},Ptr{UInt8}})
   err = ccall((:PCFactorSetMatSolverType,GridapDistributedPETScWrappers.C.petscRealDouble),
                 GridapDistributedPETScWrappers.C.PetscErrorCode,
                 (GridapDistributedPETScWrappers.C.PC{Float64},Cstring),
                 arg1,arg2)
   return err
end
function PCFactorSetUpMatSolverType(arg1::GridapDistributedPETScWrappers.C.PC{Float64})
   err = ccall((:PCFactorSetUpMatSolverType,GridapDistributedPETScWrappers.C.petscRealDouble),
                 GridapDistributedPETScWrappers.C.PetscErrorCode,
                 (GridapDistributedPETScWrappers.C.PC{Float64},),
               arg1)
   return err
end
function MatMumpsSetIcntl(arg1::GridapDistributedPETScWrappers.C.Mat{Float64},
                          arg2::GridapDistributedPETScWrappers.C.PetscInt,
                          arg3::GridapDistributedPETScWrappers.C.PetscInt)
    err = ccall((:MatMumpsSetIcntl,GridapDistributedPETScWrappers.C.petscRealDouble),
                 GridapDistributedPETScWrappers.C.PetscErrorCode,
                 (GridapDistributedPETScWrappers.C.Mat{Float64},
                  GridapDistributedPETScWrappers.C.PetscInt,
                  GridapDistributedPETScWrappers.C.PetscInt),
                 arg1,arg2,arg3)
    return err
end
function MatMumpsSetCntl(arg1::GridapDistributedPETScWrappers.C.Mat{Float64},
                         arg2::GridapDistributedPETScWrappers.C.PetscInt,
                         arg3::Cdouble)
    err = ccall((:MatMumpsSetIcntl,GridapDistributedPETScWrappers.C.petscRealDouble),
                 GridapDistributedPETScWrappers.C.PetscErrorCode,
                 (GridapDistributedPETScWrappers.C.Mat{Float64},
                  GridapDistributedPETScWrappers.C.PetscInt,
                  Cdouble),
                 arg1,arg2,arg3)
    return err
end

function solve_linear_system_petsc_mumps_low_level_API(op)
  A=op.op.matrix
  b=op.op.vector
  ksp=Ref{GridapDistributedPETScWrappers.C.KSP{Float64}}()
  pc=Ref{GridapDistributedPETScWrappers.C.PC{Float64}}()
  mumpsmat=Ref{GridapDistributedPETScWrappers.C.Mat{Float64}}()

  GridapDistributedPETScWrappers.C.KSPCreate(comm(A),ksp)
  GridapDistributedPETScWrappers.C.KSPSetOperators(ksp[],A.p,A.p)
  GridapDistributedPETScWrappers.C.KSPSetType(ksp[],GridapDistributedPETScWrappers.C.KSPPREONLY)
  GridapDistributedPETScWrappers.C.KSPGetPC(ksp[],pc)

  # If system is SPD use the following two calls
  GridapDistributedPETScWrappers.C.PCSetType(pc[],GridapDistributedPETScWrappers.C.PCCHOLESKY)
  GridapDistributedPETScWrappers.C.MatSetOption(A.p,
                                                GridapDistributedPETScWrappers.C.MAT_SPD,GridapDistributedPETScWrappers.C.PETSC_TRUE);
  # Else ... use only the following one
  # GridapDistributedPETScWrappers.C.PCSetType(pc,GridapDistributedPETScWrappers.C.PCLU)

  PCFactorSetMatSolverType(pc[],GridapDistributedPETScWrappers.C.MATSOLVERMUMPS)
  PCFactorSetUpMatSolverType(pc[])
  GridapDistributedPETScWrappers.C.PCFactorGetMatrix(pc[],mumpsmat)
  MatMumpsSetIcntl(mumpsmat[],4 ,2)     # level of printing (0 to 4)
  MatMumpsSetIcntl(mumpsmat[],28,2)     # use 1 for sequential analysis and ictnl(7) ordering,
                                      # or 2 for parallel analysis and ictnl(29) ordering
  MatMumpsSetIcntl(mumpsmat[],29,2)     # parallel ordering 1 = ptscotch, 2 = parmetis
  MatMumpsSetCntl(mumpsmat[] ,3,1.0e-6)  # threshhold for row pivot detection
  GridapDistributedPETScWrappers.C.KSPSetUp(ksp[])

  x=copy(b)
  GridapDistributedPETScWrappers.C.KSPSolve(ksp[], b.p, x.p)

  GridapDistributedPETScWrappers.C.KSPDestroy(ksp)

  uh = FEFunction(op.trial,x)
end

function generate_petsc_gamg_options()
  options=Dict{Symbol,Any}()
  options[:ksp_type]="cg"
  options[:ksp_rtol]=1.0e-06
  options[:ksp_atol]=0.0
  options[:ksp_monitor]=""
  options[:pc_type]="gamg"
  options[:pc_gamg_type]="agg"
  options[:pc_gamg_est_ksp_type]="cg"
  options[:mg_levels_esteig_ksp_type]="cg"
  options[:mg_coarse_sub_pc_type]="cholesky"
  options[:mg_coarse_sub_pc_factor_mat_ordering_type]="nd"
  options[:pc_gamg_process_eq_limit]=50
  options[:pc_gamg_square_graph]=9
  options[:pc_gamg_agg_nsmooths]=1
  options
end

function run(comm, generate_model, solver="gamg")
  @assert solver == "gamg" || solver == "mumps"

  # Manufactured solution
  u(x) = x[1] + x[2]
  f(x) = -Δ(u)(x)

  model,tmesh=generate_model()

  # FESpaces
  tfespaces = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "FESpaces")
  GridapDistributed.timer_start(tfespaces)
  order=1
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = FESpace(model=model,
              reffe=reffe,
              conformity=:H1,
              dirichlet_tags="boundary")
  U = TrialFESpace(V, u)
  GridapDistributed.timer_stop(tfespaces)

  trian=Triangulation(model)
  dΩ=Measure(trian,2*(order+1))
  function a(u,v)
    ∫(∇(v)⋅∇(u))dΩ
  end
  function l(v)
    ∫(v*f)dΩ
  end

  # FE Affine Operator
  taffineoperator = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "AffineFEOperator")
  GridapDistributed.timer_start(taffineoperator)
  op = AffineFEOperator(a,l,U,V)
  GridapDistributed.timer_stop(taffineoperator)

  # Linear system solution
  tsolve = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "Solve")
  GridapDistributed.timer_start(tsolve)
  if solver == "gamg"

    ls = PETScLinearSolver(Float64; generate_petsc_gamg_options()...)
    fels = LinearFESolver(ls)
    uh = solve(fels, op)
    GridapDistributedPETScWrappers.PetscDestroy(ls.ksp)
  elseif solver == "mumps"
    uh = solve_linear_system_petsc_mumps_low_level_API(op)
  end
  GridapDistributed.timer_stop(tsolve)

  GridapDistributedPETScWrappers.PetscDestroy(op.op.matrix)
  GridapDistributedPETScWrappers.PetscDestroy(op.op.vector)

  # Error norms and print solution
  tl2norm = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "L2Norm")
  GridapDistributed.timer_start(tl2norm)
  trian=Triangulation(OwnedCells,model)
  dΩ=Measure(trian,2*order)
  e = u-uh
  e_l2 = sum(∫(e*e)dΩ)
  GridapDistributed.timer_stop(tl2norm)

  tol = 1.0e-6
  if (i_am_master(comm)) println("$(e_l2) < $(tol)\n") end

  GridapDistributed.timer_report(tmesh)
  GridapDistributed.timer_report(tfespaces, false)
  GridapDistributed.timer_report(taffineoperator, false)
  GridapDistributed.timer_report(tsolve, false)
  GridapDistributed.timer_report(tl2norm, false)
end

MPIPETScCommunicator() do comm
   function generate_model_cartesian()
     tmesh = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "Model")
     GridapDistributed.timer_start(tmesh)
     model = CartesianDiscreteModel(comm, (1,1), (0,1,0,1), (4,4))
     GridapDistributed.timer_stop(tmesh)
     model, tmesh
   end
   function generate_model_p4est()
    tmesh = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "Model")
    GridapDistributed.timer_start(tmesh)
    coarse_model = CartesianDiscreteModel((0,1,0,1), (4,4))
    model = UniformlyRefinedForestOfOctreesDiscreteModel(comm, coarse_model, 1)
    GridapDistributed.timer_stop(tmesh)
    model, tmesh
   end
   run(comm,generate_model_cartesian, "gamg")
   run(comm,generate_model_p4est    , "gamg")
end

end # module
