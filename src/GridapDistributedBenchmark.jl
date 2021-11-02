module GridapDistributedBenchmark

using Gridap
using GridapDistributed
using GridapPETSc
using GridapP4est
using PartitionedArrays
const PArrays=PartitionedArrays

#   GridapDistributedPETScWrappers.C.PCFactorGetMatrix(pc[],mumpsmat)
#   MatMumpsSetIcntl(mumpsmat[],4 ,2)     # level of printing (0 to 4)
#   MatMumpsSetIcntl(mumpsmat[],28,2)     # use 1 for sequential analysis and ictnl(7) ordering,
#                                       # or 2 for parallel analysis and ictnl(29) ordering
#   MatMumpsSetIcntl(mumpsmat[],29,2)     # parallel ordering 1 = ptscotch, 2 = parmetis
#   MatMumpsSetCntl(mumpsmat[] ,3,1.0e-6)  # threshhold for row pivot detection
#   GridapDistributedPETScWrappers.C.KSPSetUp(ksp[])
# end

function petsc_gamg_options()
  """
    -ksp_type cg -ksp_rtol 1.0e-06 -ksp_atol 0.0
    -ksp_monitor -pc_type gamg -pc_gamg_type agg -pc_gamg_est_ksp_type cg
    -mg_levels_esteig_ksp_type cg -mg_coarse_sub_pc_type cholesky
    -mg_coarse_sub_pc_factor_mat_ordering_type nd -pc_gamg_process_eq_limit 50
    -pc_gamg_square_graph 9 pc_gamg_agg_nsmooths 1
  """
end


function run(parts, generate_model, solver="gamg")
  @assert solver == "gamg" || solver == "mumps"

  if solver == "gamg"
    options=petsc_gamg_options()
  else
    #options=petsc_mumps_options()
  end

  GridapPETSc.with(args=split(options)) do
    # Manufactured solution
    u(x) = x[1] + x[2]
    f(x) = -Δ(u,x)

    t = PArrays.PTimer(parts,verbose=true)
    PArrays.tic!(t)
    model=generate_model()
    PArrays.toc!(t,"Model")

    # FESpaces
    order=1
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
    PArrays.tic!(t)
    op = AffineFEOperator(a,l,U,V)
    PArrays.toc!(t,"AffineFEOperator")

    # Linear Solver
    ls = PETScLinearSolver()
    fels = LinearFESolver(ls)
    uh = solve(fels, op)
    PArrays.toc!(t,"Solve")

    # Error norms and print solution
    dΩ=Measure(trian,2*order)
    e = u-uh
    e_l2 = sum(∫(e*e)dΩ)
    PArrays.toc!(t,"L2Norm")

    tol = 1.0e-6
    map_parts(parts) do part
      if part == 1
        println("$(e_l2) < $(tol)\n")
      end
    end

    display(t)

    GridapPETSc.gridap_petsc_gc()
  end
end

function main_cartesian(parts)
  function generate_model_cartesian()
    CartesianDiscreteModel(parts, (0,1,0,1), (4,4))
  end
  run(parts,generate_model_cartesian,"gamg")
end

function main_p4est(parts)
  function generate_model_p4est()
    coarse_model = CartesianDiscreteModel((0,1,0,1), (4,4))
    UniformlyRefinedForestOfOctreesDiscreteModel(parts, coarse_model, 1)
  end
  run(parts,generate_model_p4est,"gamg")
end

#prun(main_cartesian,mpi,(1,1))
#prun(main_p4est,mpi,1)

end # module
