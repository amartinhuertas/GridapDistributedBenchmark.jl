module GridapDistributedBenchmark

using Gridap
using GridapDistributed

function run(comm, subdomains=(1, 1), cells=(4, 4); kwargs...)

  # Manufactured solution
  u(x) = x[1] + x[2]
  f(x) = -Δ(u)(x)

  d = length(subdomains)
  domain = Vector{Float64}(undef, 2 * d)
  for i = 1:2:2 * d
    domain[i  ] = 0
    domain[i + 1] = 1
  end
  domain = Tuple(domain)

  # Geometry
  tmesh = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "CartesianDiscreteModel")
  GridapDistributed.timer_start(tmesh)
  model = CartesianDiscreteModel(comm, subdomains, domain, cells)
  GridapDistributed.timer_stop(tmesh)


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
  ls = PETScLinearSolver(Float64; kwargs...)
  fels = LinearFESolver(ls)
  uh = solve(fels, op)
  GridapDistributed.timer_stop(tsolve)

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

# options=Dict{Symbol,Any}()
# options[:ksp_type]="cg"
# options[:ksp_rtol]=1.0e-06
# options[:ksp_atol]=0.0
# options[:ksp_monitor]=""
# options[:pc_type]="gamg"
# options[:pc_gamg_type]="agg"
# options[:pc_gamg_est_ksp_type]="cg"
# options[:mg_levels_esteig_ksp_type]="cg"
# options[:mg_coarse_sub_pc_type]="cholesky"
# options[:mg_coarse_sub_pc_factor_mat_ordering_type]="nd"
# options[:pc_gamg_process_eq_limit]=50
# options[:pc_gamg_square_graph]=9
# options[:pc_gamg_agg_nsmooths]=1

# # MPIPETScCommunicator() do comm
# #    run(comm;options...)
# # end

end # module
