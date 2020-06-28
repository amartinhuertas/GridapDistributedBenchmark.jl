module GridapDistributedBenchmark

using Gridap
using Gridap.FESpaces
using GridapDistributed
using PETSc

function run(assembly_strategy::AbstractString, subdomains=(1, 1), cells=(4, 4); kwargs...)
    T = Float64
    vector_type = PETSc.Vec{T}
    matrix_type = PETSc.Mat{T}

  # Manufactured solution
    u(x) = x[1] + x[2]
    f(x) = -Δ(u)(x)

  # Discretization
    comm = MPIPETScCommunicator()

    d = length(subdomains)
    domain = Vector{Float64}(undef, 2 * d)
    for i = 1:2:2 * d
        domain[i  ] = 0
        domain[i + 1] = 1
    end
    domain = Tuple(domain)

    tmesh = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "CartesianDiscreteModel")
    GridapDistributed.timer_start(tmesh)
    model = CartesianDiscreteModel(comm, subdomains, domain, cells)
    GridapDistributed.timer_stop(tmesh)

    tfespaces = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "FESpaces")
    GridapDistributed.timer_start(tfespaces)
  # FE Spaces
    order = 1
    V = FESpace(
    vector_type,
    valuetype=Float64,
    reffe=:Lagrangian,
    order=order,
    model=model,
    conformity=:H1,
    dirichlet_tags="boundary",
  )
    U = TrialFESpace(V, u)
    GridapDistributed.timer_stop(tfespaces)


    if (assembly_strategy == "RowsComputedLocally")
        strategy = RowsComputedLocally(V; global_dofs=false)
    elseif (assembly_strategy == "OwnedCellsStrategy")
        strategy = OwnedCellsStrategy(model, V; global_dofs=false)
    else
        @assert false "Unknown AssemblyStrategy: $(assembly_strategy)"
    end

  # Terms in the weak form
    terms = DistributedData(model, strategy) do part, (model, gids), strategy
        trian = Triangulation(strategy, model)
        degree = 2 * order
        quad = CellQuadrature(trian, degree)
        a(u, v) = ∇(v) ⋅ ∇(u)
        l(v) = v * f
        t1 = AffineFETerm(a, l, trian, quad)
        (t1,)
    end

  # # Assembler
    assem = SparseMatrixAssembler(
    matrix_type,
    vector_type,
    U,
    V,
    strategy,
  )

    taffineoperator = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "AffineFEOperator")
    GridapDistributed.timer_start(taffineoperator)
  # FE solution
    op = AffineFEOperator(assem, terms)
    GridapDistributed.timer_stop(taffineoperator)


    tsolve = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "Solve")
    GridapDistributed.timer_start(tsolve)
    ls = PETScLinearSolver(
    Float64; kwargs...
  )
    fels = LinearFESolver(ls)
    uh = solve(fels, op)
    GridapDistributed.timer_stop(tsolve)


    tl2norm = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "L2Norm")
    GridapDistributed.timer_start(tl2norm)
  # Error norms and print solution
    sums = DistributedData(model, uh) do part, (model, gids), uh
        trian = Triangulation(model)
        owned_trian = remove_ghost_cells(trian, part, gids)
        owned_quad = CellQuadrature(owned_trian, 2 * order)
        owned_uh = restrict(uh, owned_trian)
    # writevtk(owned_trian, "results_$part", cellfields = ["uh" => owned_uh])
        e = u - owned_uh
        l2(u) = u * u
        sum(integrate(l2(e), owned_trian, owned_quad))
    end
    e_l2 = sum(gather(sums))
    GridapDistributed.timer_stop(tl2norm)

    tol = 1.0e-6
    if (i_am_master(comm)) println("$(e_l2) < $(tol)\n") end

    GridapDistributed.timer_report(tmesh)
    GridapDistributed.timer_report(tfespaces, false)
    GridapDistributed.timer_report(taffineoperator, false)
    GridapDistributed.timer_report(tsolve, false)
    GridapDistributed.timer_report(tl2norm, false)
end

end # module
