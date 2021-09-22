using Gridap
using GridapDistributedBenchmark
using GridapDistributed
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--mesh", "-m"
        help = "Mesh generator"
        arg_type = String
        default="cartesian"
        "--subdomains", "-s"
        help = "Tuple with the # of subdomains per Cartesian direction"
        arg_type = Int64
        default=[1]
        nargs='+'
        "--partition", "-p"
        help = "Tuple with the # of cells per Cartesian direction"
        arg_type = Int64
        default=[4]
        nargs='+'
        "--num-uniform-refinements", "-r"
        help = "# of uniform refinements"
        arg_type = Int64
        default=1
        "--solver", "-k"
        help = "PETSc solver to use (gamg, mumps)"
        default = "gamg"
        "--nruns", "-n"
        help = "Number of times to repeat the experiment"
        arg_type = Int64
        default = 1
    end
    return parse_args(s)
end

parsed_args = parse_commandline()
mesh        = parsed_args["mesh"]
solver      = parsed_args["solver"]
subdomains  = Tuple(parsed_args["subdomains"])
partition   = Tuple(parsed_args["partition"])
numrefs     = parsed_args["num-uniform-refinements"]
n           = parsed_args["nruns"]


MPIPETScCommunicator() do comm
  function generate_model()
    @assert mesh=="cartesian" || mesh=="p4est"
    d = length(subdomains)
    domain = Vector{Float64}(undef, 2 * d)
    for i = 1:2:2 * d
         domain[i  ] = 0
         domain[i + 1] = 1
    end
    domain = Tuple(domain)
    coarse_model = CartesianDiscreteModel(domain,subdomains)
    tmesh = GridapDistributed.MPITimer{GridapDistributed.MPITimerModeMin}(comm.comm, "Model")
     GridapDistributed.timer_start(tmesh)
     if (mesh=="cartesian")
       model = CartesianDiscreteModel(comm, subdomains, domain, partition)
     else
       model = UniformlyRefinedForestOfOctreesDiscreteModel(comm,coarse_model,numrefs)
     end
     GridapDistributed.timer_stop(tmesh)
     model, tmesh
  end
  for i=1:n
    GridapDistributedBenchmark.run(comm, generate_model, solver)
  end
end
