using Gridap
using GridapPETSc
using GridapDistributed
using GridapDistributedBenchmark
using ArgParse
using GridapP4est
using PartitionedArrays
const PArrays=PartitionedArrays

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

function main_cartesian(parts)
  function generate_model_cartesian()
    d = length(subdomains)
    domain = Vector{Float64}(undef, 2*d)
    for i = 1:2:2 * d
      domain[i]=0
      domain[i+1]=1
    end
    domain = Tuple(domain)
    model = CartesianDiscreteModel(parts, domain, partition)
  end
  GridapDistributedBenchmark.run(parts,generate_model_cartesian,"gamg")
end

function main_p4est(parts)
  function generate_model_p4est()
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
  GridapDistributedBenchmark.run(parts,generate_model_p4est,"gamg")
end

@assert mesh=="cartesian" || mesh=="p4est"
if (mesh=="cartesian")
  for i=1:n
    prun(main_cartesian,mpi,subdomains)
  end
else
  for i=1:n
    prun(main_p4est,mpi,prod(subdomains))
  end
end
