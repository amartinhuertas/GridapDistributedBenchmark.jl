using GridapDistributedBenchmark
using GridapDistributed
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
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
        "--nruns", "-n"
        help = "Number of times to repeat the experiment"
        arg_type = Int64
        default = 1
    end
    return parse_args(s)
end

parsed_args = parse_commandline()
subdomains = Tuple(parsed_args["subdomains"])
partition = Tuple(parsed_args["partition"])
n = parsed_args["nruns"]

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

MPIPETScCommunicator() do comm
  for i=1:n
    GridapDistributedBenchmark.run(comm, subdomains, partition; options...)
  end
end
