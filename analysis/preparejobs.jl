using Mustache
using DrWatson

jobname(args...) = replace(savename(args...;connector="_"),"="=>"_")
driverdir(args...) = normpath(projectdir("..",args...))

function convert_nc_np_to_prod(d)
  o=Dict()
  for k in keys(d)
    if k==:nc || k==:np
     o[k]=prod(d[k])
    else
     o[k]=d[k]
    end
  end
  return o
end

function jobdict(params)
  nr = params[:nr]
  np = params[:np]
  nc = params[:nc]
  fparams=convert_nc_np_to_prod(params)
  Dict(
  "q" => "normal",
  "o" => datadir(jobname(fparams,"o.txt")),
  "e" => datadir(jobname(fparams,"e.txt")),
  "walltime" => "00:15:00",
  "ncpus" => prod(np),
  "mem" => "$(prod(np)*4)gb",
  "name" => jobname(fparams),
  "nc" => nc,
  "numrefs" => haskey(params,:numrefs) ? params[:numrefs] : -1,
  "mesh" => params[:mesh],
  "solver" => params[:solver],
  "n" => prod(np),
  "np" => np,
  "nr" => nr,
  "projectdir" => driverdir(),
  "modules" => driverdir("modules.sh"),
  "title" => datadir(jobname(fparams)),
  "sysimage" => driverdir("GridapDistributedBenchmark.so")
  )
end

function generate_2d_dicts(mesh,solver,lst_nodes,lst_ls,nr=10)
   dicts = Dict[]
   d=2
   for node in lst_nodes
     px=6*node
     py=8*node
     for ls in lst_ls
       if (mesh==:cartesian)
          nx=px*ls
          ny=py*ls
          aux=Dict(:d=>2,
                   :nc=>(nx,ny),
                   :np=>(px,py),
                   :mesh=>:cartesian,
                   :solver=>solver,
                   :nr=>nr)
          push!(dicts,aux)
       else
          nx=px*(2^ls)
          ny=py*(2^ls)
          aux=Dict(:d=>2,
                   :numrefs=>ls,
                   :nc=>(nx,ny),
                   :np=>(px,py),
                   :mesh=>:p4est,
                   :solver=>solver,
                   :nr=>nr)
          push!(dicts,aux)
       end
     end
   end
   dicts
end

function generate_3d_dicts(mesh,solver,lst_nodes,lst_ls,nr=10)
  d=3
  for node in lst_nodes
    px=4*node
    py=4*node
    pz=3*node
    for ls in lst_ls
      if (mesh==:cartesian)
        nx=px*ls
        ny=py*ls
        nz=pz*ls
        aux=Dict(:d=>3,
                 :nc=>(nx,ny,nz),
                 :np=>(px,py,pz),
                 :mesh=>:cartesian,
                 :solver=>solver,
                 :nr=>nr)
        push!(dicts,aux)
      else
        nx=px*(2^ls)
        ny=py*(2^ls)
        nz=pz*(2^ls)
        aux=Dict(:d=>3,
                 :numrefs=>ls,
                 :nc=>(nx,ny,nz),
                 :np=>(px,py,pz),
                 :mesh=>:p4est,
                 :solver=>solver,
                 :nr=>nr)
      end
    end
  end
  dicts
end

dicts_cartesian=generate_2d_dicts(:cartesian,:gamg,[1],[16,32,64,128,256,512])
dicts=generate_2d_dicts(:p4est,:gamg,[1],collect(4:9))
append!(dicts,dicts_cartesian)
template = read(projectdir("jobtemplate.sh"),String)
for params in dicts
   fparams=convert_nc_np_to_prod(params)
   jobfile = datadir(jobname(fparams,"sh"))
   open(jobfile,"w") do io
     render(io,template,jobdict(params))
   end
end
