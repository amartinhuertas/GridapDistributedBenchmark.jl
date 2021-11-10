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

dicts = Dict[]
d=2
lst_nodes = collect(1:1)
lst_ls    = [16]#,32,64,128,256,512]
for node in lst_nodes
  px=6*node
  py=8*node
  for ls in lst_ls
     nx=px*ls
     ny=py*ls
     aux=Dict(:d=>2,:nc=>(nx,ny),:np=>(px,py),:mesh=>:cartesian,:solver=>:gamg,:nr=>10)
     push!(dicts,aux)
  end
end

d=3
lst_nodes = collect(1:1)
lst_ls    = [10]#,20,30,40]
for node in lst_nodes
  px=4*node
  py=4*node
  pz=3*node
  for ls in lst_ls
     nx=px*ls
     ny=py*ls
     nz=pz*ls
     aux=Dict(:d=>3,:nc=>(nx,ny,nz),:np=>(px,py,pz),:mesh=>:cartesian,:solver=>:gamg,:nr=>10)
     push!(dicts,aux)
  end
end

dicts=[Dict(:d=>2,:nc=>(10,10),:np=>(1,1),:mesh=>:cartesian,:solver=>:gamg,:nr=>10)]

# allparams=Dict()
# allparams = Dict(
#  :npx => map(i->ceil(Int,2^(i/3)*2),[0,1,3,4,5,6,7,8,9,10]),
#  :ncx => 300,
#  :nr => 6
#  )
template = read(projectdir("jobtemplate.sh"),String)
# dicts = dict_list(allparams)
for params in dicts
   fparams=convert_nc_np_to_prod(params)
   jobfile = datadir(jobname(fparams,"sh"))
   open(jobfile,"w") do io
     render(io,template,jobdict(params))
   end
end
