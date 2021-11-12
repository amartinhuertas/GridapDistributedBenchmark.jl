using Plots

function doplot(file,ls,mesh)
  ls = string(ls)
  f = open(file)

  x=Float64[]
  ymod=Float64[]
  yfes=Float64[]
  yope=Float64[]
  ysol=Float64[]
  ynrm=Float64[]
  for line in eachline(f)
     words=split(line)
     if words[1] == "STRATEGY"
      readline(f)
      continue
     end
     if words[4] == ls
       append!(x,parse(Float64,words[2]))
       append!(ymod,parse(Float64,words[5]))
       append!(yfes,parse(Float64,words[6]))
       append!(yope,parse(Float64,words[7]))
       append!(ysol,parse(Float64,words[8]))
       append!(ynrm,parse(Float64,words[9]))
     end
  end
  m=maximum([maximum(ymod),maximum(yfes),maximum(yope),maximum(ynrm)])
  if mesh == "P4EST"
    ls = 2^parse(Int,ls)
  end
  plt = plot(x,ymod,
             thickness_scaling = 1.2,
             xaxis=("Number of cores"),
             yaxis=("Wall time [s]"),
             ylims = (0,1.2*m),
             shape=:auto,
             label="Model ($(mesh))",
             title="Gridap 0.1.0 LS=$(ls)",
             legend=:outertopright,
             markerstrokecolor=:white)
  plot!(x,yfes,label="FESpace",shape=:auto)
  plot!(x,yope,label="FEOperator",shape=:auto)
  plot!(x,ysol,label="Solver",shape=:auto)
  plot!(x,ynrm,label="L2Norm",shape=:auto)
  close(f)
  plt
end

function saveplot2D(file,ls,mesh)
  gr()
  doplot(file,ls,mesh)
  savefig("weak_scaling_0.1.0_2D_$(ls)_$(mesh).pdf")
end

function saveplot3D(file,ls,mesh)
  gr()
  doplot(file,ls,mesh)
  savefig("weak_scaling_0.1.0_3D_$(ls)_$(mesh).pdf")
end


meshes = ["CARTESIAN","P4EST"]
ls_lst_cartesian  = [2^i for i=4:9]
ls_lst_p4est      = [i for i=4:9]
for mesh in meshes
  if mesh == "CARTESIAN"
    ls_lst=ls_lst_cartesian
    file="report_2D_cartesian.txt"
  elseif mesh == "P4EST"
    ls_lst=ls_lst_p4est
    file="report_2D_p4est.txt"
  end
  for ls in ls_lst
    saveplot2D(file,ls,mesh)
  end
end

meshes = ["CARTESIAN"]
ls_lst_cartesian  = [2^i for i=2:5]
ls_lst_p4est      = [i for i=2:5]
for mesh in meshes
  if mesh == "CARTESIAN"
    ls_lst=ls_lst_cartesian
    file="report_3D_cartesian.txt"
  elseif mesh == "P4EST"
    ls_lst=ls_lst_p4est
    file="report_3D_p4est_exp1.txt"
  end
  for ls in ls_lst
    saveplot3D(file,ls,mesh)
  end
end
