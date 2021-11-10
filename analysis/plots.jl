# using Plots
# using DrWatson
# using DataFrames
# using CSV

# fn = plotsdir("summary.csv")
# df = CSV.File(fn) |> DataFrame
# x = df[:,:nparts]
# y = df[:,:wct]

# function doplot()
#   plt = plot(x,y,
#      thickness_scaling = 1.2,
#      xaxis=("Number of cores",:log),
#      yaxis=("Wall time [s]",:log),
#      shape=:square,
#      label="Measured",
#      legend=:outertopright,
#      markerstrokecolor=:white,
#     )

#   x1 = first(x)
#   y1 = first(y)
#   s = y1.*(x1./x)
#   plot!(x,s,xaxis=:log, yaxis=:log, label="Ideal")
#   l = first(df.ngdofs)/25e3
#   plot!([l,l],collect(ylims(plt)),xaxis=:log, yaxis=:log,linestyle=:auto, label="25KDOFs/core")
#   plt
# end

# unicodeplots()
# display(doplot())

# gr()
# doplot()
# savefig(plotsdir("total_scaling.pdf"))
# savefig(plotsdir("total_scaling.png"))
