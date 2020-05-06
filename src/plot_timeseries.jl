# %% Plotting packages

using Plots
using StatsPlots
import Plots.PlotMeasures: mm
using ColorSchemes

# %% Filtering results

# Split dataframes with Tuple-type columns using split_df_tuple(), e.g.
# split_df_tuple(res[:FLOW],:Arcs,[:From,:To])
# split_df_tuple(res[:N_TECH],:Nodes_Techs,[:Nodes,:Techs])

# res = dtr.results

# df = res[:G]
df_aug = resSplit[:G] # = split_df_tuple(res[:G],:Nodes_Techs,[:Nodes,:Techs])
# dfDict = dtr.data["dataframes"]
df_nodes = dfDict["nodes"]

node2DemReg = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:DemandRegion]))

df_spatial = join(df_aug,df_nodes,on=:Nodes)
df_filter = select(df_spatial,[:Nodes,:Technologies,:DemandRegion,:Hours,:G])

dfStates = Dict{String,DataFrame}()

for reg in dtr.sets[:DemandRegions]
      dfStates[reg] = @linq df_filter |>
                        where(:DemandRegion .== reg) |>
                        select(:Nodes,:Technologies,:Hours,:G) |>
                        by([:Technologies,:Hours], Level = sum(:G))
end

#=
# df_flow = res[:FLOW]
df_flow = split_df_tuple(res[:FLOW],:Arcs,[:From,:To])

# Augment with a column showing the Demand Regions of the flow's nodes:
df_flow[!,:FromRegion] = map(x -> node2DemReg[x], df_flow[!,:From])
df_flow[!,:ToRegion] = map(x -> node2DemReg[x], df_flow[!,:To])
=#
# Inter-state flow:
df_interflow = resSplit[:INTERFLOW] # = @where(df_flow, :FromRegion .!== :ToRegion)


# %% Create plots

Hours = dtr.sets[:Hours]
L = 1:168
# L = 1:336
# L = 5000:5336

gr()
# gr(size=(5000,2000))
gr(size=(3500,600))
# plotly()
# plotly(size=(3000,600))

# DemandReg = "TAS1"
# p = plot(layout=grid(5,1, height=4*[0.1,0.1,0.1,0.1,0.1]),margin=5mm);

for (count, DemandReg) in enumerate(["NSW1", "QLD1", "VIC1", "SA1", "TAS1"])
            # Load = dtr.parameters[:Load]
      Demand = @where(dfDict["load"],:DemandRegion .== DemandReg)
      df_plot = dfStates[DemandReg]
      Techs = [Symbol(i) for i in DataFrames.unique(copy(df_plot[!,:Technologies]))]
      NumTechs = length(Techs)
      # reTechs = reshape(Techs,NumTechs,1)

      df_unstack = unstack(df_plot,:Technologies,:Level)
      p = plot()
      # @df df_all[L,:] groupedbar(HOURS[L], cols(1:NumTech),  #cols(NumTech:-1:1),
      @df df_unstack[L,Techs] groupedbar!(p, Hours[L], cols(Techs),
          # subplot=count,
          margin=10mm,
          title=DemandReg,
          xlabel="Time",
          fillalpha=0.5,linealpha=0.1,
          bar_position=:stack,
          legend=:best,  # `:none`, `:best`, `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
          color_palette=:balance) # phase delta rainbow inferno darkrainbow colorwheel

      plot!(Hours[L], [Demand[L,:Load]],label="Demand",
            # subplot=count,
            line=4, linecolour=:steelblue,
            xtickfont = font(10, "Courier"),
            xlabel="Time (hr)",
            ylabel="Generation (MW)",
            )
      # plot!(p,margin=15mm)

      df_regflow = @linq df_interflow |>
                    where(:FromRegion .== DemandReg) |>
                    by(:Hours, Level = sum(:FLOW))

      plot!(df_regflow[L,:Hours], [df_regflow[L,:Level]],label="Flow",
            # subplot=count,
            line=4, linecolour=:red,
            xtickfont = font(10, "Courier"),
            xlabel="Time (hr)",
            ylabel="Generation (MW)",
            margin=5mm
            )
      display(p)
end

# %% Misc.
# color_dict = Dict()
# default_colour = ColorSchemes.leonardo
# # ColorScheme([Colors.RGB(0.0, 0.0, 0.0), Colors.RGB(1.0, 1.0, 1.0)],
#       # "custom", "twotone, black and white");
#       color_dict["Wind Offshore"] = [default_colour[10]];
#       color_dict["Wind Onshore"] = [default_colour[20]];
#       color_dict["Solar PV"] = [default_colour[30]];
#       color_dict["Storages"] = [default_colour[5]];
#       color_dict["Hydrogen"] = [default_colour[15]];
#       color_dict["Curtailment"] = [default_colour[25]];
#
# marker = (:hexagon, 5, 0.6, :green, stroke(3, 0.2, :black, :dot))
#
# plot_all(rdir,color_dict,sector=:ev,marker=marker,legend=true)
