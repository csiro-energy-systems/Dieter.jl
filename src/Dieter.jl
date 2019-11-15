# module Dieter

# export df_tech, df_load, df_storage, df_avail
# export dtr

# using Reexport
using DataFrames
using CSV
# using Tables
# @reexport using DataFrames, CSV
# using Feather
# using Arrow
# using Missings
using JuMP
# @reexport using JuMP
using ProgressMeter
# @reexport using Plots
# @reexport using ColorSchemes

# import MathOptInterface
# const _MOI = MathOptInterface

include("struct.jl")
include("util.jl")
include("base.jl")
include("ev.jl")

include("run.jl")

dfDict = parse_data_to_model!(dtr)

# include("heat.jl")
# include("p2g.jl")
# include("model.jl")
# include("analysis.jl")
# include("save.jl")
# include("merge.jl")
# include("summary_plot.jl")
#
#
#
#
# export DieterModel,
# solve_model!,
# summarize_result,
# save_results,
# post_process_results,
# plot_generation_investments,
# plot_storage_investments,
# plot_curtailment,
# plot_all,
# get_result,
# c_gradient,
# create_filter_dict

# end # module
