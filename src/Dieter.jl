module Dieter

# %%
# using Reexport
using DataFrames
using CSV
# using Tables
# @reexport using DataFrames, CSV
using Feather
using Arrow
# using Missings
using JuMP
# import MathOptInterface
# const MOI = MathOptInterface
# @reexport using JuMP
using ProgressMeter
using Plots
# @reexport using Plots
using ColorSchemes
# @reexport using ColorSchemes

# %%
include("util.jl")
include("struct.jl")
include("base.jl")

# %% Model components
include("ev.jl")
include("heat.jl")
include("p2g.jl")

# %%
include("model.jl")

# %%
include("analysis.jl")

# %%
include("save.jl")

# %%
include("merge.jl")

# %
include("summary_plot.jl")

# %%
export DieterModel,
InitialiseDieterModel,
initialise_data_file_dict!,
check_data_files_exist,
parse_data_to_model!,
build_model!,
solve_model!,
generate_results!,
summarize_result,
save_results,
post_process_results,
plot_generation_investments,
plot_storage_investments,
plot_curtailment,
plot_all,
get_result,
c_gradient,
create_filter_dict

end # module
