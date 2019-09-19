module Dieter

using Reexport
@reexport using DataFrames, CSV, Feather
using Arrow
using Missings
@reexport using JuMP
using ProgressMeter
@reexport using Plots, ColorSchemes


include("struct.jl")
include("util.jl")
include("base.jl")
include("ev.jl")
include("heat.jl")
include("p2g.jl")
include("model.jl")
include("analysis.jl")
include("save.jl")
include("merge.jl")
include("summary_plot.jl")




export DieterModel,
solve_model!,
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
