module Dieter

# %%
# using Reexport
using DataFrames
using DataFramesMeta
using CSV
using SQLite
# using Tables
# @reexport using DataFrames, CSV

# using Feather
# using Arrow  # review whether Arrow needs to be a dependency
# # using Missings

using JuMP
# import MathOptInterface
# const MOI = MathOptInterface
# @reexport using JuMP
using ProgressMeter
# using Plots
# @reexport using Plots
# using ColorSchemes
# @reexport using ColorSchemes

# %% Fundamental structures
include("util.jl")
include("struct.jl")
include("base.jl")

# include("spacetime.jl")

# %% Model components
include("ev.jl")
include("heat.jl")
include("p2g.jl")

# %%
const hoursInYear = 8760     # Units: hours/yr
const H2_energy_density = 39405.6 # Units : MWh / kg
const infeas_cost = 1000     # Currency (cost)

include("model.jl")


# %%
export DieterModel,
InitialiseDieterModel,
create_relation,
initialise_data_file_dict!,
check_files_exist,
parse_data_to_model!,
build_model!,
solve_model!,
generate_results!

end # module Dieter
