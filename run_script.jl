## Running code and tests for Dieter.jl

# %% Load packages
using JuMP
using CPLEX
using Dieter
using ColorSchemes
using Plots

# %% Instantiate data model
projectpath = joinpath(ENV["HOME"],"Documents/Projects/ESM/Dieter.jl/")
datapath = joinpath(projectpath,"testdata/")
rdir = joinpath(projectpath,"results_of_test_case/")

data_instance = Dict{String,Any}()
dtr = InitialiseDieterModel(DieterModel, data_instance)
dtr.settings[:datapath] = datapath
# data_instance["files"] = Dict{String,String}()
# fileDict = dtr.data["files"]

initialise_data_file_dict!(dtr)
check_data_files_exist(dtr.data["files"])

# Modify the :min_res setting over [0,100] and rerun to see comparison.
dtr.settings[:min_res] = 50
# dtr.settings[:ev] = missing
# dtr.settings[:heat] = missing
# dtr.settings[:h2] = missing

dfDict = parse_data_to_model!(dtr)

# %% Initialise model
# include("model.jl")

# Construct an optimizer factory
# solver = with_optimizer(Clp.Optimizer)
# solver = with_optimizer(Gurobi.Optimizer)
solver = JuMP.with_optimizer(CPLEX.Optimizer)
build_model!(dtr,solver)

# %% Specialised solver settings
# CPLEX:
JuMP.MOI.set(dtr.model, JuMP.MOI.RawParameter("CPX_PARAM_THREADS"), 4)        #  number of threads
JuMP.MOI.set(dtr.model, JuMP.MOI.RawParameter("CPX_PARAM_PARALLELMODE"), 1)   #   -1: Opportunistic parallel, 0: Automatic parallel, 1: Deterministic
JuMP.MOI.set(dtr.model, JuMP.MOI.RawParameter("CPX_PARAM_LPMETHOD"), 4)       #  0: auto, 1: primal simplex, 2: dual, 3: network, 4: barrier
JuMP.MOI.set(dtr.model, JuMP.MOI.RawParameter("CPX_PARAM_BARCROSSALG"), 2)    #  0: automatic, 1: primal, 2: dual
JuMP.MOI.set(dtr.model, JuMP.MOI.RawParameter("CPX_PARAM_SOLUTIONTYPE"), 2)   #  Specifies type of solution (basic or non basic) that CPLEX produces
JuMP.MOI.set(dtr.model, JuMP.MOI.RawParameter("CPX_PARAM_BAREPCOMP"), 1e-8)   # Sets the tolerance on complementarity for convergence; default: 1e-8.

# %% Solve the model and generate results
solve_model!(dtr)
generate_results!(dtr)

# %% Analysis
# include("analysis.jl")
df_summ = summarize_result(dtr,del_zeros=false)

# %% Save results to file
# include("save.jl")

save_results(dtr, rdir)

# %% Merge results with other runs
# include("merge.jl")
post_process_results(rdir)

# %% Plotting
color_dict = Dict()
default_colour = ColorSchemes.leonardo
# ColorScheme([Colors.RGB(0.0, 0.0, 0.0), Colors.RGB(1.0, 1.0, 1.0)],
      # "custom", "twotone, black and white");
      color_dict["Wind Offshore"] = [default_colour[10]];
      color_dict["Wind Onshore"] = [default_colour[20]];
      color_dict["Solar PV"] = [default_colour[30]];
      color_dict["Storages"] = [default_colour[5]];
      color_dict["Hydrogen"] = [default_colour[15]];
      color_dict["Curtailment"] = [default_colour[25]];

marker = (:hexagon, 5, 0.6, :green, stroke(3, 0.2, :black, :dot))

plot_all(rdir,color_dict,sector=:ev,marker=marker,legend=true)
