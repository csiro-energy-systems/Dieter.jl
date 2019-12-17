To run the model you can use this script:
```julia
using Dieter
using JuMP

# select your favourite solver
using Gurobi
solver = with_optimizer(Gurobi.Optimizer)
# using CPLEX
# solver = with_optimizer(CPLEX.Optimizer)
# using Clp
# solver = with_optimizer(Clp.Optimizer)

datapath = joinpath(pwd(),"testdata")

# set the Dieter parameters
dtr = InitialiseDieterModel(DieterModel, Dict{String,Any}())
dtr.settings[:datapath] = datapath
dtr.settings[:min_res] = 50

initialise_data_file_dict!(dtr)
check_data_files_exist(dtr.data["files"])

dfDict = parse_data_to_model!(dtr)

build_model!(dtr,solver)
solve_model!(dtr)

generate_results!(dtr)
save_results(dtr, "results_of_test_case")
```
See the file `run_script.jl` for a fuller example workflow.
