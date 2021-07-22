
# DIETER.jl

[DIETER](http://www.diw.de/dieter) is the _Dispatch and Investment Evaluation Tool with Endogenous Renewables_, an open energy system model developed by the German Institute for Economic Research ([DIW Berlin](https://www.diw.de/en)). The model is structured as a cost-minimising [linear program](https://en.wikipedia.org/wiki/Linear_programming) that when solved to optimality will produce an hourly operation profile of an interconnected transmission system over a single year, including aspects such as
- generator power output,
- storage technology charge and discharge,
- transmission flows,
- variable renewable energy utilisation, and
- hydrogen electrolyser operation.

The model also decides on cost-optimal capacity expansion of relevant model technologies relative to the year under consideration. A higher level module for programmatically constructing the required data model is enabled by [**STABLE.jl**](https://github.com/jd-foster/STABLE.jl).

This code repository holds a derivative version of DIETER written using the [Julia](https://julialang.org/) language and the [JuMP](https://github.com/jump-dev/JuMP.jl) algebraic modelling language (a Julia package). While the original version of DIETER is implemented in [GAMS](https://www.gams.com/), Julia/JuMP enables a free and open-source implementation. Linear programming (LP) solvers are required to compute optimal solutions, and both free/open-source and commercial solver can be freely interchanged using interfaces via JuMP (commercial solvers such as the Gurobi LP solver may required for obtaining reasonable run times when solving large-scale problems).

The **DIETER.jl** Julia/JuMP package has been developed at [CSIRO](https://www.csiro.au) within the Energy Systems research program. Thanks to Mario Kendziorski for an earlier version, as well as Wolf-Peter Schill.

In keeping with the original version, **DIETER.jl** is open source and may be freely used and modified by anyone. Contributions are welcome. The code is licensed under the [MIT License](http://opensource.org/licenses/MIT). Input data is licensed under the [Creative Commons Attribution-ShareAlike 4.0 International Public License](http://creativecommons.org/licenses/by-sa/4.0/).

## Model code

For a quick overview, the model code can be directly viewed in the [model.jl](https://github.com/jd-foster/Dieter.jl/blob/dev/src/model.jl) file.

Note that the model is (currently) lacking some features which exist in the GAMS model such as the representation of reserves and prosumage.

General utility functions for handling `DataFrame`s, `Dict`s and translations from JuMP data structures, useful independently of DIETER.jl, are to be found in the [util.jl](https://github.com/jd-foster/Dieter.jl/blob/dev/src/util.jl) file.

## Installation

To install this as a Julia module, download this repository to a sub-folder in, for example, your home folder e.g. 
```julia 
dieter_path=joinpath(ENV["HOME],"Dieter")
```

In the `julia` REPL, run
```julia
import Pkg
Pkg.add(Pkg.PackageSpec(path=dieter_path; rev="dev"))
```
The argument `rev="dev"` selects the git branch to use; change this to use your own alternative branch as required.

## Example of a run script

To run the model, one might use a script similar to the following (_this requires other user modification to work_):
```julia
using Dieter
using JuMP

# select your favourite solver
using Gurobi
solver = Gurobi.Optimizer
# using CPLEX
# solver = CPLEX.Optimizer
# using Clp
# solver = Clp.Optimizer

current_dir = @__DIR__
datapath = joinpath(current_dir,"test","testdata")

# set the Dieter parameters
dtr = InitialiseDieterModel(DieterModel, Dict{String,Any}())
dtr.settings[:datapath] = datapath

initialise_data_file_dict!(dtr)
check_files_exist(dtr.data["files"])

dbfilepath = joinpath(datapath,"DTR_regions_test_data/DTR.db")
dfDict = parse_data_to_model!(dtr;dataname=dbfilepath)

build_model!(dtr)

JuMP.set_optimizer(dtr.model, solver)
solve_model!(dtr)

generate_results!(dtr)
```

