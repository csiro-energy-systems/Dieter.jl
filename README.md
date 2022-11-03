# Dieter.jl

[![In Development](https://img.shields.io/badge/docs-main-blue.svg)](https://csiro-energy-systems.github.io/Dieter.jl/main/)

[DIETER](http://www.diw.de/dieter) is the _Dispatch and Investment Evaluation Tool with Endogenous Renewables_, an open energy system model originally developed by the German Institute for Economic Research ([DIW Berlin](https://www.diw.de/en)). The model is structured as a cost-minimising [linear program](https://en.wikipedia.org/wiki/Linear_programming) that when solved to optimality will produce an hourly operation profile of an interconnected transmission system over a single year, including aspects such as
- generator power output,
- storage technology charge and discharge,
- transmission flows,
- variable renewable energy utilisation, and
- hydrogen electrolyser operation.

The model also decides on cost-optimal capacity expansion of relevant model technologies relative to the year under consideration.

This code repository holds a derivative version of DIETER written using the [Julia](https://julialang.org/) language and the [JuMP](https://github.com/jump-dev/JuMP.jl) algebraic modelling language (a Julia package). While the original version of DIETER is implemented in [GAMS](https://www.gams.com/), Julia/JuMP enables a free and open-source implementation. Linear programming (LP) solvers are required to compute optimal solutions, and both free/open-source and commercial solver can be freely interchanged using interfaces via JuMP (commercial solvers such as the Gurobi LP solver may required for obtaining reasonable run times when solving large-scale problems).

CSIRO uses this version of the DIETER model in the STABLE model of Australia's National Electricity Market (NEM) in the east and southern states, and the South West Interconnected System (SWIS) and North West Interconnected System (NWIS) in Western Australia. This package is the core model for STABLE but is not supplied with the required data to run the STABLE model or any other real-world system.

## Documentation

The latest documentation can be found [here](https://csiro-energy-systems.github.io/Dieter.jl/main/).

## Acknowledgements

The *Dieter.jl* Julia/JuMP package has been developed at [CSIRO](https://www.csiro.au) within the Energy Systems research program.

Thanks to Mario Kendziorski for an earlier initial version, as well as Wolf-Peter Schill of DIW Berlin.

The CSIRO maintainer of this package is James Foster.


## License

See [LICENSE](LICENSE) for this software's license terms. If using this software, you must abide by this license agreement. Contributions and issues are welcome but no guarantees of support are made.

In keeping with the original version, *Dieter.jl* is open source and may be freely used and modified by anyone. The code is licensed under the CSIRO Open Source Software Licence, a variation of the BSD / [MIT License](http://opensource.org/licenses/MIT). Any input data stored in this repository is licensed under the [Creative Commons Attribution-ShareAlike 4.0 International Public License](http://creativecommons.org/licenses/by-sa/4.0/).

## Model code

For a quick overview, the model code can be directly viewed in the [model.jl](https://github.com/csiro-energy-systems/Dieter.jl/blob/main/src/model.jl) file.

Note that this model differs from the original GAMS model in some respects, such as no direct representation of reserves and prosumage, while some new features exist, such as hydrogen systems and more detailed transmission constraints.

General utility functions for handling `DataFrame`s, `Dict`s and translations from JuMP data structures, useful independently of *Dieter.jl*, are to be found in the [util.jl](https://github.com/csiro-energy-systems/Dieter.jl/blob/main/src/util.jl) file.

## Installation

This package is not registered on the Julia package system. To install this as a Julia module in a Julia REPL session, you can run the following package installation commands to install directly from the model code repository:
```julia
import Pkg
Pkg.add(Pkg.PackageSpec(url="https://github.com/csiro-energy-systems/Dieter.jl.git"; rev="main"))
```
The argument `rev="main"` selects the git branch to use; change this to use your own alternative branch as required.


