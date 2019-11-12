To run the model you can use this script:
```julia
using Dieter
# select your favourite solver
using Gurobi
#using CPLEX #using Clp

datapath = joinpath(pwd(),"testdata")

# set the Dieter parameters
resshare = 100
ev = 0
heat = 0
h2 = 0

dieter = DieterModel(datapath, ev=ev, heat=heat, h2=h2, res_share=resshare)
solve_model!(dieter, with_optimizer(Gurobi.Optimizer))

save_results(dieter, "results_of_test_case")
```
