## Running code and tests for Dieter.jl

projectpath = "/Users/fos08b/Documents/Projects/ESM"

datapath = joinpath(projectpath,"Dieter.jl/testdata/")

fileDict = Dict{String,String}()

# %% Base
fileDict["tech"] = joinpath(datapath,"base","technologies.csv")
fileDict["storage"] = joinpath(datapath,"base","storages.csv")
fileDict["load"] = joinpath(datapath,"base","load.csv")
fileDict["avail"] = joinpath(datapath,"base","availability.csv")
# %% EV
fileDict["ev"] = joinpath(datapath,"ev","ev.csv")
fileDict["ev_demand"] = joinpath(datapath,"ev","ev_demand.csv")
fileDict["ev_power"] = joinpath(datapath,"ev","ev_power.csv")

# %% Initialise model
data_instance = Dict{String,Any}()

data_instance["files"] = fileDict

dtr = InitialiseDieterModel(DieterModel, data_instance)

dtr.settings[:datapath] = datapath
