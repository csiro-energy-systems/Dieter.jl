## Running code and tests for Dieter.jl

# %% Load packages
using Dieter
using JuMP
import CPLEX

using LibPQ
using Tables
import DataFrames
using DataFramesMeta
# import SQLite
# import CSV
using TimeSeries
# import XLSX

# using Plots
# using ColorSchemes

# %% In-line functions
getTableAsDF(conn,string) = DataFrame(Tables.columntable(LibPQ.execute(conn,string)))

# %% Year parameters
Demand_Year = 2018 # Financial year 2019-2019

# %% Data paths and connections
trace_read_path = joinpath(ENV["HOME"],"Documents/Projects/NEM_modelling/Data_Tables/REZ_Data/STABLE_input_traces/")
wind_traces_path = joinpath(trace_read_path,"REZ_Wind_Traces_RefYear2019_FYE2030.csv")
solar_traces_path = joinpath(trace_read_path,"REZ_Solar_Traces_RefYear2019_FYE2030.csv")

conn_to_pgdb = "dbname=nemcheck host=/tmp port=5432 user=fos08b"
conn = LibPQ.Connection(conn_to_pgdb)

# %% Instantiate data model
projectpath = joinpath(ENV["HOME"],"Documents/Projects/ESM/Dieter.jl/")
# datapath = joinpath(projectpath,"testdata/")
datapath = joinpath(projectpath,"STABLE_run_data/")
rdir = joinpath(projectpath,"results_stable/")

# Dieter.initialise_data_dir_structure(datapath)

data_instance = Dict{String,Any}()
dtr = InitialiseDieterModel(DieterModel, data_instance)
dtr.settings[:datapath] = datapath
# data_instance["files"] = Dict{String,String}()
# fileDict = dtr.data["files"]

dtr.settings[:interest] = 0.06
# Modify the :min_res setting over [0,100] and rerun to see comparison.
dtr.settings[:min_res] = 10
dtr.settings[:ev] = missing
dtr.settings[:heat] = missing
dtr.settings[:h2] = missing

initialise_data_file_dict!(dtr,"sql")
# dtr.data["files"]
check_files_exist(dtr.data["files"])

sql_db_path = joinpath(datapath,"STABLE_run_data.db")

# %% Load data into the model

# dfDict = parse_data_to_model!(dtr,dataname=sql_db_path)
# dfDict acccesible as dtr.data["dataframes"]
dataname=sql_db_path
datapath = dtr.settings[:datapath]

fileDict = dtr.data["files"]

dtr.data["dataframes"] = Dict{String,DataFrame}()

dfDict = dtr.data["dataframes"]

# %% Base data

import Dieter: parse_file, parse_nodes!, parse_base_technologies!, parse_storages!, parse_load!, parse_availibility!
import Dieter: initialise_set_relation_data!, parse_set_relations!,parse_arcs!, calc_base_parameters!

dfDict["nodes"] = parse_file(fileDict["nodes"]; dataname=dataname)
parse_nodes!(dtr,dfDict["nodes"])

# e.g. fileDict["tech"] = joinpath(datapath,"base","technologies.csv")
dfDict["tech"] = parse_file(fileDict["tech"]; dataname=dataname)

# %% Modify data in-frame:

# Make HydRoR (Hydro Run of River) plants Dispatchable
dt = dfDict["tech"]
@linq dt |> where(:Technologies .== "HydRoR_Exi") |> select(:Dispatchable)
@linq dt |> where(:Technologies .== "HydRoR_Exi") |> select(:Dispatchable)
dt = @byrow! dt if :Technologies == "HydRoR_Exi"; :Dispatchable = 1 end
dt = @byrow! dt if :Technologies == "HydRoR_New"; :Dispatchable = 1 end
@linq dt |> where(:Technologies .== "HydRoR_Exi") |> select(:Dispatchable)
@linq dfDict["tech"] |> where(:Technologies .== "HydRoR_Exi") |> select(:Dispatchable)

# Add a zero Carbon Content column
insertcols!(dt, size(dt)[2], CarbonContent=zeros(size(dt)[1]))

parse_base_technologies!(dtr, dt)
# parse_base_technologies!(dtr, dfDict["tech"])

# e.g. fileDict["storage"] = joinpath(datapath,"base","storages.csv")
dfDict["storage"] = parse_file(fileDict["storage"]; dataname=dataname)
parse_storages!(dtr, dfDict["storage"])

# %% Relations (i.e set-to-set correspondences)
dfDict["map_node_demand"] = parse_file(fileDict["map_node_demand"]; dataname=dataname)
dfDict["map_node_tech"] = parse_file(fileDict["map_node_tech"]; dataname=dataname)
dfDict["map_node_storages"] = parse_file(fileDict["map_node_storages"]; dataname=dataname)
dfDict["arcs"] = parse_file(fileDict["arcs"]; dataname=dataname)
initialise_set_relation_data!(dtr)
parse_set_relations!(dtr)
parse_arcs!(dtr)

# %% Construct demand load data
# Format of columns: TimeIndex, DemandRegion, Load (value)
# e.g. fileDict["load"] = joinpath(datapath,"base","load.csv")
# dfDict["load"] = parse_file(fileDict["load"]; dataname=dataname)

df_OpDem = getTableAsDF(conn,"select * from nem_op_demand where getfinyr(datetime)=$(Demand_Year)")
# The operational demand table has time stamps for keys and one column per region.

# Transformations of the raw data:
# ta_OpDem = TimeArray(df_OpDem, timestamp = :datetime)
dropmissing!(df_OpDem)
insertcols!(df_OpDem, 1, TimeIndex=collect(1:nrow(df_OpDem)))
select!(df_OpDem, Not(:datetime))

reg_rename_dict = Dict(:nsw1 => :NSW1,
                       :qld1 => :QLD1,
                       :sa1 => :SA1,
                       :tas1 => :TAS1,
                       :vic1 => :VIC1)
region_names = sort(collect(values(reg_rename_dict)))
DataFrames.rename!(df_OpDem, reg_rename_dict...)

df_OpDem_stack = stack(df_OpDem, region_names, :TimeIndex)
                        # variable_name=:DemandRegion, value_name=:Load)
dfDict["load"] = @linq df_OpDem_stack |>
                  transform(DemandRegion = String.(:variable), Load = :value) |>
                  select(:TimeIndex,:DemandRegion,:Load)
parse_load!(dtr, dfDict["load"])

# %% Construct availability traces
# Format of columns: TimeIndex, RenewRegionID, TechTypeID, Availability
# e.g. fileDict["avail"] = joinpath(datapath,"base","availability.csv")
# dfDict["avail"] = parse_file(fileDict["avail"]; dataname=dataname)

dfDict["avail"] = DataFrame()

wind_traces = readtimearray(wind_traces_path)
df_wind_traces = DataFrame(wind_traces)
len_traces = nrow(df_wind_traces)
insertcols!(df_wind_traces, 1, TimeIndex=collect(1:len_traces))
select!(df_wind_traces, Not(:timestamp))

# tech_name = "WindOn_Exi"
# if tech_name in ["BattInvWind_Exi", "BattInvWind_New"]
for tech_name in ["WindOn_Exi", "WindOn_New", "WindOff_New"]
      df_trace_mod = copy(df_wind_traces)
      insertcols!(df_trace_mod, 1, TechTypeID=repeat([tech_name],len_traces))
      df_trace_stack = stack(df_trace_mod, Not([:TimeIndex,:TechTypeID]))

       df_trace_select = @linq df_trace_stack |>
                              transform(RenewRegionID = String.(:variable),
                                        Availability = :value) |>
                              select(:TimeIndex,:RenewRegionID,:TechTypeID,:Availability)

      append!(dfDict["avail"],df_trace_select)
end

solar_traces = readtimearray(solar_traces_path)
df_solar_traces = DataFrame(solar_traces)
len_traces = nrow(df_solar_traces)
insertcols!(df_solar_traces, 1, TimeIndex=collect(1:len_traces))
select!(df_solar_traces, Not(:timestamp))

# tech_name = "WindOn_Exi"
# if tech_name in ["BattInvWind_Exi", "BattInvWind_New"]
for tech_name in ["SolFixedPV_Exi", "SolLargePV_Exi", "SolLargePV_New", "SolThermal_New"]
      df_trace_mod = copy(df_solar_traces)
      insertcols!(df_trace_mod, 1, TechTypeID=repeat([tech_name],len_traces))
      df_trace_stack = stack(df_trace_mod, Not([:TimeIndex,:TechTypeID]))

       df_trace_select = @linq df_trace_stack |>
                              transform(RenewRegionID = String.(:variable),
                                        Availability = :value) |>
                              select(:TimeIndex,:RenewRegionID,:TechTypeID,:Availability)

      append!(dfDict["avail"],df_trace_select)
end

parse_availibility!(dtr,dfDict["avail"])

df_REZ_to_Techs = unique(dfDict["avail"][!,[:RenewRegionID,:TechTypeID]])
insertcols!(df_REZ_to_Techs,3,IncludeFlag=ones(Int,nrow(df_REZ_to_Techs)))
rel_node_tech = Dieter.create_relation(df_REZ_to_Techs,:RenewRegionID,:TechTypeID,:IncludeFlag)
Nodes = dtr.sets[:Nodes]
NonDispatchable = dtr.sets[:NonDispatchable]
dtr.sets[:Nodes_NonDispatch] = Dieter.tuple2_filter(rel_node_tech, Nodes, NonDispatchable)
# NND = Dieter.tuple2_filter(rel_node_tech, Nodes, NonDispatchable)

# %% Calculated base parameters
calc_base_parameters!(dtr)

# %% Extensions
Dieter.parse_extensions!(dtr,dataname=sql_db_path)

## Filter tuples examples:
# [x[1] for x in filter(x -> (x[2]=="REZone"),s[:Nodes_Types])]
# [x[1] for x in filter(x -> (x[2]=="TxZone"),s[:Nodes_Types])]

# %% Initialise model

# Construct an optimizer factory
# solver = JuMP.with_optimizer(Clp.Optimizer)
# solver = JuMP.with_optimizer(Gurobi.Optimizer)
solver = JuMP.with_optimizer(CPLEX.Optimizer)
build_model!(dtr,solver)

# %% Specialised solver settings
# CPLEX:
JuMP.set_parameter(dtr.model,"CPX_PARAM_THREADS", 4)        #  number of threads
JuMP.set_parameter(dtr.model, "CPX_PARAM_PARALLELMODE", 1)   #   -1: Opportunistic parallel, 0: Automatic parallel, 1: Deterministic
JuMP.set_parameter(dtr.model, "CPX_PARAM_LPMETHOD", 4)       #  0: auto, 1: primal simplex, 2: dual, 3: network, 4: barrier
JuMP.set_parameter(dtr.model, "CPX_PARAM_BARCROSSALG", 2)    #  0: automatic, 1: primal, 2: dual
JuMP.set_parameter(dtr.model, "CPX_PARAM_SOLUTIONTYPE", 2)   #  Specifies type of solution (basic or non basic) that CPLEX produces
JuMP.set_parameter(dtr.model, "CPX_PARAM_BAREPCOMP", 1e-8)   # Sets the tolerance on complementarity for convergence; default: 1e-8.

# %% Solve the model and generate results
solve_model!(dtr)
generate_results!(dtr)

# %% Analysis
# include("analysis.jl")
df_summ = summarize_result(dtr,del_zeros=false)

# %% Save results to file
# include("save.jl")

# rdir = joinpath(projectpath,resultspath)
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
