## Running code and tests for Dieter.jl

# %% Load packages
using Dieter

import Dieter: parse_file, parse_nodes!, parse_base_technologies!, parse_storages!, parse_load!, parse_availibility!
import Dieter: initialise_set_relation_data!, parse_set_relations!,parse_arcs!, calc_base_parameters!

using JuMP
import MathOptInterface
import CPLEX

# using LibPQ
using Tables
using DataFrames
using DataFramesMeta
import SQLite
import CSV
using TimeSeries
# import XLSX

# %% Misc. utility functions
# Filter a Dictionary for a matching string in the value or key, resp.
dvalmatch(d::Dict,s::Regex) = filter(x -> !ismissing(x.second) && occursin(s,x.second),d)
dkeymatch(d::Dict,s::Regex) = filter(x -> !ismissing(x.second) && occursin(s,x.first),d)
dkeymatch(d::Dict,s::Regex,pos::Int) = filter(x -> !ismissing(x.second) && occursin(s,x.first[pos]),d)
# dkeymatch(d::Dict{Tuple{String,String},Float64},s::Regex,pos::Int) = dkeymatch(d,s,pos)

# %% Scenario Settings (to customise by modeller)

# Year parameters
Demand_Year = 2019 # Financial year 2018-2019
Reference_Year = Demand_Year # 2019  # Year of data set to use for renewable traces
Trace_Year = 2030 # Which year to use from the Reference_Year trace dataset

# Technology
scen_settings =Dict{Symbol,Any}()

scen_settings[:interest] = 0.06
# Modify the :min_res setting over [0,100] and rerun to see comparison.
scen_settings[:min_res] = 10
scen_settings[:ev] = missing
scen_settings[:heat] = missing
scen_settings[:h2] = missing

# Set the scaling of existing coal
scen_settings[:coal_adjust] = 0.5;
# Set the maximum allowed contribution to peak demand by a single technology
scen_settings[:peak_factor] = 0.9
# %% Data paths and connections

# projectpath = joinpath(ENV["HOME"],"Documents/Projects/ESM/Dieter.jl/")
projectpath = pwd()
# datapath = joinpath(projectpath,"testdata/")
datapath = joinpath(projectpath,"STABLE_run_data")
# resultsdir = joinpath(ENV["HOME"],"Documents/Projects/ESM/results_stable")

trace_read_path = joinpath(datapath,"STABLE_input_traces")
wind_traces_path = joinpath(trace_read_path,"REZ_Wind_Traces_RefYear$(Reference_Year)_FYE$(Trace_Year).csv")
solar_traces_path = joinpath(trace_read_path,"REZ_Solar_Traces_RefYear$(Reference_Year)_FYE$(Trace_Year).csv")

# %% Instantiate data model
# Dieter.initialise_data_dir_structure(datapath)

data_instance = Dict{String,Any}()
dtr = InitialiseDieterModel(DieterModel, data_instance)
dtr.settings[:datapath] = datapath
# data_instance["files"] = Dict{String,String}()
# fileDict = dtr.data["files"]

dtr.settings[:interest] = scen_settings[:interest]
# Modify the :min_res setting over [0,100] and rerun to see comparison.
dtr.settings[:min_res] = scen_settings[:min_res]
dtr.settings[:ev] = scen_settings[:ev]
dtr.settings[:heat] = scen_settings[:heat]
dtr.settings[:h2] = scen_settings[:h2]

dtr.settings[:coal_adjust] = scen_settings[:coal_adjust]
dtr.settings[:peak_factor] = scen_settings[:peak_factor]

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

dfDict["nodes"] = parse_file(fileDict["nodes"]; dataname=dataname)
parse_nodes!(dtr,dfDict["nodes"])

# e.g. fileDict["tech"] = joinpath(datapath,"base","technologies.csv")
dfDict["tech"] = parse_file(fileDict["tech"]; dataname=dataname)

# %% Modify data in-frame:

# Make Hydro (Hydro Gravity and Run of River) plants Dispatchable
# dt = dfDict["tech"]
# @linq dt |> where(:Technologies .== "Hydro_Exi") |> select(:Dispatchable)
dfDict["tech"] = @byrow! dfDict["tech"] if :Technologies == "Hydro_Exi"; :Dispatchable = 1 end
dfDict["tech"] = @byrow! dfDict["tech"] if :Technologies == "Hydro_New"; :Dispatchable = 1 end
# @linq dt |> where(:Technologies .== "Hydro_Exi") |> select(:Dispatchable)
@linq dfDict["tech"] |> where(:Technologies .== "Hydro_Exi") |> select(:Dispatchable)

# Filter out new technologies and overwrite the DataFrame
# dfDict["tech"] = @where(dfDict["tech"], occursin.(r"_Exi",:Technologies ))

# Add a Carbon Content column with zero values
insertcols!(dfDict["tech"], size(dfDict["tech"])[2], CarbonContent=zeros(size(dfDict["tech"])[1]))

parse_base_technologies!(dtr, dfDict["tech"])
# parse_base_technologies!(dtr, dfDict["tech"])

# e.g. fileDict["storage"] = joinpath(datapath,"base","storages.csv")
dfDict["storage"] = parse_file(fileDict["storage"]; dataname=dataname)
parse_storages!(dtr, dfDict["storage"])

# %% Relations (i.e set-to-set correspondences)
dfDict["map_node_demand"] = parse_file(fileDict["map_node_demand"]; dataname=dataname)
dfDict["map_node_tech"] = parse_file(fileDict["map_node_tech"]; dataname=dataname)
dfDict["map_node_storages"] = parse_file(fileDict["map_node_storages"]; dataname=dataname)
dfDict["arcs"] = parse_file(fileDict["arcs"]; dataname=dataname)

# %% Modify included technologies

#=
# If the technology is not existing, do not include it.
dfDict["map_node_tech"] =
@byrow! dfDict["map_node_tech"] begin
      if !occursin(r"_Exi",:Technologies)
            :IncludeFlag = 0
      end
end
=#
# @linq dfDict["map_node_tech"] |> where(:IncludeFlag .== 1)
#=
@byrow! dfDict["map_node_storages"] begin
      if !occursin(r"_Exi",:Technologies)
            :IncludeFlag = 0
      end
end
=#
#=
# Get the nodes for existing pumped hydro "HydPump_Exi":
df_phn = @linq dfDict["map_node_storages"] |>
where(:Technologies .== "HydPump_Exi") |>
select(:Nodes)
# Convert to a simpler array for checking against:
hydpump_nodes = df_phn[!,:Nodes]
# Filter out locations without existing pumped hydro, and "BattInvWind_New"
dfDict["map_node_storages"] = @byrow! dfDict["map_node_storages"] begin
      if ( ( :Technologies == "HydPump_New" && !(:Nodes in hydpump_nodes) )  || :Technologies == "BattInvWind_New")
            :IncludeFlag = 0
      end
end
=#


# %% Create relational sets

initialise_set_relation_data!(dtr)
parse_set_relations!(dtr)
parse_arcs!(dtr,dfDict["arcs"])

# %% Construct demand load data
# Format of columns: TimeIndex, DemandRegion, Load (value)
# e.g. fileDict["load"] = joinpath(datapath,"base","load.csv")
# dfDict["load"] = parse_file(fileDict["load"]; dataname=dataname)

load_trace_datapath = joinpath(trace_read_path,"Load_FYE$(Demand_Year).csv")

# In-line function:
# getTableAsDF(conn,string) = DataFrame(Tables.columntable(LibPQ.execute(conn,string)))

# conn_to_pgdb = "dbname=nemcheck host=/tmp port=5432 user=xxx"
# conn = LibPQ.Connection(conn_to_pgdb)
# df_OpDem = getTableAsDF(conn,"select * from nem_op_demand where getfinyr(datetime)=$(Demand_Year)")
# CSV.write(load_trace_datapath,df_OpDem)

# The operational demand table has time stamps for keys and one column per region.
df_OpDem = DataFrame(CSV.File(load_trace_datapath))

# Transformations of the raw data:

# Obtain hourly approximation:
df_OpDem = df_OpDem[1:2:end,names(df_OpDem)]

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

# %%
parse_load!(dtr, dfDict["load"])

# %% Calculate peak demand
DemandRegions = dtr.sets[:DemandRegions]
Peaks = Dict{String,Float64}()
for dr in DemandRegions
      Peaks[dr] = maximum(values(dkeymatch(dtr.parameters[:Load],Regex(dr),1)))
end

dtr.parameters[:Peaks] = Peaks

# %% Construct availability traces
# Format of columns: TimeIndex, RenewRegionID, TechTypeID, Availability
# e.g. fileDict["avail"] = joinpath(datapath,"base","availability.csv")
# dfDict["avail"] = parse_file(fileDict["avail"]; dataname=dataname)

dfDict["avail"] = DataFrame()

# Read data on assignment of trace names to a `Dict`ionary
sqlquery_rez_trace = SQLite.Query(SQLite.DB(sql_db_path), "SELECT * FROM REZ_Trace_Map"; stricttypes=false);
trace_corr = Dieter.SQLqueryToDict(sqlquery_rez_trace)

wind_traces = readtimearray(wind_traces_path)
df_wind_traces = DataFrame(wind_traces)

# Obtain hourly approximation:
df_wind_traces = df_wind_traces[1:2:end,names(df_wind_traces)]

len_traces = nrow(df_wind_traces)
insertcols!(df_wind_traces, 1, TimeIndex=collect(1:len_traces))
select!(df_wind_traces, Not(:timestamp))

# tech_name = "WindOn_Exi"
# if tech_name in ["BattInvWind_Exi", "BattInvWind_New"]
for tech_name in ["WindOn_Exi", "WindOn_New", "WindOff_New"]
      df_trace_mod = copy(df_wind_traces)

      if tech_name == "WindOn_Exi"
            for rz in Symbol.(dtr.sets[:REZones])
                  if rz in names(df_trace_mod) && ismissing(trace_corr[String(rz)]["Wind_Exi"])
                        select!(df_trace_mod, Not(rz))
                  end
            end
      end

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

# Obtain hourly approximation:
df_solar_traces = df_solar_traces[1:2:end,names(df_solar_traces)]

len_traces = nrow(df_solar_traces)
insertcols!(df_solar_traces, 1, TimeIndex=collect(1:len_traces))
select!(df_solar_traces, Not(:timestamp))

# tech_name = "WindOn_Exi"
# if tech_name in ["BattInvWind_Exi", "BattInvWind_New"]
# for tech_name in ["SolFixedPV_Exi", "SolLargePV_Exi", "SolLargePV_New", "SolThermal_New"]
for tech_name in ["SolarPV_Exi", "SolarPV_New", "SolThermal_New"]
      df_trace_mod = copy(df_solar_traces)

      if tech_name in ["SolarPV_Exi"]
            for rz in Symbol.(dtr.sets[:REZones])
                  if rz in names(df_trace_mod) && ismissing(trace_corr[String(rz)]["Solar_Exi"])
                        select!(df_trace_mod, Not(rz))
                  end
            end
      end

      insertcols!(df_trace_mod, 1, TechTypeID=repeat([tech_name],len_traces))
      df_trace_stack = stack(df_trace_mod, Not([:TimeIndex,:TechTypeID]))

       df_trace_select = @linq df_trace_stack |>
                              transform(RenewRegionID = String.(:variable),
                                        Availability = :value) |>
                              select(:TimeIndex,:RenewRegionID,:TechTypeID,:Availability)

      append!(dfDict["avail"],df_trace_select)
end

parse_availibility!(dtr,dfDict["avail"])

## Existing capacity
# sqlquery_exi_cap = SQLite.Query(SQLite.DB(sql_db_path), "SELECT * FROM Existing_Cap"; stricttypes=false);
# exi_cap = Dieter.SQLqueryToDict(sqlquery_exi_cap)

# Set the overnight cost of existing capacity to 0:
OvernightCostPower = dtr.parameters[:OvernightCostPower]
for t in keys(OvernightCostPower)
      if occursin(r"_Exi",t[2])
            OvernightCostPower[t] = 0
      end
end

# Assumption: there is a direct relationship between REZs and NonDispatchable techs.
df_REZ_to_Techs = unique(dfDict["avail"][!,[:RenewRegionID,:TechTypeID]])
insertcols!(df_REZ_to_Techs,3,IncludeFlag=ones(Int,nrow(df_REZ_to_Techs)))
rel_rez_tech = Dieter.create_relation(df_REZ_to_Techs,:RenewRegionID,:TechTypeID,:IncludeFlag)
Nodes = dtr.sets[:Nodes]
NonDispatchable = dtr.sets[:NonDispatchable]
dtr.sets[:Nodes_NonDispatch] = Dieter.tuple2_filter(rel_rez_tech, Nodes, NonDispatchable)
# Reduce set Nodes_NonDispatch by intersecting with
dtr.sets[:Nodes_NonDispatch] = intersect(dtr.sets[:Nodes_NonDispatch],dtr.sets[:Nodes_Techs])
# NND = Dieter.tuple2_filter(rel_rez_tech, Nodes, NonDispatchable)

# %% Calculated base parameters
calc_base_parameters!(dtr)

# %% Extensions
Dieter.parse_extensions!(dtr,dataname=sql_db_path)

# # Filter tuples examples:
# [x[1] for x in filter(x -> (x[2]=="REZone"),s[:Nodes_Types])]
# [x[1] for x in filter(x -> (x[2]=="TxZone"),s[:Nodes_Types])]


# %% Additional parameters

# # Minimum stable generation levels:
fileDict["min_stable_gen"] = joinpath(datapath,"base","min_stable_gen.sql")

dfDict["min_stable_gen"] = parse_file(fileDict["min_stable_gen"]; dataname=dataname)

params_msg = Dieter.map_idcol(dfDict["min_stable_gen"], [:Region, :Technologies], skip_cols=Symbol[])
for (k,v) in params_msg Dieter.update_dict!(dtr.parameters, k, v) end

# # REZ build parameters and bounds
fileDict["rez_build"] = joinpath(datapath,"base","rez_build.sql")
dfDict["rez_build"] = parse_file(fileDict["rez_build"]; dataname=dataname)

params_rzb = Dieter.map_idcol(dfDict["rez_build"], [:Nodes], skip_cols=Symbol[])
for (k,v) in params_rzb Dieter.update_dict!(dtr.parameters, k, v) end

# %% Initialise model

# Construct an optimizer factory
# solver = JuMP.with_optimizer(Clp.Optimizer)
# solver = JuMP.with_optimizer(Gurobi.Optimizer)
# solver = JuMP.with_optimizer(CPLEX.Optimizer)
solver = CPLEX.Optimizer
# build_model!(dtr,solver; nthhour=-1)
build_model!(dtr,solver)

# %% Fix necessary variables for this scenario:

sets = dtr.sets
par = dtr.parameters

N_TECH = dtr.model.obj_dict[:N_TECH]
N_STO_P = dtr.model.obj_dict[:N_STO_P]
N_STO_E = dtr.model.obj_dict[:N_STO_E]
STO_IN = dtr.model.obj_dict[:STO_IN]

MaxEtoP_ratio = dtr.parameters[:MaxEnergyToPowerRatio]

ExistingCapDict = filter(x -> !ismissing(x.second), dtr.parameters[:ExistingCapacity])

Coal_ExistingCapDict = Dict([x for x in par[:ExistingCapacity]
      if x.first in keys(dvalmatch(par[:FuelType],r"Coal")) && par[:Status][x.first] == "GenericExisting"])

# # Fix the capacity of existing generation and storage technologies
for (n,t) in keys(ExistingCapDict)
      if (n,t) in dtr.sets[:Nodes_Techs]
            JuMP.fix(N_TECH[(n,t)], ExistingCapDict[(n,t)]; force=true)
      end
      if (n,t) in dtr.sets[:Nodes_Storages]
            JuMP.fix(N_STO_P[(n,t)], ExistingCapDict[(n,t)] ; force=true)
            JuMP.fix(N_STO_E[(n,t)], MaxEtoP_ratio[n,t]*ExistingCapDict[(n,t)]; force=true)
      end
end

# # (Re)Fix the capacity of existing coal generation at a certain fraction `coal_adjust`
for (n,t) in keys(Coal_ExistingCapDict)
      if (n,t) in dtr.sets[:Nodes_Techs]
            JuMP.fix(N_TECH[(n,t)], dtr.settings[:coal_adjust]*ExistingCapDict[(n,t)]; force=true)
      end
end

# No storage inflow in first period
for (n,sto) in dtr.sets[:Nodes_Storages]
      JuMP.fix(STO_IN[(n,sto),1],0; force=true)
end

df_tech = dfDict["tech"][!,[:Technologies, :Status]] |> dropmissing
tech_status = Dict(zip(df_tech[!,:Technologies],df_tech[!,:Status]))

df_nodes = dfDict["nodes"] |> dropmissing
node_DemReg = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:DemandRegion]))

df_node_tech = @linq dfDict["map_node_tech"] |>
                  where(:IncludeFlag .== 1) |>
                  select(:Nodes, :Technologies) |>
                  dropmissing

df_node_tech[!,:Status] = map(x -> tech_status[x], df_node_tech[!,:Technologies])
df_node_tech[!,:DemandRegion] = map(x -> node_DemReg[x], df_node_tech[!,:Nodes])

# Set an upper bound on each technology as a certain fraction `peak_factor` of the peak demand
for x in eachrow(df_node_tech)
      if x.Status == "NewEntrant"
            # println(copy(x))
            JuMP.set_upper_bound(N_TECH[(x.Nodes,x.Technologies)],
                  dtr.settings[:peak_factor]*Peaks[x.DemandRegion])
      end
end

# NewCapKeys = keys(filter(x -> ismissing(x.second), dtr.parameters[:ExistingCapacity]))


# ExistingStoDict = filter(
#                      x -> (!ismissing(x.second) && (x.first in dtr.sets[:Nodes_Techs])),
#                      dtr.parameters[:ExistingCapacity]
#                    )
# ExistingStoDict = filter(
#                      x -> (!ismissing(x.second) && (x.first in dtr.sets[:Nodes_Storages])),
#                      dtr.parameters[:ExistingCapacity]
#                    )

# NewCapDict = filter(x -> ismissing(x.second), dtr.parameters[:ExistingCapacity])
# for (n,t) in keys(NewCapDict)
#       if t in ["RecipCmpD_New","RecipSpkG_New","OCGTd_New"]
#             # println("Fixing N_TECH[($n, $t)]=0.")
#             JuMP.fix(N_TECH[(n,t)], 0; force=true)
#       end
#
#       if t in ["BattInvWind_New"]
#             # println("Fixing N_STO_P[($n, $t)]=0.")
#             JuMP.fix(N_STO_P[(n,t)], 0; force=true)
#             # println("Fixing N_STO_E[($n, $t)]=0.")
#             JuMP.fix(N_STO_E[(n,t)], 0; force=true)
#       end
# end

# %% Specialised solver settings
# CPLEX:
# JuMP.set_optimizer_attribute(dtr.model,"CPX_PARAM_THREADS", 4)        #  number of threads
JuMP.set_optimizer_attribute(dtr.model, "CPX_PARAM_PARALLELMODE", 0)   #   -1: Opportunistic parallel, 0: Automatic parallel, 1: Deterministic
JuMP.set_optimizer_attribute(dtr.model, "CPX_PARAM_LPMETHOD", 6)       #  0: auto, 1: primal simplex, 2: dual, 3: network, 4: barrier, 6: Concurrent
# JuMP.set_optimizer_attribute(dtr.model, "CPX_PARAM_BARCROSSALG", 2)    #  0: automatic, 1: primal, 2: dual
# JuMP.set_optimizer_attribute(dtr.model, "CPX_PARAM_SOLUTIONTYPE", 2)   #  Specifies type of solution (basic or non basic) that CPLEX produces
JuMP.set_optimizer_attribute(dtr.model, "CPX_PARAM_BAREPCOMP", 1e-6)   # Sets the tolerance on complementarity for convergence; default: 1e-8.

# %% Solve the model and generate results
solve_model!(dtr)
generate_results!(dtr)

# %% Analysis
# # include("analysis.jl")
# df_summ = summarize_result(dtr,del_zeros=false)
#

# %% Save results to file
# # include("save.jl")
#
# # rdir = joinpath(projectpath,resultspath)

# save_results(dtr, resultsdir)
#

# using Serialization

# solved_dtr = copy(dtr)
# solved_dtr.model = []
# Serialization.serialize(joinpath(resultsdir,"dtr_saved.jl"), solved_dtr)

# %% Merge results with other runs
# # include("merge.jl")
# post_process_results(rdir)
#
# %% Filtering results

res = dtr.results

df = res[:G]
df_aug = @byrow! df begin
              @newcol Nodes::Array{String}
              @newcol Techs::Array{String}
              :Nodes = :Nodes_Techs[1]
              :Techs = :Nodes_Techs[2]
       end
# dfDict = dtr.data["dataframes"]
df_nodes = dfDict["nodes"]

node2DemReg = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:DemandRegion]))

df_spatial = join(df_aug,df_nodes,on=:Nodes)
df_filter = select(df_spatial,[:Nodes,:Techs,:DemandRegion,:Hours,:Value])

dfStates = Dict{String,DataFrame}()

for reg in dtr.sets[:DemandRegions]
      dfStates[reg] = @linq df_filter |>
                        where(:DemandRegion .== reg) |>
                        select(:Nodes,:Techs,:Hours,:Value) |>
                        by([:Techs,:Hours], Level = sum(:Value))
end

df_flow = res[:FLOW]
df_flow = @byrow! df_flow begin
            @newcol From::Array{String}
            @newcol To::Array{String}
            :From = :Arcs[1]
            :To   = :Arcs[2]
      end

df_flow[!,:FromRegion] = map(x -> node2DemReg[x], df_flow[!,:From])
df_flow[!,:ToRegion] = map(x -> node2DemReg[x], df_flow[!,:To])

# Inter-state flow:
df_interflow = @where(df_flow, :FromRegion .!== :ToRegion)

# %% Plotting

using Plots
using StatsPlots
import Plots.PlotMeasures: mm
using ColorSchemes

# %% Plotting

Hours = dtr.sets[:Hours]
L = 1:336
# L = 5000:5336

gr()
gr(size=(5000,2000))
# gr(size=(3000,600))
# plotly()
# plotly(size=(3000,600))

# DemandReg = "TAS1"
p = plot(layout=grid(5,1, height=4*[0.1,0.1,0.1,0.1,0.1]),margin=5mm);

for (count, DemandReg) in enumerate(["NSW1", "QLD1", "VIC1", "SA1", "TAS1"])
            # Load = dtr.parameters[:Load]
      Demand = @where(dfDict["load"],:DemandRegion .== DemandReg)
      df_plot = dfStates[DemandReg]
      Techs = [Symbol(i) for i in DataFrames.unique(copy(df_plot[!,:Techs]))]
      NumTechs = length(Techs)
      # reTechs = reshape(Techs,NumTechs,1)

      df_unstack = unstack(df_plot,:Techs,:Level)

      # @df df_all[L,:] groupedbar(HOURS[L], cols(1:NumTech),  #cols(NumTech:-1:1),
      @df df_unstack[L,Techs] groupedbar!(p, Hours[L], cols(Techs),
          subplot=count,
          margin=10mm,
          title=DemandReg,
          xlabel="Time",
          fillalpha=0.5,linealpha=0.1,
          bar_position=:stack,
          legend=:none,  # `:none`, `:best`, `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
          color_palette=:colorwheel) # delta rainbow inferno darkrainbow colorwheel

      plot!(Hours[L], [Demand[L,:Load]],label="Demand",
            subplot=count,
            line=4, linecolour=:steelblue,
            xtickfont = font(10, "Courier"),
            xlabel="Time (hr)",
            ylabel="Generation (MW)",
            )
      # plot!(p,margin=15mm)

      df_regflow = @linq df_interflow |>
                    where(:FromRegion .== DemandReg) |>
                    by(:Hours, Level = sum(:Value))

      plot!(df_regflow[L,:Hours], [df_regflow[L,:Level]],label="Flow",
            subplot=count,
            line=4, linecolour=:red,
            xtickfont = font(10, "Courier"),
            xlabel="Time (hr)",
            ylabel="Generation (MW)",
            margin=5mm
            )
      # display(q)
end
display(p)

# %% Misc.
# color_dict = Dict()
# default_colour = ColorSchemes.leonardo
# # ColorScheme([Colors.RGB(0.0, 0.0, 0.0), Colors.RGB(1.0, 1.0, 1.0)],
#       # "custom", "twotone, black and white");
#       color_dict["Wind Offshore"] = [default_colour[10]];
#       color_dict["Wind Onshore"] = [default_colour[20]];
#       color_dict["Solar PV"] = [default_colour[30]];
#       color_dict["Storages"] = [default_colour[5]];
#       color_dict["Hydrogen"] = [default_colour[15]];
#       color_dict["Curtailment"] = [default_colour[25]];
#
# marker = (:hexagon, 5, 0.6, :green, stroke(3, 0.2, :black, :dot))
#
# plot_all(rdir,color_dict,sector=:ev,marker=marker,legend=true)
