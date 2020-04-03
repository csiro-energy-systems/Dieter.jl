## Running code and tests for Dieter.jl

# %% Load packages
using Dieter
import Dieter: parse_file, parse_nodes!, parse_base_technologies!, parse_storages!, parse_load!, parse_availibility!
import Dieter: initialise_set_relation_data!, parse_set_relations!,parse_arcs!, calc_base_parameters!, parse_extensions!
import Dieter: dvalmatch, dkeymatch, split_df_tuple

using JuMP
import MathOptInterface
import CPLEX

using DataFrames
using DataFramesMeta
import DBInterface
import SQLite
import CSV
import TimeSeries
using Dates
# using Tables
# import XLSX

# %% Scenario Settings (to customise by modeller)

run_timestamp = "$(Date(Dates.now()))-H$(hour(now()))"

# ScenarioName = "Scen1_BAU"
ScenarioName = "Scen2_DDC"
Scen_Map = Dict("Scen1_BAU" => "4deg", "Scen2_DDC" => "2deg")
# Specfied Year for the scenario setting:
ScenarioYear = 2030
ScYr_Sym = Symbol("FYE$ScenarioYear")

NoNewGas = false
Note = "Testing"

BattEnergyType = "N_BattEnergy"
HydPumpEnergyType = "N_HydPumpEnergy"

# Year parameters
Demand_Year = 2019 # Financial year 2018-2019
Reference_Year = Demand_Year # 2019  # Year of data set to use for renewable traces
Trace_Year = 2030 # Which year to use from the Reference_Year trace dataset

# Technology
scen_settings =Dict{Symbol,Any}()

scen_settings[:scen] = run_timestamp*"-$(ScenarioName)-ScYr$(ScenarioYear)-$(Note)"
scen_settings[:interest] = 0.06
scen_settings[:cost_scaling] = 1 # 1.0e-6
# Modify the :min_res setting over [0,100] and rerun to see comparison.
scen_settings[:min_res] = 10
scen_settings[:ev] = missing
scen_settings[:heat] = missing
scen_settings[:h2] = 0  # missing -> H2 not included, any number -> H2 included


# Half-hourly: timestep=1, Hourly: timestep=2,
# If timestep = 2, we obtain an hourly approximation of half-hour data by sampling on every second data point.
scen_settings[:timestep] = 2
timestep = scen_settings[:timestep]

# Set the scaling of existing coal
scen_settings[:coal_adjust] = 1;
# Set the maximum allowed contribution to peak demand by a single technology
# scen_settings[:peak_factor] = 2.5;

# CarbonBudgetDict by ScenarioName and ScenarioYear in millions of t-CO2 (Mt-CO2e)
# Construction should normally be automated from a carbon scenario table...
# CarbonBudgetDict["Scen1_BAU",2030] = 3500
# CarbonBudgetDict["Scen1_BAU",2050] = 2800
# CarbonBudgetDict["Scen2_DDC",2030] = 1500
# CarbonBudgetDict["Scen2_DDC",2050] = 0
#
# scen_settings[:carbon_budget] = 1000*CarbonBudgetDict[ScenarioName,ScenarioYear]


# %% Data paths and connections (to customise by modeller)

# projectpath = joinpath(ENV["HOME"],"Documents/Projects/ESM/Dieter.jl/")
projectpath = pwd()
# datapath = joinpath(projectpath,"testdata/")
datapath_STABLE = joinpath(projectpath,"STABLE_run_data")

if Base.Sys.isapple()
      resultsdir = joinpath(ENV["HOME"],"Documents/Projects/ESM/","results_STABLE")
elseif Base.Sys.iswindows()
      resultsdir = joinpath("F:\\STABLE\\","results_STABLE")
end

results_filename = scen_settings[:scen]*"-results-Julia_Serial.dat"

sql_db_path = joinpath(datapath_STABLE,"STABLE_run_data.db")
dataname=sql_db_path

trace_read_path = joinpath(datapath_STABLE,"STABLE_input_traces")

# %% Instantiate data model
# # Initialise a basic file structure for accessing the model data:
# Dieter.initialise_data_dir_structure(datapath)

data_instance = Dict{String,Any}()
dtr = InitialiseDieterModel(DieterModel, data_instance; datapath=datapath_STABLE, settings=scen_settings)

datapath = dtr.settings[:datapath]

dtr.data["dataframes"] = Dict{String,DataFrame}()
dfDict = dtr.data["dataframes"]

# Create dtr.data["files"] Dict:
initialise_data_file_dict!(dtr,"sql")

check_files_exist(dtr.data["files"])
fileDict = dtr.data["files"]

# # Generically we run this to instantiate the data:
# Dieter.parse_data_to_model!(dtr; dataname=dataname)
# # else:
# %% Load data into the model: Base data

dfDict["nodes"] = parse_file(fileDict["nodes"]; dataname=dataname)
parse_nodes!(dtr,dfDict["nodes"])

df_nodes = dfDict["nodes"]
node2DemReg = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:DemandRegion]))

dtr.parameters[:node_demreg_map] = node2DemReg

# e.g. fileDict["tech"] = joinpath(datapath,"base","technologies.csv")
dfDict["tech"] = parse_file(fileDict["tech"]; dataname=dataname)

# e.g. fileDict["storage"] = joinpath(datapath,"base","storages.csv")
dfDict["storage"] = parse_file(fileDict["storage"]; dataname=dataname)

# %% Relations (i.e set-to-set correspondences)
dfDict["map_node_demand"] = parse_file(fileDict["map_node_demand"]; dataname=dataname)
dfDict["map_node_tech"] = parse_file(fileDict["map_node_tech"]; dataname=dataname)
dfDict["map_node_storages"] = parse_file(fileDict["map_node_storages"]; dataname=dataname)
dfDict["arcs"] = parse_file(fileDict["arcs"]; dataname=dataname)

# %% Renewable Energy Targets

# Load in the complete RET scenario data sets:
fileDict["re_targets"] = joinpath(datapath,"base","re_targets.sql")
dfDict["re_targets"] = parse_file(fileDict["re_targets"]; dataname=dataname)

# Filter the data set for the current Scenario, and put in the model settings:
dtr.settings[:min_res] = Dict(eachrow(
      @linq dfDict["re_targets"] |>
      where(:ScenarioName .== ScenarioName, :ScenarioYear .== ScenarioYear) |>
      select(:Region, :MinRET)
      ))


# %% Additional parameters

# Remove new gas if specified:
if NoNewGas
      GasTech_New = @linq dfDict["tech"] |>
                  where(:Status .== "NewEntrant", occursin.(r"Gas",:FuelType)) |>
                  select(:Technologies)
      GasTech_New = (unique(GasTech_New,:Technologies))[!,:Technologies]

      dfDict["map_node_tech"] = @byrow! dfDict["map_node_tech"] begin
            if in(:Technologies,GasTech_New)
                  :IncludeFlag = 0
            end
      end
end

# # New capital / overnight costs
fileDict["capital_costs"] = joinpath(datapath,"base","capital_costs.sql")
dfDict["capital_costs"] = parse_file(fileDict["capital_costs"]; dataname=dataname)

# Rename the ScenarioYear column with ScenCostPower:
df_cap_costs = DataFrames.rename!(
      dfDict["capital_costs"][!,[:Technologies, :Scenario, ScYr_Sym]],
      Dict(ScYr_Sym => :ScenCostPower)
      )
# Filter by the current Scenario:
df_cap_costs = @linq df_cap_costs |> where(:Scenario .== Scen_Map[ScenarioName])

params_cap_costs = Dieter.map_idcol(df_cap_costs, [:Technologies], skip_cols=Symbol[:Scenario])

for (k,v) in params_cap_costs Dieter.update_dict!(dtr.parameters, k, v) end

fileDict["cost_energy_storage"] = joinpath(datapath,"base","cost_energy_storage.sql")
dfDict["cost_energy_storage"] = parse_file(fileDict["cost_energy_storage"]; dataname=dataname)

# Rename the ScenarioYear column with ScenCostEnergy:
df_costES = DataFrames.rename!(
      dfDict["cost_energy_storage"][!,[:Technologies, :Scenario, ScYr_Sym]],
      Dict(ScYr_Sym => :ScenCostEnergy)
      )

df_costES = @linq df_costES |> where(:Scenario .== Scen_Map[ScenarioName])

params_costES = Dieter.map_idcol(df_costES, [:Technologies], skip_cols=Symbol[:Scenario])
for (k,v) in params_costES Dieter.update_dict!(dtr.parameters, k, v) end

# %% Synchronous condensers
# Synchronous condenser with flywheel - cost in $ per inertial units of `MWs`

SynConOvernightCost =  round(1.0E6*185.2/4400, digits=2) # Units:  $/MWs
# ElectraNet project of $185.2 million, for 4400 MWs of inertia via 4 sync. cons.
# ElectraNet, "Main grid system strength project: Contingent project application", p. 21, 28 June 2019, p. 19.
# All dollar amounts in this document are in real, $2017–18 in line with the ElectraNet's revenue determination unless otherwise stated.

LifetimeSynCon = 40   # Source: GHD; Economic life for ElectraNet synchronous condensers - ElectraNet 28 June 2019
i = dtr.settings[:interest]
SynConCapCost = SynConOvernightCost*Dieter.annuity(i, LifetimeSynCon)
dtr.parameters[:SynConCapCost] = Dict("SynConNew" => SynConCapCost)
# SynCon_New_Cost =

# %% Additional parameters
# # Minimum stable generation levels:
fileDict["min_stable_gen"] = joinpath(datapath,"base","min_stable_gen.sql")
dfDict["min_stable_gen"] = parse_file(fileDict["min_stable_gen"]; dataname=dataname)

params_msg = Dieter.map_idcol(dfDict["min_stable_gen"], [:Region, :Technologies], skip_cols=Symbol[])
for (k,v) in params_msg Dieter.update_dict!(dtr.parameters, k, v) end

# # REZ build parameters and bounds
fileDict["rez_build"] = joinpath(datapath,"base","rez_build.sql")
dfDict["rez_build"] = parse_file(fileDict["rez_build"]; dataname=dataname)

params_rzb = Dieter.map_idcol(dfDict["rez_build"], [:Region], skip_cols=Symbol[])
for (k,v) in params_rzb Dieter.update_dict!(dtr.parameters, k, v) end

# # RE Zone connections parameters and bounds
fileDict["rez_connect"] = joinpath(datapath,"base","rez_connect.sql")
dfDict["rez_connect"] = parse_file(fileDict["rez_connect"]; dataname=dataname)

params_rzc = Dieter.map_idcol(dfDict["rez_connect"], [:Region, :Technologies], skip_cols=Symbol[])
for (k,v) in params_rzc Dieter.update_dict!(dtr.parameters, k, v) end

# # Transmission Zone connections parameters and bounds
fileDict["TxZ_connect"] = joinpath(datapath,"base","TxZ_connect.sql")
dfDict["TxZ_connect"] = parse_file(fileDict["TxZ_connect"]; dataname=dataname)

params_txc = Dieter.map_idcol(dfDict["TxZ_connect"], [:Region, :Technologies], skip_cols=Symbol[])
for (k,v) in params_txc Dieter.update_dict!(dtr.parameters, k, v) end

# %% Technology capacity overwrites by scenario:
# This creates a parameter called ScenarioCapacity indexed by Region, Technologies
fileDict["tech_scenario"] = joinpath(datapath,"base","tech_scenario.sql")
dfDict["tech_scenario"] = parse_file(fileDict["tech_scenario"]; dataname=dataname)

df_techscen = DataFrames.rename!(
      dfDict["tech_scenario"][!,[:Region, :TechID, ScYr_Sym]],
      Dict(:TechID => :Technologies, ScYr_Sym => :ScenarioCapacity)
      )

for (k,v) in Dieter.map_idcol(
            df_techscen,
            [:Region, :Technologies],
             skip_cols=Symbol[]
            )
      Dieter.update_dict!(dtr.parameters, k, v)
end

ScenarioCapacityDict = dtr.parameters[:ScenarioCapacity]
# OverwriteCapDict = Dict([(n,t) => ScenarioCapacityDict[n,t,y,sc]
#                         for (n,t,y,sc) in keys(ScenCapacityDict)
#                         if (y == ScenarioYear && sc == ScenarioName)])

# This will be used later to overwrite and fix certain capacity in the model.
# It may be necessary to check technologies included here;
# particularly certain fossil fuel tech. types like OCGT w. diesel/distillate.

# %% Carbon

fileDict["carbon_param"] = joinpath(datapath,"base","carbon_param.sql")
dfDict["carbon_param"] = parse_file(fileDict["carbon_param"]; dataname=dataname)

Scen_co2 = @where(dfDict["carbon_param"], :ScenarioName .== ScenarioName, :ScenarioYear .== ScenarioYear)
# Scen_co2[!, :CarbonPrice][1]
dtr.settings[:co2] = Scen_co2[!,:CarbonPrice][1]
dtr.settings[:carbon_budget] = Scen_co2[!,:CarbonBudget][1]

# Add a Carbon Content column with zero values
if !in(:CarbonContent, names(dfDict["tech"]))
      insertcols!(dfDict["tech"], size(dfDict["tech"])[2], CarbonContent=zeros(size(dfDict["tech"])[1]))
end

# Note: CarbonContent is generally assumed to be in Units t-CO2/MWh-thermal
fileDict["carbon_content"] = joinpath(datapath,"base","carbon_content.sql")
dfDict["carbon_content"] = parse_file(fileDict["carbon_content"]; dataname=dataname)

fueltype_cc_dict = Dieter.map_idcol(dfDict["carbon_content"], [:FuelType], skip_cols=Symbol[])
fuel_to_cc_dict = fueltype_cc_dict[:CarbonContent]
dfDict["tech"] = @byrow! dfDict["tech"] begin
            :CarbonContent = fuel_to_cc_dict[:FuelType]
      end

# dtr.parameters[:CarbonBudget] = scen_settings[:carbon_budget]

# %% Modify data in-frame:

# Relation to determine which (Nodes,Technologies) pairs are included in the model:
rel_node_tech = create_relation(dfDict["map_node_tech"],:Nodes,:Technologies,:IncludeFlag)
rel_node_storages = create_relation(dfDict["map_node_storages"],:Nodes,:Technologies,:IncludeFlag)

# Existing capacity
# sqlquery_exi_cap = SQLite.Query(SQLite.DB(sql_db_path), "SELECT * FROM Existing_Cap"; stricttypes=false);
# exi_cap = Dieter.SQLqueryToDict(sqlquery_exi_cap)

# Set the overnight cost of existing capacity to 0:
dfDict["tech"] = @byrow! dfDict["tech"] if :Status == "GenericExisting"; :OvernightCostPower = 0 end
dfDict["storage"] = @byrow! dfDict["storage"] if :Status == "GenericExisting"; :OvernightCostPower = 0 end
dfDict["storage"] = @byrow! dfDict["storage"] if :Status == "GenericExisting"; :OvernightCostEnergy = 0 end

# Set the overnight cost of new capacity based on scenario:
ScenCostPower = dtr.parameters[:ScenCostPower]
dfDict["tech"] = @byrow! dfDict["tech"] begin
       if :Status == "NewEntrant" && rel_node_tech(:Region,:Technologies) == true
             :OvernightCostPower = ScenCostPower[:Technologies]
       end
       # if :Status == "NewEntrant" && rel_node_tech(:Region,:Technologies) == false
       #       :OvernightCostPower = missing
       # end
 end

BattEnergyCost = dtr.parameters[:ScenCostEnergy][BattEnergyType]
HydPumpEnergyCost = dtr.parameters[:ScenCostEnergy][HydPumpEnergyType]
dfDict["storage"] = @byrow! dfDict["storage"] begin
       if :Status == "NewEntrant" && rel_node_storages(:Region,:Storages) == true
            if occursin(r"Battery",:TechType)
                 :OvernightCostEnergy = BattEnergyCost
           elseif occursin(r"Pump Storage",:TechType)
                 :OvernightCostEnergy = HydPumpEnergyCost
           end
       end
 end

# # Hornsdale
dfDict["storage"] = @byrow! dfDict["storage"] begin
      if :Storages == "BattInvWind_Exi" && :Region == "NSA"
             :MaxEnergyToPowerRatio = 1.0
      end
end

# # Shoalhaven (NSW)
dfDict["storage"] = @byrow! dfDict["storage"] begin
      if :Storages == "HydPump_Exi" && :Region == "NCEN"
             :MaxEnergyToPowerRatio = 64.0
      end
end

# # Snowy 2.0
Snowy2_Capacity = 2000.0 # MW
Snowy2_HoursStorage = 168.0
dfDict["storage"] = @byrow! dfDict["storage"] begin
      if :Storages == "HydPump_Exi" && :Region == "SWNSW"
             :ExistingCapacity = :ExistingCapacity + Snowy2_Capacity
             :MaxEnergyToPowerRatio = Snowy2_HoursStorage
      end
end
@linq dfDict["storage"] |> where(:Storages .== "HydPump_Exi") |> select(:Region, :ExistingCapacity, :MaxEnergyToPowerRatio)

# # Battery of the Nation...

ConnectCost = dtr.parameters[:ConnectCost]

# Add connection costs specific to regions to the capital cost of plants:
dfDict["tech"] = @byrow! dfDict["tech"] begin
      # @newcol IncCost::Array{Float64}
      if :Status == "NewEntrant" && rel_node_tech(:Region,:Technologies) == true
            :OvernightCostPower = :OvernightCostPower + ConnectCost[(:Region,:Technologies)]
      # else
            # :IncCost = :OvernightCostPower
      end
end

# Make Hydro (Hydro Gravity and Run of River) plants Dispatchable
dfDict["tech"] = @byrow! dfDict["tech"] if :TechType .== "Hydro"; :Dispatchable = 1 end
# dfDict["tech"] = @byrow! dfDict["tech"] if :Technologies == "Hydro_Exi"; :Dispatchable = 1 end
# dfDict["tech"] = @byrow! dfDict["tech"] if :Technologies == "Hydro_New"; :Dispatchable = 1 end

@linq dfDict["tech"] |> where(:TechType .== "Hydro") |> select(:Dispatchable)

# Turn a missing ExisitingCapacity value into a 0 value:
# alt: @linq dfDict["tech"] |> where(:Status .== NewEntrant)
dfDict["tech"] = @byrow! dfDict["tech"] if :ExistingCapacity |> ismissing; :ExistingCapacity = 0 end
# dropmissing!(dfDict["tech"],:ExistingCapacity)
dfDict["storage"] = @byrow! dfDict["storage"] if :ExistingCapacity |> ismissing; :ExistingCapacity = 0 end
# dropmissing!(dfDict["storage"],:ExistingCapacity)

dfExistingCap = @linq dfDict["tech"] |>
                  select(:Region, :Technologies, :TechType, :ExistingCapacity) |>
                  where(:ExistingCapacity .> 0)

# Modify (if necessary) upper bound on REZones that incorporates existing capacity:
TotalBuildCap = dtr.parameters[:TotalBuildCap]
for rez in keys(TotalBuildCap)
      existing_rez_capacity = @where(dfExistingCap, :Region .== rez, :TechType .!== "Hydro")
      rez_sum = sum(existing_rez_capacity[!,:ExistingCapacity])
      new_rez_bound = max(TotalBuildCap[rez],rez_sum)
      if TotalBuildCap[rez] !== new_rez_bound
            @info "TotalBuildCap parameter overwritten for REZone $(rez)."
      end
      TotalBuildCap[rez] = new_rez_bound
end

# %% Parse data into model structure:

parse_base_technologies!(dtr, dfDict["tech"])
# parse_base_technologies!(dtr, dfDict["tech"])

parse_storages!(dtr, dfDict["storage"])

# %% Create relational sets

initialise_set_relation_data!(dtr)
parse_set_relations!(dtr)
parse_arcs!(dtr,dfDict["arcs"])

# %% Construct demand load data
# Format of columns: TimeIndex, DemandRegion, Load (value)
# e.g. fileDict["load"] = joinpath(datapath,"base","load.csv")
# dfDict["load"] = parse_file(fileDict["load"]; dataname=dataname)

load_trace_datapath = joinpath(trace_read_path,"Load_FYE$(Demand_Year).csv")

# The operational demand table has time stamps for keys and one column per region.
df_OpDem = DataFrame(CSV.File(load_trace_datapath))

# Transformations of the raw data:

# If timestep = 2, we obtain an hourly approximation of half-hour data by sampling on every second data point:
df_OpDem = df_OpDem[1:timestep:end,names(df_OpDem)]

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
# CSV.write(joinpath(resultsdir,"$(run_timestamp)-Load.csv"),dfDict["load"])
parse_load!(dtr, dfDict["load"])

# %% Calculate peak demand
DemandRegions = dtr.sets[:DemandRegions]
Peaks = Dict{String,Float64}()
for dr in DemandRegions
      Peaks[dr] = maximum(values(dkeymatch(dtr.parameters[:Load],Regex(dr),1)))
end

dtr.parameters[:Peaks] = Peaks

# %% Inertia -  data and constraints

fileDict["inertia_tech"] = joinpath(datapath,"base","inertia_tech.sql")
dfDict["inertia_tech"] = parse_file(fileDict["inertia_tech"]; dataname=dataname)

for (k,v) in Dieter.map_idcol(
            dfDict["inertia_tech"], [:Technologies], skip_cols=Symbol[])
                  Dieter.update_dict!(dtr.parameters, k, v)
end

fileDict["inertia_require"] = joinpath(datapath,"base","inertia_require.sql")
dfDict["inertia_require"] = parse_file(fileDict["inertia_require"]; dataname=dataname)

for (k,v) in Dieter.map_idcol(
            dfDict["inertia_require"], [:Region], skip_cols=Symbol[])
                  Dieter.update_dict!(dtr.parameters, k, v)
end

# FYE2018 reference demand peaks:
Peaks_BaseYear = Dict("QLD1" => 9383, "NSW1" => 14226,"VIC1" => 9886, "SA1" => 3181, "TAS1" => 1389)

RequireRatio = Dict([dr => Peaks[dr]/Peaks_BaseYear[dr] for dr in DemandRegions])
dtr.parameters[:RequireRatio] = RequireRatio
# SynCon_New_Cost = 37727.0  # $/MWs

# Nodes_Techs = dtr.sets[:Nodes_Techs]
# Nodes_Storages = dtr.sets[:Nodes_Storages]
# Nodes_Dispatch = dtr.sets[:Nodes_Dispatch]
# dtr.parameters[:InertialCoeff] = InertialSecs
# dtr.parameters[:InertiaMinThreshold] = InertiaMinThreshold
# dtr.parameters[:InertiaMinSecure] = InertiaMinSecure

# Define a synchronous condenser capacity for each demand region:
# @variable(dtr.model, N_SYNC[DemandRegions], lower_bound=0)

# for dr in DemandRegions
#       conref = constraint_by_name(dtr.model,"InertiaNormalThreshold[$(dr)]")
#       if is_valid(dtr.model,conref)
#             println(dr)
#             # delete(dtr.model, constraint_by_name(dtr.model,"InertiaNormalThreshold[$(dr)]"))
#       end
# end



# %% Construct availability traces
# Format of columns: TimeIndex, RenewRegionID, TechTypeID, Availability
# e.g. fileDict["avail"] = joinpath(datapath,"base","availability.csv")
# dfDict["avail"] = parse_file(fileDict["avail"]; dataname=dataname)

dfDict["avail"] = DataFrame()

# Read data to a `Dict`ionary that creates an assignment of trace names:
sqlquery_rez_trace = DBInterface.execute(SQLite.DB(sql_db_path), "SELECT * FROM REZ_Trace_Map") #; stricttypes=false
trace_corr = Dieter.SQLqueryToDict(sqlquery_rez_trace)

wind_traces_path = joinpath(trace_read_path,"REZ_Wind_Traces_RefYear$(Reference_Year)_FYE$(Trace_Year).csv")
solar_traces_path = joinpath(trace_read_path,"REZ_Solar_Traces_RefYear$(Reference_Year)_FYE$(Trace_Year).csv")

dtr.data["files"]["Wind"] = wind_traces_path
dtr.data["files"]["Solar"] = solar_traces_path

tech_subset = Dict()
tech_subset["Wind"]  = ["WindOn_Exi", "WindOn_New", "WindOff_New"]
tech_subset["Solar"] = ["SolarPV_Exi", "SolarPV_New", "SolThermal_New"]
# alt.: ["SolFixedPV_Exi", "SolLargePV_Exi", "SolLargePV_New", "SolThermal_New"]

for tech_class in ["Wind", "Solar"]
      traces_path = dtr.data["files"][tech_class]
      # Read trace data:
      df_traces = DataFrame(TimeSeries.readtimearray(traces_path))

      # If timestep = 2, we obtain an hourly approximation of half-hour data by sampling on every second data point:
      df_traces = df_traces[1:timestep:end,names(df_traces)]

      # Construct an integer OneTo(n)-like time indexing instead of timestamps:
      len_traces = nrow(df_traces)
      insertcols!(df_traces, 1, TimeIndex=collect(1:len_traces))
      select!(df_traces, Not(:timestamp))

      for tech_name in tech_subset[tech_class]
            df_trace_mod = copy(df_traces)

            nodes_for_tech = [x[1] for x in dtr.sets[:Nodes_Techs] if x[2] == tech_name] |> sort #  println(nodes_for_tech)

            # Require all necessary traces to exist:
            for rz in dtr.sets[:REZones]
                  if rz in nodes_for_tech && !(Symbol(rz) in names(df_trace_mod))
                        error("Require trace in $(rz) for technology $(tech_name)")
                  end
            # Remove unnecessary traces from data before compiling into model:
                  if Symbol(rz) in names(df_trace_mod) && !(rz in nodes_for_tech)
                        @debug "Tech. not in zone: Removing $(tech_name) trace in $(rz)"
                        select!(df_trace_mod, Not(Symbol(rz)))
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
end

parse_availibility!(dtr,dfDict["avail"])

# # Assumption: there is a direct relationship between REZs and NonDispatchable techs.
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

parse_extensions!(dtr,dataname=sql_db_path)

# # Hydrogen - code for testing sub-functions
# fileDict["h2_technologies"] = joinpath(datapath,"h2","h2_technologies.sql")
# dfDict["h2_technologies"] = parse_file(fileDict["h2_technologies"]; dataname=dataname)
#
# Dieter.parse_h2_technologies!(dtr, dfDict["h2_technologies"])
#
# rel_node_h2tech = create_relation(dfDict["h2_technologies"],:Region,:H2Technologies,:Efficiency)
# Nodes = dtr.sets[:Nodes]
# H2Technologies = dtr.sets[:H2Technologies]
# dtr.sets[:Nodes_H2Tech] = Dieter.tuple2_filter(rel_node_h2tech, Nodes, H2Technologies)
# Dieter.calc_inv_gas!(dtr)
#
# temp_h2_set = dtr.settings[:h2]
# dtr.settings[:h2] = missing
# Dieter.parse_extensions!(dtr,dataname=sql_db_path)
# dtr.settings[:h2] = temp_h2_set

# %% Hydrogen parameters

# dfDict["h2_technologies"] = parse_file(fileDict["h2_technologies"]; dataname=dataname)
# Dieter.parse_h2_technologies!(dtr, dfDict["h2_technologies"])


# %% Initialise model

# Construct an optimizer factory
# solver = JuMP.with_optimizer(Clp.Optimizer)
# solver = JuMP.with_optimizer(Gurobi.Optimizer)
# solver = JuMP.with_optimizer(CPLEX.Optimizer)
solver = CPLEX.Optimizer
# build_model!(dtr,solver; timestep=-1)
build_model!(dtr,solver,timestep=timestep)

# %% Access model objects for further development:

sets = dtr.sets
par = dtr.parameters

G = dtr.model.obj_dict[:G]
N_TECH = dtr.model.obj_dict[:N_TECH]
N_STO_P = dtr.model.obj_dict[:N_STO_P]
N_STO_E = dtr.model.obj_dict[:N_STO_E]
STO_IN = dtr.model.obj_dict[:STO_IN]
N_RES_EXP = dtr.model.obj_dict[:N_RES_EXP]
N_SYNC = dtr.model.obj_dict[:N_SYNC]

# %% Fix necessary variables for this scenario:

# No storage inflow in first period
for (n,sto) in dtr.sets[:Nodes_Storages]
    JuMP.fix(STO_IN[(n,sto),1],0; force=true)
end


NoExpansionREZones = keys(filter(x -> ismissing(x[2]), dtr.parameters[:TransExpansionCost]))
for rez in NoExpansionREZones
      JuMP.fix(N_RES_EXP[rez], 0; force=true)
end

MaxEtoP_ratio = dtr.parameters[:MaxEnergyToPowerRatio]

ExistingCapDict = filter(x -> x.second > 0, dtr.parameters[:ExistingCapacity])
# ExistingCapDict = filter(x -> !ismissing(x.second), dtr.parameters[:ExistingCapacity])

# Overwrite with any imposed scenario capacities read via "tech_scenario" files:
FixCapDict = merge(ExistingCapDict,ScenarioCapacityDict)

# # Fix the capacity of existing generation and storage technologies
for (n,t) in keys(FixCapDict)
      if (n,t) in dtr.sets[:Nodes_Techs]
            JuMP.fix(N_TECH[(n,t)], FixCapDict[(n,t)]; force=true)
      end
      if (n,t) in dtr.sets[:Nodes_Storages]
            JuMP.fix(N_STO_P[(n,t)], FixCapDict[(n,t)] ; force=true)
            JuMP.fix(N_STO_E[(n,t)], MaxEtoP_ratio[n,t]*FixCapDict[(n,t)]; force=true)
      end
end

# Coal_ExistingCapDict = Dict([x for x in par[:ExistingCapacity]
#       if x.first in keys(dvalmatch(par[:FuelType],r"Coal")) && par[:Status][x.first] == "GenericExisting"])
#
# # # (Re)Fix the capacity of existing coal generation at a certain fraction `coal_adjust`
# for (n,t) in keys(Coal_ExistingCapDict)
#       if (n,t) in dtr.sets[:Nodes_Techs]
#             JuMP.fix(N_TECH[(n,t)], dtr.settings[:coal_adjust]*ExistingCapDict[(n,t)]; force=true)
#       end
# end

# %% Hydro reservoir constraints:

inflows_table = CSV.read(joinpath(trace_read_path,"HydroInflowsGWh.csv"))
inflows = stack(inflows_table, variable_name=:Month)

HydroInflow = Dict{Tuple{String,Int64},Float64}()
for (sch,mth,v) in eachrow(inflows[!,[:Scheme,:Month,:value]])
      HydroInflow[sch,parse(Int,String(mth))] = 1000*v  # NOTE: Conversion from GWh to MWh!
end

# Create a correspondence of Hours to the month of the year.
Hours = dtr.sets[:Hours]
HtoM = Dict(eachrow(CSV.read(joinpath(trace_read_path,"HoursToMonths.csv"))));
dtr.parameters[:HtoM] = HtoM

SchemeToRegion = Dict{String,Vector{Tuple{String,String}}}()

Schemes = ["Snowy","Tasmania","FarNorthQLD","OvensMurray"]

SchemeToRegion["Snowy"] = [("N8", "Hydro_Exi")]
SchemeToRegion["Tasmania"] = [("T2", "Hydro_Exi"),("T3", "Hydro_Exi")]
SchemeToRegion["FarNorthQLD"] = [("Q1", "Hydro_Exi")]
SchemeToRegion["OvensMurray"] = [("V1", "Hydro_Exi")]

@constraint(dtr.model, HydroInflowsLimit[sch=Schemes, month=1:12],
      sum(G[(z,t),h] for (z,t) in SchemeToRegion[sch]
                     for h in Hours if HtoM[h] == month)
            <= HydroInflow[sch,month]
);


# Set an upper bound on each technology as a certain fraction `peak_factor` of the peak demand
# df_tech = dfDict["tech"][!,[:Technologies, :Status]] |> dropmissing
# tech_status = Dict(zip(df_tech[!,:Technologies],df_tech[!,:Status]))
#
# df_nodes = dfDict["nodes"] |> dropmissing
# node_DemReg = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:DemandRegion]))
#
# df_node_tech = @linq dfDict["map_node_tech"] |>
#                   where(:IncludeFlag .== 1) |>
#                   select(:Nodes, :Technologies) |>
#                   dropmissing
#
# df_node_tech[!,:Status] = map(x -> tech_status[x], df_node_tech[!,:Technologies])
# df_node_tech[!,:DemandRegion] = map(x -> node_DemReg[x], df_node_tech[!,:Nodes])
#
# for x in eachrow(df_node_tech)
#       if x.Status == "NewEntrant"
#             # println(copy(x))
#             JuMP.set_upper_bound(N_TECH[(x.Nodes,x.Technologies)],
#                   dtr.settings[:peak_factor]*Peaks[x.DemandRegion])
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
println(dtr.settings)
solve_model!(dtr)
resultsIndex = generate_results!(dtr)

sets = dtr.sets
mod = dtr.model
par = dtr.parameters
res = dtr.results

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

using Serialization

# solved_dtr = copy(dtr)
# solved_dtr.model = []
Serialization.serialize(joinpath(resultsdir,results_filename), res)

include("src/output.jl")

# %% Merge results with other runs
# # include("merge.jl")
# post_process_results(rdir)
#
# %% Filtering results

# Split dataframes with Tuple-type columns using split_df_tuple(), e.g.
# split_df_tuple(res[:FLOW],:Arcs,[:From,:To])
# split_df_tuple(res[:N_TECH],:Nodes_Techs,[:Nodes,:Techs])

# df = res[:G]
df_aug = split_df_tuple(res[:G],:Nodes_Techs,[:Nodes,:Techs])
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

# df_flow = res[:FLOW]
df_flow = split_df_tuple(res[:FLOW],:Arcs,[:From,:To])

# Augment with a column showing the Demand Regions of the flow's nodes:
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
L = 1:168
# L = 1:336
# L = 5000:5336

gr()
# gr(size=(5000,2000))
gr(size=(3000,600))
# plotly()
# plotly(size=(3000,600))

# DemandReg = "TAS1"
# p = plot(layout=grid(5,1, height=4*[0.1,0.1,0.1,0.1,0.1]),margin=5mm);

for (count, DemandReg) in enumerate(["NSW1", "QLD1", "VIC1", "SA1", "TAS1"])
            # Load = dtr.parameters[:Load]
      Demand = @where(dfDict["load"],:DemandRegion .== DemandReg)
      df_plot = dfStates[DemandReg]
      Techs = [Symbol(i) for i in DataFrames.unique(copy(df_plot[!,:Techs]))]
      NumTechs = length(Techs)
      # reTechs = reshape(Techs,NumTechs,1)

      df_unstack = unstack(df_plot,:Techs,:Level)
      p = plot()
      # @df df_all[L,:] groupedbar(HOURS[L], cols(1:NumTech),  #cols(NumTech:-1:1),
      @df df_unstack[L,Techs] groupedbar!(p, Hours[L], cols(Techs),
          # subplot=count,
          margin=10mm,
          title=DemandReg,
          xlabel="Time",
          fillalpha=0.5,linealpha=0.1,
          bar_position=:stack,
          legend=:best,  # `:none`, `:best`, `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
          color_palette=:balance) # phase delta rainbow inferno darkrainbow colorwheel

      plot!(Hours[L], [Demand[L,:Load]],label="Demand",
            # subplot=count,
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
            # subplot=count,
            line=4, linecolour=:red,
            xtickfont = font(10, "Courier"),
            xlabel="Time (hr)",
            ylabel="Generation (MW)",
            margin=5mm
            )
      display(p)
end

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
