## Running code and tests for Dieter.jl

# %% Load packages
using Dieter
import Dieter: parse_file, parse_nodes!, parse_base_technologies!, parse_storages!, parse_load!, parse_availibility!
import Dieter: initialise_set_relation_data!, parse_set_relations!,parse_arcs!, calc_base_parameters!, parse_extensions!
import Dieter: dvalmatch, dkeymatch, split_df_tuple

using JuMP
# import MathOptInterface
import CPLEX

using DataFrames
using DataFramesMeta
import DBInterface
import SQLite
import CSV
import TimeSeries
import Dates
import Serialization
# using Tables
# import XLSX

# %% Mappings from Scenario designations to parameters
Scenario_Number_Dict = Dict(
      "Scen1_BAU" => 1,  # Current Trends
      "Scen2_DDC" => 2,  # Deep Decarbonisation
      "Scen3_PRP" => 3,  # Prosumer Power
      "Scen4_DID" => 4,  # De-Industrialisation Death Spiral
      "Scen5_MAG" => 5,  # Make Australia Great Again
      "Scen6_NGA" => 6,  # NSW Goes It Alone
      "Scen7_DCU" => 7,  # The DC Universe
      "Scen8_LCE" => 8   # Australia's Low Cost Energy Advantage
)

Scen_ISP_Map = Dict("Scen1_BAU" => "Central", "Scen2_DDC" => "HighVRE")

Scen_FuelCost_Map = Dict("Scen1_BAU" => "Neutral", "Scen2_DDC" => "Fast")

Scenario_BattVPP_Dict = Dict(
      "Scen1_BAU" => 2,  # Current Trends
      "Scen2_DDC" => 3,  # Deep Decarbonisation
      "Scen3_PRP" => 4,  # Prosumer Power
      "Scen4_DID" => 1,  # De-Industrialisation Death Spiral
      "Scen5_MAG" => 1,  # Make Australia Great Again
      "Scen6_NGA" => 3,  # NSW Goes It Alone
      "Scen7_DCU" => 2,  # The DC Universe
      "Scen8_LCE" => 3   # Australia's Low Cost Energy Advantage
)

Scen_h2_setting = Dict("Scen1_BAU" => missing, "Scen2_DDC" => 0)
# %% Scenario Settings (to customise by modeller)

run_timestamp = "$(Dates.Date(Dates.now()))-H$(Dates.hour(Dates.now()))"

ScenarioName = "Scen1_BAU"
# ScenarioName = "Scen2_DDC"

# Specfied Year for the scenario setting:
# ScenarioYear = 2030
ScenarioYear = 2050

ScYr_Sym = Symbol("FYE$ScenarioYear")  # Scenario year symbol
ScenarioNumber = Scenario_Number_Dict[ScenarioName]

# # If FixExistingCapFlag is `true`, then fix existing capacity to given values,
# otherwise if `false` just use the existing capacity as an upper bound.
FixExistingCapFlag = true

NoNewGas = false
NoNewDistillate = false

Note = "Testing"

scen_settings = Dict{Symbol,Any}()

scen_types = Dict{Symbol, String}()
BattEnergyType = "N_BattEnergy"
HydPumpEnergyType = "N_HydPumpEnergy"

H2ElectrolyserType = "N_Electrolyser"
H2RecipEngType = "N_RecipH2"

scen_types[:BattEnergyType] = BattEnergyType
scen_types[:HydPumpEnergyType] = HydPumpEnergyType
scen_types[:H2ElectrolyserType] = H2ElectrolyserType
scen_types[:H2RecipEngType] = H2RecipEngType

scen_settings[:scen_types] = scen_types

# Year parameters
WeatherYear = 2019 # e.g if 2019, this is Financial year 2018-2019
ReferenceYear = WeatherYear # Year of data set to use for renewable traces
TraceYear = 2030 # Which year to use from the ReferenceYear trace dataset

scen_settings[:weather_year] = WeatherYear
scen_settings[:trace_year] = TraceYear
# Technology

scen_settings[:scen_name] = ScenarioName
scen_settings[:scen_year] = ScenarioYear
scen_settings[:scen] = run_timestamp*"-$(ScenarioName)-ScYr$(ScenarioYear)-$(Note)"

scenario_timestamp = scen_settings[:scen]

scen_settings[:interest] = 0.059
scen_settings[:cost_scaling] = 1 # 1.0e-6
# Modify the :min_res setting over [0,100] and rerun to see comparison.
scen_settings[:min_res] = 10
scen_settings[:ev] = missing
scen_settings[:heat] = missing
scen_settings[:h2] = Scen_h2_setting[ScenarioName]  # missing -> H2 not included, any number -> H2 included


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

TxZones = dtr.sets[:TxZones]
REZones = dtr.sets[:REZones]

DemandZones = dtr.sets[:DemandZones]
DemandRegions = dtr.sets[:DemandRegions]

node2DemReg = dtr.parameters[:node_demreg_dict]

# e.g. fileDict["tech"] = joinpath(datapath,"base","technologies.csv")
dfDict["tech"] = parse_file(fileDict["tech"]; dataname=dataname)

# e.g. fileDict["storage"] = joinpath(datapath,"base","storages.csv")
dfDict["storage"] = parse_file(fileDict["storage"]; dataname=dataname)

# %% Relations (i.e set-to-set correspondences)
dfDict["map_node_demand"] = parse_file(fileDict["map_node_demand"]; dataname=dataname)
dfDict["map_node_tech"] = parse_file(fileDict["map_node_tech"]; dataname=dataname)
dfDict["map_node_storages"] = parse_file(fileDict["map_node_storages"]; dataname=dataname)
dfDict["arcs"] = parse_file(fileDict["arcs"]; dataname=dataname)

# %% Add VPP battery type (included as existing capacity)

# Split VPP state capacities evenly across TxZones in that state:
RegionTxSplit = Dict("NSW1" => 4, "VIC1" => 4, "QLD1" => 4, "SA1" => 3, "TAS1" => 1)

#  Virtual power plant (VPP) battery uptake data:
fileDict["vpp_uptake"] = joinpath(datapath,"base","vpp_uptake.sql")
dfDict["vpp_uptake"] = parse_file(fileDict["vpp_uptake"]; dataname=dataname)

df_vpp_data = @linq dfDict["vpp_uptake"] |>
            where(
                  :Year .== ScenarioYear,
                  :VPP_Scenario .== Scenario_BattVPP_Dict[ScenarioName]
            ) |>
            select(:Region, :ValueType, :Data)

df_vpp_data = unstack(df_vpp_data,:ValueType,:Data)

df_vpp_data = @byrow! df_vpp_data begin
           @newcol Ratio::Array{Float64}
           :Ratio = :Energy/:Power
       end

VPP_Power_Dict = Dict(zip(df_vpp_data[!,:Region],df_vpp_data[!,:Power]))
VPP_Ratio_Dict = Dict(zip(df_vpp_data[!,:Region],df_vpp_data[!,:Ratio]))

# Copy off on the New Entrant Battery Storage
df_new_vpp = @where(dfDict["storage"], :Storages .== "BattStor_New")

df_add_vpp = @byrow! df_new_vpp begin
                  # Annotate with the Demand Region:
                  @newcol DemandRegion::Array{String}
                  :DemandRegion = node2DemReg[:Region]
                  # Change the names and capacties
                  :Storages = "BattVPP_Exi"
                  :Status = "GenericExisting"
                  :ExistingCapacity = VPP_Power_Dict[:DemandRegion]/RegionTxSplit[:DemandRegion]
                  :MaxEnergyToPowerRatio = VPP_Ratio_Dict[:DemandRegion]
            end

df_add_vpp_map = @byrow! df_add_vpp begin
                  # Include all new storage tech.:
                  @newcol IncludeFlag::Array{Int64}
                  :IncludeFlag = 1
            end

rename!(df_add_vpp_map, Dict(:Region => :Nodes, :Storages => :Technologies))
select(df_add_vpp_map, :Nodes,:Technologies, :IncludeFlag)

# Add into data:
append!(dfDict["storage"], select(df_add_vpp, Not(:DemandRegion)))
append!(dfDict["map_node_storages"], select(df_add_vpp_map, :Nodes,:Technologies, :IncludeFlag))

# %% Additional modifications to included technologies:

# # Remove technologies if directed:
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

# Remove new distillate if specified:
if NoNewDistillate
      DistillateTech_New = @linq dfDict["tech"] |>
                  where(:Status .== "NewEntrant", occursin.(r"Distillate",:FuelType)) |>
                  select(:Technologies)
      DistillateTech_New = (unique(DistillateTech_New,:Technologies))[!,:Technologies]

      dfDict["map_node_tech"] = @byrow! dfDict["map_node_tech"] begin
            if in(:Technologies,DistillateTech_New)
                  :IncludeFlag = 0
            end
      end
end

# %%  Relations to determine which (Nodes,Technologies) pairs are included in the model:
rel_node_tech = create_relation(dfDict["map_node_tech"],:Nodes,:Technologies,:IncludeFlag)
rel_node_storages = create_relation(dfDict["map_node_storages"],:Nodes,:Technologies,:IncludeFlag)

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

if dtr.settings[:min_res]["TAS1"] > 100
      dtr.settings[:min_res]["TAS1"] = 100
end
# %% New capital / overnight costs
fileDict["capital_costs"] = joinpath(datapath,"base","capital_costs.sql")
dfDict["capital_costs"] = parse_file(fileDict["capital_costs"]; dataname=dataname)

# Rename the ScenarioYear column as `:ScenCostPower`:
df_cap_costs = DataFrames.rename!(
      dfDict["capital_costs"][!,[:Technologies, :Scenario, ScYr_Sym]],
      Dict(ScYr_Sym => :ScenCostPower)
      )
# Filter by the current Scenario:
df_cap_costs = @linq df_cap_costs |> where(:Scenario .== Scen_ISP_Map[ScenarioName])
# Read into model parameters:
params_cap_costs = Dieter.map_idcol(df_cap_costs, [:Technologies], skip_cols=Symbol[:Scenario])
for (k,v) in params_cap_costs Dieter.update_dict!(dtr.parameters, k, v) end

# # Energy storage technologies:
fileDict["cost_energy_storage"] = joinpath(datapath,"base","cost_energy_storage.sql")
dfDict["cost_energy_storage"] = parse_file(fileDict["cost_energy_storage"]; dataname=dataname)

# Rename the ScenarioYear column with ScenCostEnergy:
df_costES = DataFrames.rename!(
      dfDict["cost_energy_storage"][!,[:Technologies, :Scenario, ScYr_Sym]],
      Dict(ScYr_Sym => :ScenCostEnergy)
      )
# Filter by the current Scenario:
df_costES = @linq df_costES |> where(:Scenario .== Scen_ISP_Map[ScenarioName])
# Read into model parameters:
params_costES = Dieter.map_idcol(df_costES, [:Technologies], skip_cols=Symbol[:Scenario])
for (k,v) in params_costES Dieter.update_dict!(dtr.parameters, k, v) end

# # Hydrogen technologies:
fileDict["cost_h2"] = joinpath(datapath,"h2","cost_h2.sql")
dfDict["cost_h2"] = parse_file(fileDict["cost_h2"]; dataname=dataname)

df_cost_h2 = DataFrames.rename!(
      dfDict["cost_h2"][!,[:Technologies, :Scenario, ScYr_Sym]],
      Dict(ScYr_Sym => :ScenCostPower)
      )
# Filter by the current Scenario:
df_cost_h2 = @linq df_cost_h2 |> where(:Scenario .== Scen_ISP_Map[ScenarioName])
# Read into model parameters:
params_cost_h2 = Dieter.map_idcol(df_cost_h2, [:Technologies], skip_cols=Symbol[:Scenario])
for (k,v) in params_cost_h2 Dieter.update_dict!(dtr.parameters, k, v) end

# The above should create parameters in dtr.parameters[:ScenCostPower] for H2 tech.
# %% Synchronous condensers
# Synchronous condenser with flywheel - cost in $ per inertial units of `MWs`

SynConOvernightCost =  round(1.0E6*185.2/4400, digits=2) # Units:  $/MWs
# ElectraNet project of $185.2 million, for 4400 MWs of inertia via 4 sync. cons.
# ElectraNet, "Main grid system strength project: Contingent project application", p. 21, 28 June 2019, p. 19.
# All dollar amounts in this document are in real, $2017–18 in line with the ElectraNet's revenue determination unless otherwise stated.

LifetimeSynCon = 30 # InvestRecovery # InvestLifetime 40 # Source: GHD; Economic life for ElectraNet synchronous condensers - ElectraNet 28 June 2019
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

df_techscen = @linq dfDict["tech_scenario"] |> where(:ScenarioName .== ScenarioName)

df_techscen = DataFrames.rename!(
      df_techscen[!,[:Region, :TechID, ScYr_Sym]],
      Dict(:TechID => :Technologies, ScYr_Sym => :ScenarioCapacity)
      )


# Filter the technology capacities if they are included in our region/technology list:
df_techscen_inc = @linq df_techscen |>
                    where( rel_node_tech.(:Region, :Technologies) .== true )
# Filter the technology capacities if they are not in our region/technology list:
df_techscen_exc = @linq df_techscen |>
                    where( rel_node_tech.(:Region, :Technologies) .== false )

# TODO: check/review the excluded list of technology capacities to make sure none are excluded unnecessary
# This includes expanding this to also include storage technology overwrites if needed.

# The following adds the parameter `ScenarioCapacity` to `dtr.parameters`:
for (k,v) in Dieter.map_idcol(
            df_techscen_inc,
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
dtr.settings[:co2] = Scen_co2[!,:CarbonPrice][1]  # Units in $/t-CO2
Mt_To_t = 1000 #  convert Mega-tonnes to tonnes
dtr.settings[:carbon_budget] = Mt_To_t*Scen_co2[!,:CarbonBudget][1]  # Data units in Mt-CO2, convert to t-CO2

# Add a Carbon Content column with zero values
if !in(:CarbonContent, propertynames(dfDict["tech"]))
      insertcols!(dfDict["tech"], size(dfDict["tech"])[2], :CarbonContent => zeros(size(dfDict["tech"])[1]))
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

BattEnergyType = dtr.settings[:scen_types][:BattEnergyType]
HydPumpEnergyType =  dtr.settings[:scen_types][:HydPumpEnergyType]

BattEnergyCost = dtr.parameters[:ScenCostEnergy][BattEnergyType]
HydPumpEnergyCost = dtr.parameters[:ScenCostEnergy][HydPumpEnergyType]
dfDict["storage"] = @byrow! dfDict["storage"] begin
       if :Status == "NewEntrant" && rel_node_storages(:Region,:Storages) == true
             :OvernightCostPower = ScenCostPower[:Storages]
                  if occursin(r"Battery",:TechType)
                       :OvernightCostEnergy = BattEnergyCost
                       :MaxEnergyToPowerRatio = 12.0
                  elseif occursin(r"Pump Storage",:TechType)
                       :OvernightCostEnergy = HydPumpEnergyCost
                 end
       end
 end

 dfDict["storage"] = @byrow! dfDict["storage"] begin
       if :Storages == "HydPump_New"
              :MaxEnergyToPowerRatio = 200.0
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
@linq dfDict["storage"] |> where(:Storages .== "HydPump_New") |> select(:Region, :ExistingCapacity, :MaxEnergyToPowerRatio)
@linq dfDict["storage"] |> where(occursin.(r"Battery",:TechType)) |> select(:Region, :ExistingCapacity, :MaxEnergyToPowerRatio)

# # Battery of the Nation...

ConnectCost = dtr.parameters[:ConnectCost]

# Add connection costs specific to regions to the capital cost of plants:
dfDict["tech"] = @byrow! dfDict["tech"] begin
      # @newcol IncCost::Array{Float64}
      if :Status == "NewEntrant" && rel_node_tech(:Region,:Technologies) == true
            :OvernightCostPower = :OvernightCostPower + ConnectCost[(:Region,:Technologies)]
      end
end

dfDict["storage"] = @byrow! dfDict["storage"]  begin
      # @newcol IncCost::Array{Float64}
      if :Status == "NewEntrant" && rel_node_storages(:Region,:Storages) == true
            :OvernightCostPower = :OvernightCostPower + ConnectCost[(:Region,:Storages)]
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

# load_trace_datapath = joinpath(trace_read_path,"Demand_Traces","Load_FYE$(Demand_Year).csv")
load_trace_datapath = joinpath(trace_read_path,"Demand_Traces")
load_trace_filename = "wy_$(WeatherYear)_scen_$(ScenarioNumber)_AllRegions.csv"
load_trace_filepath = joinpath(load_trace_datapath,load_trace_filename)

@assert isfile(load_trace_filepath) "File $(load_trace_filename) not found in $(load_trace_datapath)"

# The operational demand table has time-stamps for keys and one column per region.
# This can be read using TimeSeries.jl as TimeSeries.readtimearray()
ta_OpDem = TimeSeries.readtimearray(load_trace_filepath)
# This could also be read using CSV.jl as CSV.File(...)
# df_OpDem = DataFrame(CSV.File(load_trace_filepath))

# Filter by Scenario Year
Financial_Year_Array = [Dates.DateTime(ScenarioYear-1,07,01,00,00,00):Dates.Minute(30):Dates.DateTime(ScenarioYear,06,30,23,30,00)]

ta_OpDem = ta_OpDem[Financial_Year_Array...]

df_OpDem = DataFrame(ta_OpDem)

# Transformations of the raw data:

# If timestep = 2, we obtain an hourly approximation of half-hour data by sampling on every second data point:
df_OpDem = df_OpDem[1:timestep:end,propertynames(df_OpDem)]

# ta_OpDem = TimeArray(df_OpDem, timestamp = :datetime)
dropmissing!(df_OpDem)
insertcols!(df_OpDem, 1, :TimeIndex => collect(1:nrow(df_OpDem)))
select!(df_OpDem, Not(:timestamp))

# reg_rename_dict = Dict(:nsw1 => :NSW1,
#                        :qld1 => :QLD1,
#                        :sa1 => :SA1,
#                        :tas1 => :TAS1,
#                        :vic1 => :VIC1)
# region_names = sort(collect(values(reg_rename_dict)))
# DataFrames.rename!(df_OpDem, reg_rename_dict...)

df_OpDem_stack = stack(df_OpDem, Symbol.(DemandRegions), :TimeIndex)
                        # variable_name=:DemandRegion, value_name=:Load)
dfDict["load"] = @linq df_OpDem_stack |>
                  transform(DemandRegion = String.(:variable), Load = :value) |>
                  select(:TimeIndex,:DemandRegion,:Load)

df_load_neg = @where(dfDict["load"], :Load .< 0)

if size(df_load_neg)[1] !== 0
      @info "Negative demand values found: $(size(df_load_neg)[1]) values."
      @info "Setting negative values to 0."
      df_threshold = @byrow! dfDict["load"] begin
                    @newcol NegLoad::Array{Float64}
                    if :Load .< 0
                       :NegLoad = -(:Load)
                       :Load = 0
                    else
                        :NegLoad = 0
                    end
       end
       dfDict["load"] = select(df_threshold, Not(:NegLoad))
end


# %% Calculate total and peak demand
# DemandRegions = dtr.sets[:DemandRegions]
YearlyEnergy = Dict{String,Float64}()
Peaks = Dict{String,Float64}()
for dr in DemandRegions
      YearlyEnergy[dr] = 1e-3*sum(@where(dfDict["load"], :DemandRegion .== dr)[!,:Load])
      # 1e-3*sum(values(dkeymatch(dtr.parameters[:Load],Regex(dr),1)))
      Peaks[dr] = maximum(@where(dfDict["load"], :DemandRegion .== dr)[!,:Load])
      # maximum(values(dkeymatch(dtr.parameters[:Load],Regex(dr),1)))
end

dtr.parameters[:Peaks] = Peaks

# %% Scenario modification of Demand :

# Transform from Demand Region to a sub-regional Demand Zone with a certain split:
fileDict["demand_split"] = joinpath(datapath,"base","demand_split.sql")
dfDict["demand_split"] = parse_file(fileDict["demand_split"]; dataname=dataname)
# Expected columns: | DemandZone | DemandRegion | DemandShare |

# Join the Demand Zones along the superset of Demand Regions:
df_load_join = innerjoin(dfDict["load"],dfDict["demand_split"], on = :DemandRegion)
# Combine the demand with the shared split:
dfDict["load_share"] = @linq df_load_join |>
                  transform(LoadShare = :Load .* :DemandShare) |>
                  select(:TimeIndex, :DemandZone, :LoadShare)
rename!(dfDict["load_share"], Dict(:DemandZone => :Nodes, :LoadShare => :Load))


#  Note : not needed if input data is correct without modification
# fileDict["demand_scenario"] = joinpath(datapath,"base","demand_scenario.sql")
# dfDict["demand_scenario"] = parse_file(fileDict["demand_scenario"]; dataname=dataname)

# If dfDict["demand_scenario"] corresponds to a scaling of demand,
# │ Row │ ScenarioName │ Region │ FYE2020 │ FYE2021  │ FYE2022  │ FYE2023  │ FYE2024  │
# │     │ String       │ String │ Float64 │ Float64  │ Float64  │ Float64  │ Float64  │
# ├─────┼──────────────┼────────┼─────────┼──────────┼──────────┼──────────┼──────────┤
# │ 1   │ Scen1_BAU    │ NSW1   │ 1.0     │ 0.99591  │ 0.991982 │ 0.985298 │ 0.979469 │
# │ 2   │ Scen1_BAU    │ QLD1   │ 1.0     │ 1.01197  │ 1.02199  │ 1.03051  │ 1.03456  │
# ...
# │ 6   │ Scen2_DDC    │ NSW1   │ 1.0     │ 0.99911  │ 0.996888 │ 0.992015 │ 0.988109 │
# │ 7   │ Scen2_DDC    │ QLD1   │ 1.0     │ 1.01334  │ 1.02574  │ 1.03758  │ 1.04457  │
# ...
# we would use code like this:
#=
df_ds = @where(dfDict["demand_scenario"], :ScenarioName .== ScenarioName)
ds_Dict = Dict(eachrow(select(df_ds,:Region,ScYr_Sym)))

dfDict["load"] = @byrow! dfDict["load"] begin
                        :Load = ds_Dict[:DemandRegion]*:Load
                  end
=#
# %%

parse_load!(dtr, dfDict["load_share"])
# parse_load!(dtr, dfDict["load"])


# %% Inertia -  data and constraints

fileDict["inertia_tech"] = joinpath(datapath,"base","inertia_tech.sql")
dfDict["inertia_tech"] = parse_file(fileDict["inertia_tech"]; dataname=dataname)

for (k,v) in Dieter.map_idcol(
            dfDict["inertia_tech"], [:Technologies], skip_cols=Symbol[])
                  Dieter.update_dict!(dtr.parameters, k, v)
end

fileDict["inertia_storage"] = joinpath(datapath,"base","inertia_storage.sql")
dfDict["inertia_storage"] = parse_file(fileDict["inertia_storage"]; dataname=dataname)

for (k,v) in Dieter.map_idcol(
            dfDict["inertia_storage"], [:Region,:Technologies], skip_cols=Symbol[])
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

# RequireRatio is used in the inertia constraint to adapt variables demand profiles to the reference year for inertia levels.
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

wind_traces_path = joinpath(trace_read_path,"REZ_Traces","REZ_Wind_Traces_RefYear$(ReferenceYear)_FYE$(TraceYear).csv")
solar_traces_path = joinpath(trace_read_path,"REZ_Traces","REZ_Solar_Traces_RefYear$(ReferenceYear)_FYE$(TraceYear).csv")

dtr.data["files"]["Wind"] = wind_traces_path
dtr.data["files"]["Solar"] = solar_traces_path

tech_subset = Dict()
tech_subset["Wind"]  = ["WindOn_Exi", "WindOn_New", "WindOff_New"]
tech_subset["Solar"] = ["SolarPV_Exi", "SolarPV_New", "SolThermal_New"]
# alt.: ["SolFixedPV_Exi", "SolLargePV_Exi", "SolLargePV_New", "SolThermal_New"]

for tech_class in ["Wind", "Solar"]
      traces_path = dtr.data["files"][tech_class]
      # Read trace data:
      # df_traces = DataFrame(CSV.File(traces_path))
      # This could also be read using TimeSeries.jl as TimeSeries.readtimearray()
      df_traces = DataFrame(TimeSeries.readtimearray(traces_path))


      # If timestep = 2, we obtain an hourly approximation of half-hour data by sampling on every second data point:
      df_traces = df_traces[1:timestep:end,propertynames(df_traces)]

      # Construct an integer OneTo(n)-like time indexing instead of timestamps:
      len_traces = nrow(df_traces)
      insertcols!(df_traces, 1, :TimeIndex => collect(1:len_traces))
      select!(df_traces, Not(:timestamp))

      for tech_name in tech_subset[tech_class]
            df_trace_mod = copy(df_traces)

            nodes_for_tech = [x[1] for x in dtr.sets[:Nodes_Techs] if x[2] == tech_name] |> sort #  println(nodes_for_tech)

            # Require all necessary traces to exist:
            for rz in dtr.sets[:REZones]
                  if rz in nodes_for_tech && !(Symbol(rz) in propertynames(df_trace_mod))
                        error("Require trace in $(rz) for technology $(tech_name)")
                  end
            # Remove unnecessary traces from data before compiling into model:
                  if Symbol(rz) in propertynames(df_trace_mod) && !(rz in nodes_for_tech)
                        @debug "Tech. not in zone: Removing $(tech_name) trace in $(rz)"
                        select!(df_trace_mod, Not(Symbol(rz)))
                  end
            end

            insertcols!(df_trace_mod, 1, :TechTypeID => repeat([tech_name],len_traces))
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
insertcols!(df_REZ_to_Techs,3, :IncludeFlag => ones(Int,nrow(df_REZ_to_Techs)))
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

# %% Initialise model

# Construct an optimizer factory
# solver = JuMP.with_optimizer(Clp.Optimizer)
# solver = JuMP.with_optimizer(Gurobi.Optimizer)
# solver = JuMP.with_optimizer(CPLEX.Optimizer)
solver = CPLEX.Optimizer
# build_model!(dtr,solver; timestep=-1)
build_model!(dtr,solver,timestep=timestep)

# %% Access model objects for further development:

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
if isempty(NoExpansionREZones)
      @info "No REZ expansions fixed to 0."
else
      for rez in NoExpansionREZones
            @info "Fixing REZ expansion to 0 for $(rez)."
            JuMP.fix(N_RES_EXP[rez], 0; force=true)
      end
end

MaxEtoP_ratio = dtr.parameters[:MaxEnergyToPowerRatio]

ExistingCapDict = filter(x -> x.second > 0, dtr.parameters[:ExistingCapacity])
# ExistingCapDict = filter(x -> !ismissing(x.second), dtr.parameters[:ExistingCapacity])

# Overwrite with any imposed scenario capacities read via "tech_scenario" files:
FixCapDict = merge(ExistingCapDict,ScenarioCapacityDict)

# if FixExistingCapFlag == true
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
# elseif FixExistingCapFlag == false
#       # # Fix the capacity of existing generation and storage technologies
#       for (n,t) in keys(FixCapDict)
#             if (n,t) in dtr.sets[:Nodes_Techs]
#                   JuMP.set_upper_bound(N_TECH[(n,t)], FixCapDict[(n,t)])
#             end
#             if (n,t) in dtr.sets[:Nodes_Storages]
#                   JuMP.set_upper_bound(N_STO_P[(n,t)], FixCapDict[(n,t)])
#                   JuMP.set_upper_bound(N_STO_E[(n,t)], MaxEtoP_ratio[n,t]*FixCapDict[(n,t)])
#             end
#       end
# else
#       @warn "The treatment of existing capacity is undefined."
# end

# Coal_ExistingCapDict = Dict([x for x in par[:ExistingCapacity]
#       if x.first in keys(dvalmatch(par[:FuelType],r"Coal")) && par[:Status][x.first] == "GenericExisting"])
#
# # # (Re)Fix the capacity of existing coal generation at a certain fraction `coal_adjust`
# for (n,t) in keys(Coal_ExistingCapDict)
#       if (n,t) in dtr.sets[:Nodes_Techs]
#             JuMP.fix(N_TECH[(n,t)], dtr.settings[:coal_adjust]*ExistingCapDict[(n,t)]; force=true)
#       end
# end

# %% Compatibility check between Minimum Stable Generation constraint `MinStableGeneration`,
# and the maximum carbon allowance / CarbonBudget constraint `CarbonBudgetLimit`:

# periods = round(Int,Dieter.hoursInYear*(2/timestep))
# time_ratio = Dieter.hoursInYear//periods
MinStableGen = dtr.parameters[:MinStableGen]
FixedStableCap = filter(x -> x[1] in keys(MinStableGen), FixCapDict)
MinStableGenTotal = sum(MinStableGen[k]*FixCapDict[k] for k in keys(MinStableGen))

CarbonBudget = dtr.settings[:carbon_budget]
CarbonContent = dtr.parameters[:CarbonContent]
Efficiency = dtr.parameters[:Efficiency]

MinStableCarbon = sum((CarbonContent[k]/Efficiency[k])*MinStableGen[k]*FixCapDict[k]
                        for k in keys(MinStableGen) )

if MinStableCarbon >= CarbonBudget
      error("Dieter:STABLE-data-check: The minimum level of carbon from min. stable generation exceeds the carbon budget.")
end

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
# df_node_tech = @linq dfDict["map_node_tech"] |>
#                   where(:IncludeFlag .== 1) |>
#                   select(:Nodes, :Technologies) |>
#                   dropmissing
#
# df_node_tech[!,:Status] = map(x -> tech_status[x], df_node_tech[!,:Technologies])
# df_node_tech[!,:DemandRegion] = map(x -> node2DemReg[x], df_node_tech[!,:Nodes])
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
display(dtr.settings)
solve_model!(dtr)
resultsIndex = generate_results!(dtr)

# %% Abbreviations:
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

# solved_dtr = copy(dtr)
# solved_dtr.model = []
Serialization.serialize(joinpath(resultsdir,results_filename), dtr.results)
Serialization.serialize(joinpath(resultsdir,scenario_timestamp*".settings"),dtr.settings)

# res = Serialization.deserialize(joinpath(resultsdir,results_filename))

# timestep = 2
# timestep = dtr.settings[:timestep]
# periods = round(Int,Dieter.hoursInYear*(2/timestep))
# time_ratio = Dieter.hoursInYear//periods

# %% Get the results:

include("src/results_transform.jl")

include("src/results_summary.jl")
# include("src/results_plot.jl")

# resultsdir = outputdir
include("src/write.jl")

# %% Merge results with other runs
# # include("merge.jl")
# post_process_results(rdir)
#

# %% Plot timeseries

# include("src/plot_timeseries.jl")
