## Code to load DieterModel names (sets, parameters) into local scope.

# periods = dtr.settings[:periods]

FuelCost = dtr.parameters[:FuelCost]      # FuelCost in currency/MWh-thermal-input
Efficiency = dtr.parameters[:Efficiency]   # Efficiency is unitless (MWh-output/MWh-thermal-input)
CarbonContent = dtr.parameters[:CarbonContent] # CarbonContent in units of t-CO2/MWh-thermal-input
CarbonPrice = dtr.settings[:co2]            # Price on carbon in currency/t-CO2
VariableCost = dtr.parameters[:VariableCost]  # VariableCost in currency/MWh-output


## The following should (largely) match the declaration code in the file "model.jl"

hoursInYear = 8760

Timestep = dtr.settings[:timestep]

periods = round(Int,hoursInYear*(2/Timestep))
Hours = Base.OneTo(periods)
Hours2 = Hours[2:end]
# time_ratio relates generation levels in one time-step (e.g 1/2-hourly) to energy in MWh on an hourly basis
# Used to compare energy to capacity in MW. A 1/2-hourly resolution means time_ratio = 1//2
time_ratio = hoursInYear//length(Hours)
dtr.settings[:time_ratio] = time_ratio
# corr_factor = 1//time_ratio  # length(Hours)//hoursInYear

# dtr.sets[:Hours] = Hours

cost_scaling = dtr.settings[:cost_scaling]

# Set definitions

# # Available sets, but not used directly in model:
    # Technologies = dtr.sets[:Technologies]  # Generation technologies
    # Storages = dtr.sets[:Storages]          # Storage technologies
    # Renewables = dtr.sets[:Renewables]           # Renewable generation technologies
    # Conventional = dtr.sets[:Conventional]       # Conventional generation technologies
    # Dispatchable = dtr.sets[:Dispatchable]       # Dispatchable generation technologies
    # NonDispatchable = dtr.sets[:NonDispatchable] # Non-dispatchable generation technologies

    # Nodes = dtr.sets[:Nodes]

TxZones = dtr.sets[:TxZones]
REZones = dtr.sets[:REZones]

DemandZones = dtr.sets[:DemandZones]
DemandRegions = dtr.sets[:DemandRegions]
node2DemReg = dtr.parameters[:node_demreg_dict]

# Mapping set definitions
Arcs = dtr.sets[:Arcs]
Arcs_From = dtr.sets[:Arcs_From]
Arcs_REZones = dtr.sets[:Arcs_REZones]  # Connections between REZones and TxZones

ArcREZpairing = dtr.parameters[:ArcREZpairing]
Arcs_REZ_Pairs = intersect(Arcs,keys(ArcREZpairing))  # Filter pairings by allowed arcs

Nodes_Techs = dtr.sets[:Nodes_Techs]
Nodes_Storages = dtr.sets[:Nodes_Storages]

Nodes_Renew = dtr.sets[:Nodes_Renew]
# Nodes_Conven = dtr.sets[:Nodes_Conven]
Nodes_Dispatch = dtr.sets[:Nodes_Dispatch]
Nodes_NonDispatch = dtr.sets[:Nodes_NonDispatch]

# Subset of Nodes_Techs with availability traces
Nodes_Avail_Techs = dtr.sets[:Nodes_Avail_Techs]

# Test (TODO - relocate elsewhere):
if !isempty(intersect(Nodes_Dispatch,Nodes_Avail_Techs))
    @warn "The sets Nodes_Dispatch and Nodes_Avail_Techs are not disjoint, and are required to be for correctness."
end

Nodes_Types = dtr.sets[:Nodes_Types]
Nodes_Promotes = dtr.sets[:Nodes_Promotes]

Nodes_Demand = dtr.sets[:Nodes_Demand]

## Electric Vehicles
EV = dtr.sets[:ElectricVehicles]

## Building heating
BU = dtr.sets[:BuildingType]  # Building archtypes
HP = dtr.sets[:HeatingType]   # Heating combination type

## Power-to-gas with Hydrogen
    # P2G = dtr.sets[:P2G]
    # G2P = dtr.sets[:G2P]
    # GasStorages = dtr.sets[:GasStorages]
Nodes_P2G = dtr.sets[:Nodes_P2G]
Nodes_G2P = dtr.sets[:Nodes_G2P]
Nodes_GasStorages = dtr.sets[:Nodes_GasStorages]

# Parameter definitions

MarginalCost = dtr.parameters[:MarginalCost] # Units: currency/MWh; Marginal (variable) cost per unit of generated energy
InvestmentCost = dtr.parameters[:InvestmentCost] # Units: currency/MW; Investment cost per unit of generation power capacity
InvestmentCostPower = dtr.parameters[:InvestmentCostPower] # Units: currency/MW; Annualised investment cost per unit of storage power capacity
InvestmentCostEnergy = dtr.parameters[:InvestmentCostEnergy] # Units: currency/MWh; Annualised investment cost per unit of storage energy capacity
FixedCost = dtr.parameters[:FixedCost] # Units: currency/MW; Fixed cost per unit of power capacity
Availability = dtr.parameters[:Availability] # Units: none; Available energy from renewables as fraction in [0,1] of installed capacity within a time-interval
CapacityDerating = dtr.parameters[:CapacityDerating] # Units: [0,1]; Reduction in nominal capacity as fraction in [0,1] of installed capacity within a time-interval
PeakContribution = dtr.parameters[:PeakContribution] # Units: [0,1]; Fraction of nominal capacity available from renewable tech. at peak demand.
MaxCapacity = dtr.parameters[:MaxCapacity] # Units: MW; Maximum installable capacity - power
MaxEnergy = dtr.parameters[:MaxEnergy] # Units: MWh; Maximum installable storage capacity energy
MaxEtoP_ratio = dtr.parameters[:MaxEnergyToPowerRatio] # Units: hours; Maximum ratio of stored energy to power delivery
Efficiency = dtr.parameters[:Efficiency] # Units: [0,1]; Combustion/Storage roundtrip efficiency
StartLevel = dtr.parameters[:StartLevel] # Units: [0,1]; Initial storage level as fraction of storage energy installed

CarbonContent = dtr.parameters[:CarbonContent] # Units: t-CO2/MWh-thermal; CO2 equivalent content per unit fuel used by tech.
CarbonBudget = dtr.settings[:carbon_budget] # Units: t-CO2

InertialSecs = dtr.parameters[:InertialCoeff]
InertialSecsSto = dtr.parameters[:InertialCoeffSto]
InertiaMinThreshold = dtr.parameters[:InertiaMinThreshold]
InertiaMinSecure = dtr.parameters[:InertiaMinSecure]
RequireRatio = dtr.parameters[:RequireRatio]
SynConCapCost = dtr.parameters[:SynConCapCost]

OperatingReserve = dtr.parameters[:OperatingReserve] # Units: MW; Regional operating reserve requirement (e.g. within state)
MinStableGen = dtr.parameters[:MinStableGen] # Units: [0,1]; Minimum stable operational level for generation as fraction of Capacity in MW.

# WindLimit = dtr.parameters[:WindLimit]
# SolarLimit = dtr.parameters[:SolarLimit]
TotalBuildCap = dtr.parameters[:TotalBuildCap]

ExpansionLimit_REZ = dtr.parameters[:ExpansionLimit]
ExpansionLimit_Tx = dtr.parameters[:ExpansionLimit_Tx]
TransferCapacity = dtr.parameters[:TransferCapacity] # Units: MW; Interconnector power transfer capability

# REZoneExpCost = dtr.parameters[:REZoneExpansionCost] # Units: currency/MW; REZ connection power transfer expansion cost
# TxConnExpCost = dtr.parameters[:TxZoneExpansionCost] # Units: currency/MW; Interconnector power transfer expansion cost

InvestmentCostTransExp = filter(x -> !ismissing(x[2]), dtr.parameters[:InvestmentCostTransExp])
InvestmentCostREZ_Exp = filter(x -> !ismissing(x[2]) && x[1] in Arcs_REZones, dtr.parameters[:InvestmentCostTransExp])

Load = dtr.parameters[:Load] # Units: MWh per time-interval; wholesale energy demand within a time-interval (e.g. hourly or 1/2-hourly)
NegOpDemand = dtr.parameters[:NegOpDemand] # Units: MWh per time-interval; wholesale energy produced from behind-the-meter generation (e.g. PV or EV)

CurtailmentCost = dtr.settings[:cu_cost] # Units: currency/MWh; Cost per unit of generated energy that is curtailed

LoadIncreaseCost = dtr.parameters[:LoadIncreaseCost] # Units: $/MW; Load change costs for changing generation upward
LoadDecreaseCost = dtr.parameters[:LoadDecreaseCost] # Units: $/MW; Load change costs for changing generation downward

MinimumRenewShare = dtr.settings[:min_res_system]

SolarLimit = dtr.parameters[:SolarLimit]
WindLimitOnshore = dtr.parameters[:WindLimitOnshore]