## This file is part of Dieter.jl : Model definition

# const hoursInYear = 8760     # Units: hours/yr
# const H2_energy_density = 39405.6 # Units : MWh / kg
# const infeas_cost = 1000     # Currency (cost)

"""
Build the JuMP model describing the optimization problem, specifying the `solver` to use.
The `Timestep` parameter should be 1 for half-hourly steps, 2 for hourly steps,
and 4 for 2-hourly steps. The data must match the time-steps of the `Timestep` parameter.
"""
function build_model!(dtr::DieterModel)

    dtr.model = Model()

    m = dtr.model

    ## Sync. the following with the file "load_names.jl"

    Timestep = dtr.settings[:timestep]

    periods = round(Int,hoursInYear*(2/Timestep))
    Hours = Base.OneTo(periods)
    Hours2 = Hours[2:end]
    # time_ratio relates generation levels in one time-step (e.g 1/2-hourly) to energy in MWh on an hourly basis
    # Used to compare energy to capacity in MW. A 1/2-hourly resolution means time_ratio = 1//2
    time_ratio = hoursInYear//length(Hours)
    dtr.settings[:time_ratio] = time_ratio
    # corr_factor = 1//time_ratio  # length(Hours)//hoursInYear

    dtr.sets[:Hours] = Hours

    cost_scaling = dtr.settings[:cost_scaling]
    loss_factor_tx = dtr.settings[:loss_factor_tx] # e.g 1.06 means a 6% increase in given Demand to account for transmission losses.

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

    # Nodes_Techs_Exi = dtr.sets[:Nodes_Techs_Exi]
    Nodes_Techs_New = dtr.sets[:Nodes_Techs_New]
    Nodes_Storages_New = dtr.sets[:Nodes_Storages_New]

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

    Nodes_RampingTechs = dtr.sets[:Nodes_RampingTechs]

    # Synchronous condensers
    Nodes_SynCons = dtr.sets[:Nodes_SynCons]

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

    InvestmentCostSynCon = dtr.parameters[:InvestmentCostSynCon]

    OperatingReserve = dtr.parameters[:OperatingReserve] # Units: MW; Regional operating reserve requirement (e.g. within state)
    MinStableGen = dtr.parameters[:MinStableGen] # Units: [0,1]; Minimum stable operational level for generation as fraction of Capacity in MW.

    # WindLimit = dtr.parameters[:WindLimit]
    # SolarLimit = dtr.parameters[:SolarLimit]
    TotalBuildCap = dtr.parameters[:TotalBuildCap]

    ExpansionLimit_REZ = dtr.parameters[:ExpansionLimit]
    ExpansionLimit_Tx = dtr.parameters[:ExpansionLimit_Tx]
    TransferCapacity = dtr.parameters[:TransferCapacity] # Units: MW; Interconnector power transfer capability

    InvestmentCostTransExp = filter(x -> !ismissing(x[2]) && x[1] in Arcs, dtr.parameters[:InvestmentCostTransExp]) # Units:  currency/MW; ; REZ capacity expansion costs, amortised to year
    InvestmentCostREZ_Exp  = filter(x -> !ismissing(x[2]) && x[1] in Arcs_REZones, dtr.parameters[:InvestmentCostTransExp]) # Units:  currency/MW; ; Transmission interconnector capacity expansion costs, amortised to year

    ArcWeight = dtr.parameters[:ArcWeight]

    Load = dtr.parameters[:Load] # Units: MWh per time-interval; wholesale energy demand within a time-interval (e.g. hourly or 1/2-hourly)
    NegOpDemand = dtr.parameters[:NegOpDemand] # Units: MWh per time-interval; wholesale energy produced from behind-the-meter generation (e.g. PV or EV)

    CurtailmentCost = dtr.settings[:cu_cost] # Units: currency/MWh; Cost per unit of generated energy that is curtailed

    LoadIncreaseCost = dtr.parameters[:LoadIncreaseCost] # Units: $/MW; Load change costs for changing generation upward
    LoadDecreaseCost = dtr.parameters[:LoadDecreaseCost] # Units: $/MW; Load change costs for changing generation downward

    MaxRampUpPerHour = dtr.parameters[:MaxRampUpPerHour]      # Units MW/h; Maximum ramping up of technology
    MaxRampDownPerHour = dtr.parameters[:MaxRampDownPerHour]  # Units MW/h; Maximum ramping down of technology

    MinimumRenewShare = dtr.settings[:min_res_system]

    CapAdd = dtr.parameters[:CapAdd]

    ## End of sync. with "load_names.jl"

    # Transformations on the minimum renewable share settings:
    min_res_dict = Dict{String,Float64}()
    if isa(dtr.settings[:min_res],Number)
        for n in DemandRegions
            min_res_dict[n] = dtr.settings[:min_res]
        end
    elseif isa(dtr.settings[:min_res],Dict)
        min_res_dict = dtr.settings[:min_res]
        dtr.settings[:min_res] = minimum(collect(values(min_res_dict)))  # overwrite as a number, the smallest share.
        println("The `min_res` setting is now equal to $(dtr.settings[:min_res])")
    else
        error("The parameter `min_res_dict` is not of type `Number` or `Dict`.")
    end
    dtr.settings[:min_res_dict] = min_res_dict
    @info "The `min_res_dict` setting is currently equal to "
    @info min_res_dict
    for (k,v) in min_res_dict
        if !(0 <= v <= 100)
            error("The setting in `min_res_dict` of value $v for Demand Region $k is not in the interval [0,100].")
        elseif (0 < v <= 1)
            @warn "Note: The setting in `min_res_dict` of value $v for Demand Region $k is a percentage, not necessarily a fraction in the interval (0,1]."
        end
    end

    @info "Start of model building:"

    prog = Progress(4, dt=0.01, desc="Building Model...         \n", barlen=30)

# %% * ----------------------------------------------------------------------- *
#    ***** Variable definitions *****
#    * ----------------------------------------------------------------------- *

    shorter =  Dict("Total_cost_objective" => "Z",
                    "Generation_level" => "G",
                    "Generation_upshift" => "G_UP",
                    "Generation_downshift" => "G_DO",
                    "Generation_renewable" => "G_REZ",
                    "Generation_transmission_zones" => "G_TxZ",
                    "Storage_inflow" => "STO_IN",
                    "Storage_outflow" => "STO_OUT",
                    "Storage_level" => "STO_L",
                    "Technology_capacity" => "N_TECH",
                    "Renewable_capacity_expand" => "N_REZ_EXP",
                    "REZ_expand_by_transmission" => "N_REZ_EXP_TX",
                    "Storage_build_energy" => "N_STO_E",
                    "Storage_capacity" => "N_STO_P",
                    "SynCon_capacity" => "N_SYNC",
                    "Inertia_level" => "IN_TxZ",
                    "Internodal_flow" => "FLOW",
                    "Internodal_flow_expansion" => "N_IC_EXP",
                    "H2_power_to_gas" => "H2_P2G",
                    "H2_gas_to_power" => "H2_G2P",
                    "H2_storage_level" => "H2_GS_L",
                    "H2_storage_inflow" => "H2_GS_IN",
                    "H2_storage_outflow" => "H2_GS_OUT",
                    "H2_P2G_capacity" => "N_P2G",
                    "H2_G2P_capacity" => "N_G2P",
                    "H2_storage_capacity" => "N_GS",
                    "EV_charging" => "EV_CHARGE",
                    "EV_discharging" => "EV_DISCHARGE",
                    "EV_charge_level" => "EV_L",
                    "EV_PHEV_fuel_use" => "EV_PHEVFUEL",
                    "EV_infeasible" => "EV_INF",
                    "Heat_storage_level" => "HEAT_STO_L",
                    "Heat_heat_pump_" => "HEAT_HP_IN",
                    "Heat_infeasible" => "HEAT_INF")

    @info "Variable definitions."
    @variables(m, begin
            Z, (base_name=shorter["Total_cost_objective"], lower_bound=0)
            G[Nodes_Techs, Hours], (base_name=shorter["Generation_level"], lower_bound=0) # Units: MWh per time-interval; Generation level - all generation tech.
            G_UP[Nodes_Dispatch, Hours] , (base_name=shorter["Generation_upshift"], lower_bound=0)  # Units: MWh per time-interval; Generation level change up
            G_DO[Nodes_Dispatch, Hours], (base_name=shorter["Generation_downshift"], lower_bound=0) # Units: MWh per time-interval; Generation level change down
            # G_INF[Nodes, Hours], (base_name=shorter["Generation_infeasible"], lower_bound=0) # Units: MWh per time-interval; Infeasibility term for Energy Balance
            # CU[Nodes_NonDispatch, Hours], (base_name=shorter["Curtailment_renewables"], lower_bound=0) # Units: MWh per time-interval; Non-dispatchable curtailment
            G_REZ[REZones,Hours], (base_name=shorter["Generation_renewable"], lower_bound=0) # Units: MWh per time-interval; Generation level - renewable energy zone tech. & storage
            G_TxZ[TxZones,Hours], (base_name=shorter["Generation_transmission_zones"], lower_bound=0) # Units: MWh per time-interval; Generation level - transmission zone tech. & storage
            STO_IN[Nodes_Storages, Hours], (base_name=shorter["Storage_inflow"], lower_bound=0) # Units: MWh per time-interval; Storage energy inflow
            STO_OUT[Nodes_Storages, Hours], (base_name=shorter["Storage_outflow"], lower_bound=0) # Units: MWh per time-interval; Storage energy outflow
            STO_L[Nodes_Storages, Hours], (base_name=shorter["Storage_level"], lower_bound=0) # Units: MWh at a given time-interval; Storage energy level
            N_TECH[Nodes_Techs], (base_name=shorter["Technology_capacity"], lower_bound=0) # Units: MW; Technology capacity built
            N_REZ_EXP[REZones], (base_name=shorter["Renewable_capacity_expand"], lower_bound=0) # Units: MW; Renewable technology transmission capacity built
            N_REZ_EXP_TX[REZones], (base_name=shorter["REZ_expand_by_transmission"], lower_bound=0) # Units: MW; Free REZ transmission capacity built when tranmission network is expanded
            N_STO_E[Nodes_Storages], (base_name=shorter["Storage_build_energy"], lower_bound=0) # Units: MWh; Storage energy technology built
            N_STO_P[Nodes_Storages], (base_name=shorter["Storage_capacity"], lower_bound=0) # Units: MW; Storage loading and discharging power capacity built
            N_SYNC[Nodes_SynCons], (base_name=shorter["SynCon_capacity"], lower_bound=0) # Units: MW; synchronous condenser capacity
            FLOW[Arcs,Hours], (base_name=shorter["Internodal_flow"]) # Units: MWh; Power flow between nodes in topology
            N_IC_EXP[Arcs], (base_name=shorter["Internodal_flow_expansion"], lower_bound=0) # Units: MW; Power flow expansion betweeen TxZones, and REZones to TxZones
            #
            EV_CHARGE[EV, Hours], (base_name=shorter["EV_charging"], lower_bound=0) # Units: MWh per time-interval; Electric vehicle charge for vehicle profile in set EV
            EV_DISCHARGE[EV, Hours], (base_name=shorter["EV_discharging"], lower_bound=0) # Units: MWh per time-interval; Electric vehicle dischargw for vehicle profile in set EV
            EV_L[EV, Hours], (base_name=shorter["EV_charge_level"], lower_bound=0) # Units: MWh at a given time-interval; Electric vehicle charging level for vehicle profile in set EV
            EV_PHEVFUEL[EV, Hours], (base_name=shorter["EV_PHEV_fuel_use"], lower_bound=0) #  Plug-in hybrid electric vehicle conventional fuel use
            EV_INF[EV, Hours], (base_name=shorter["EV_infeasible"], lower_bound=0) # Units: MWh per time-interval; Infeasibility term for Electric vehicle energy balance
            #
            H2_P2G[Nodes_P2G, Hours], (base_name=shorter["H2_power_to_gas"], lower_bound=0) # Units: MWh per time-interval; Power-to-gas energy conversion
            H2_G2P[Nodes_G2P, Hours], (base_name=shorter["H2_gas_to_power"], lower_bound=0) # Units: MWh per time-interval; Gas-to-power energy conversion
            H2_GS_L[Nodes_GasStorages, Hours], (base_name=shorter["H2_storage_level"], lower_bound=0) # Units: tonne-H2 at a given time-interval; Current gas storage level
            H2_GS_IN[Nodes_GasStorages, Hours], (base_name=shorter["H2_storage_inflow"], lower_bound=0) # Units: tonne-H2 at a given time-interval; Current gas storage input
            H2_GS_OUT[Nodes_GasStorages, Hours], (base_name=shorter["H2_storage_outflow"], lower_bound=0) # Units: tonne-H2 at a given time-interval; Current gas storage output
            N_P2G[Nodes_P2G], (base_name=shorter["H2_P2G_capacity"], lower_bound=0) # Units: MW; Power-to-gas capacity
            N_G2P[Nodes_G2P], (base_name=shorter["H2_G2P_capacity"], lower_bound=0) # Units: MW; Gas-to-power capacity
            N_GS[Nodes_GasStorages], (base_name=shorter["H2_storage_capacity"], lower_bound=0) # Units; tonne-H2 ; Gas storage (energy) capacity
            #
            HEAT_STO_L[BU,HP,Hours], (base_name=shorter["Heat_storage_level"], lower_bound=0) # Units: MWh at a given time-interval; Heating: storage level
            HEAT_HP_IN[BU,HP, Hours], (base_name=shorter["Heat_heat_pump_"], lower_bound=0)   # Units: MWh per time-interval; Heating: electricity demand from heat pump
            HEAT_INF[BU,HP, Hours], (base_name=shorter["Heat_infeasible"], lower_bound=0)  # Units: MWh per time-interval; Heating: Infeasibility term for Electric vehicle energy balance
    end)

# %% * ----------------------------------------------------------------------- *
#    ***** Objective function *****
#    * ----------------------------------------------------------------------- *

    @objective(m, Min, Z)

    @info "Objective function."
    @constraint(m, ObjectiveFunction,
cost_scaling*(sum(MarginalCost[n,t] * G[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours)

            + LoadIncreaseCost[n,t] * sum(G_UP[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours2)
            + LoadDecreaseCost[n,t] * sum(G_DO[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours)

            # + sum(CurtailmentCost * CU[(n,t),h] for (n,t) in Nodes_NonDispatch, h in Hours)

            # + sum(infeas_cost * G_INF[n,h] for n in Nodes, h in Hours)

            + sum(MarginalCost[n,sto] * (STO_OUT[(n,sto),h] + STO_IN[(n,sto),h])
                for (n,sto) in Nodes_Storages, h in Hours)

            + sum(MarginalCost[ev] * EV_DISCHARGE[ev, h]
                  + EvFuel[ev] * EV_PHEVFUEL[ev, h]
                  + infeas_cost * EV_INF[ev, h] for ev in EV, h in Hours)

            + sum(MarginalCost[n,p2g] * H2_P2G[(n,p2g),h] for (n,p2g) in Nodes_P2G, h in Hours)
            + sum(MarginalCost[n,g2p] * H2_G2P[(n,g2p),h] for (n,g2p) in Nodes_G2P, h in Hours)

            + sum(infeas_cost * HEAT_INF[bu,hp,h] for bu in BU, hp in HP, h in Hours)
            )

        +
cost_scaling*(sum(InvestmentCost[n,t] * N_TECH[(n,t)] for (n,t) in Nodes_Techs_New)
            + sum(InvestmentCostPower[n,sto] * N_STO_P[(n,sto)] for (n,sto) in Nodes_Storages_New)
            + sum(InvestmentCostEnergy[n,sto] * N_STO_E[(n,sto)] for (n,sto) in Nodes_Storages_New)

            + sum(InvestmentCostSynCon[n,syn] * N_SYNC[(n,syn)] for (n,syn) in Nodes_SynCons)

            + sum(FixedCost[n,t] * N_TECH[(n,t)] for (n,t) in Nodes_Techs)
            + sum(FixedCost[n,sto] * N_STO_P[(n,sto)] for (n,sto) in Nodes_Storages)

            + sum(InvestmentCostREZ_Exp[rez,txz] * N_REZ_EXP[rez] for (rez,txz) in keys(InvestmentCostREZ_Exp) )
            + sum(ArcWeight[from,to]*InvestmentCostTransExp[from,to] * N_IC_EXP[(from,to)] for (from,to) in Arcs)

            + sum(InvestmentCost[n,p2g] * N_P2G[(n,p2g)] for (n,p2g) in Nodes_P2G)
            + sum(InvestmentCost[n,g2p] * N_G2P[(n,g2p)] for (n,g2p) in Nodes_G2P)
            + sum(InvestmentCost[n,gs] * N_GS[(n,gs)] for (n,gs) in Nodes_GasStorages)

            + sum(FixedCost[n,p2g] * N_P2G[(n,p2g)] for (n,p2g) in Nodes_P2G)
            + sum(FixedCost[n,g2p] * N_G2P[(n,g2p)] for (n,g2p) in Nodes_G2P)
            + sum(FixedCost[n,gs] * N_GS[(n,gs)] for (n,gs) in Nodes_GasStorages)
           )
        == Z
    );

    next!(prog)
    println("\n")

# %% * ----------------------------------------------------------------------- *
#    ***** Energy balance and load levels *****
#    * ----------------------------------------------------------------------- *

## Toy example to understand index sets made from tuples that are used to index constraints:
# D = ["a","b"]; T = [3,4]; NR = [("a",1),("a",2),("b",1)]; mod = JuMP.Model();
# @variable(mod, x[NR,T])  # Note NR is a set of tuples.
# This constraint form is correct for using the constraint index 'for n in D'
# and the first tuple entry is then also `n`:
# @constraint(mod, [n=D], sum(x[(i,j),t] for (i,j) in NR if i == n for t in T) <= 1)
# but NOT this one (a new scope is created inside the sum for `n`):
# @constraint(mod, Test[n=D], sum(x[(n,j)] for (n,j) in NR) <= 1)
# and not this one (do not use commas in comprehension statements, use `for` to concatenate):
# @constraint(mod, [n=D], sum(x[(i,j),t] for (i,j) in NR if i == n, t in T) <= 1)

    # The implicit design assumption of this constraint is that DemandRegions are directly related to the Dispatch region level.
    # We may overcome this with nested sums from lower levels.
    # In particular, G_REZ is an aggregating variable for generation from Renewables and used directly in other constraints
    # while G_TxZ aggregates Dispatchable generation from Transmission regions.

    @info "Definition of REZone generation book-keeping variables"
    @constraint(m, REZoneGen[rez=REZones,h=Hours],
        G_REZ[rez,h] ==
            sum(G[(z,t),h] for (z,t) in Nodes_Techs if z == rez)
         +  sum(STO_OUT[(q,sto),h] - STO_IN[(q,sto),h] for (q,sto) in Nodes_Storages if q == rez)
    );

    @info "Definition of TxZone generation book-keeping variables"
    @constraint(m, TxZoneGen[zone=TxZones,h=Hours],
        G_TxZ[zone,h] ==
            sum(G[(z,t),h] for (z,t) in Nodes_Techs if z == zone)
         +  sum(G_REZ[rez,h] for (rez, z) in Nodes_Promotes if z == zone)
         +  sum(STO_OUT[(z,sto),h] - STO_IN[(z,sto),h] for (z,sto) in Nodes_Storages if z == zone)
         -  sum(FLOW[(from,to),h] for (from,to) in Arcs if from == zone)
    );

# # TODO : An explicit assumption in EnergyBalance is that H2_G2P and H2_P2G are
# indexed over DemandZones; we should add additional checks to verify this holds.

    # Energy balance at each demand node:
    @info "Energy balance at each demand node: supply equals or exceeds demand."
    @constraint(m, EnergyBalance[n=DemandZones,h=Hours],
      sum(G_TxZ[zone,h] for (zone, d) in Nodes_Demand if d == n)
        # + sum(EV_DISCHARGE[ev,h] for ev in EV)
        + sum(H2_G2P[(zone,g2p),h] for (zone,g2p) in Nodes_G2P if zone == n)
        - sum(NegOpDemand[zone,hour] for (zone,hour) in keys(NegOpDemand) if (zone == n && hour == h))
        >=
      loss_factor_tx*(Load[n,h]
                    + sum(H2_P2G[(zone,p2g),h] for (zone,p2g) in Nodes_P2G if zone == n))
      );
      # + sum(HEAT_HP_IN[bu,hp,h] for bu in BU for hp in HP)
      # + sum(EV_CHARGE[ev,h] for ev in EV)

    # Energy flow reflexive constraint:
    @info "Energy flow reflexive constraint."
    @constraint(m, FlowEnergyReflex[(from,to)=Arcs,h=Hours; from in Arcs_From],
        FLOW[(from,to),h] + FLOW[(to,from),h] == 0
    );

    # Energy flow bounds:
    @info "Energy flow bounds."
    @constraint(m, FlowEnergyUpperBound[(from,to)=Arcs,h=Hours],
        FLOW[(from,to),h] <= time_ratio * ( TransferCapacity[(from,to)] + N_IC_EXP[(from,to)] + CapAdd[:N_IC_EXP][(from,to)])
    );

    # @info "Transmission expansion upper bounds."
    # @constraint(m, FlowExpandUpperBound[(from,to)=Arcs],
    #     N_IC_EXP[(from,to)] + trunc(CapAdd[:N_IC_EXP][(from,to)]) <= ExpansionLimit_Tx[(from,to)]
    # );

    @info "Energy flow expansion symmetry."
    @constraint(m, FlowEnergySymmetry[(from,to)=Arcs; from in Arcs_From],
        N_IC_EXP[(from,to)] == N_IC_EXP[(to,from)]
    );

    # Generation level start - con2b_loadlevelstart
    @info "Generation level start."
    @constraint(m, GenLevelStart[(n,t)=Nodes_Dispatch],
               G[(n,t),Hours[1]] == G_UP[(n,t),Hours[1]]
    );

    # Generation level - con2a_loadlevel
    @info "Generation level."
    @constraint(m, GenLevelUpdate[(n,t)=Nodes_Dispatch,h=Hours2],
        G[(n,t),h] == G[(n,t),h-1] + G_UP[(n,t),h] - G_DO[(n,t),h]
    );

    next!(prog)
    println("\n")

    # Variable upper bound on dispatchable generation by capacity
    @info "Variable upper bound on dispatchable generation by capacity."
    @constraint(m, MaxGenerationDisp[(n,t)=Nodes_Dispatch,h=Hours],
        G[(n,t),h] <= CapacityDerating[n,t,h] * time_ratio * (N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)])
    );

    @info "Variable upper bound on non-dispatchable generation by capacity."
    @constraint(m, MaxGenerationNonDisp[(n,t)=Nodes_Avail_Techs,h=Hours],
        # G[(n,t),h] + CU[(n,t),h] == Availability[n,t,h] * time_ratio * N_TECH[(n,t)]
        G[(n,t),h] <= Availability[n,t,h] * time_ratio * (N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)])
    );

    # Maximum capacity allowed
    if !(isDictAllMissing(MaxCapacity))
        @info "Maximum capacity allowed."
        @constraint(m, MaxCapacityBound[(n,t)=Nodes_Techs; !(MaxCapacity[n,t] |> ismissing)],
            N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)] <= MaxCapacity[n,t]
        );
    end

    # Maximum generated energy allowed.
    if !(isDictAllMissing(MaxEnergy))
        @info "Maximum generated energy allowed."
        @constraint(m, MaxEnergyGenerated[(n,t)=Nodes_Techs; !(MaxEnergy[n,t] |> ismissing)],
            sum(G[(n,t),h] for h in Hours) <= MaxEnergy[n,t]
        );
    end

    # Renewable energy zone build limits.
    @info "Renewable energy zone build limits."
    @constraint(m, REZBuildLimits[rez=REZones,h=Hours],
        # sum(G[(z,t),h] for (z,t) in Nodes_Avail_Techs if z == rez)
        sum(G[(z,t),h] for (z,t) in Nodes_Techs if (z == rez && !occursin(r"Hydro_",t))) ## TODO: remove this hard-coding!
            <= time_ratio*(TotalBuildCap[rez] +  N_REZ_EXP[rez] + CapAdd[:N_REZ_EXP][rez]) # + N_REZ_EXP_TX[rez]
    );

    @info "Renewable energy zone expansion limits."
    @constraint(m, REZExpansionBound[rez=REZones],
        N_REZ_EXP[rez] + trunc(CapAdd[:N_REZ_EXP][rez]) <= ExpansionLimit_REZ[rez]
    );

    @info "Renewable energy zone expansion link to transmission expansion."
    @constraint(m, REZpairExpansion[rez=REZones],
        N_REZ_EXP_TX[rez] == 0  # NOTE: This constraint is now inactive.
            # sum(N_IC_EXP[(from,to)]/0.75 for (from,to) in Arcs_REZ_Pairs if ArcREZpairing[(from,to)] == rez)
    );

    next!(prog)
    println("\n")

# %% * ----------------------------------------------------------------------- *
#    ***** Operational conditions and constraints *****
#    * ----------------------------------------------------------------------- *

    # Ramping Limits
    
    if !(isDictAllMissing(MaxRampUpPerHour))
        @info "Maximum ramp up rates"
        @constraint(m, RampingUpLimits[(n,t)=Nodes_RampingTechs, h=Hours2],
            G_UP[(n,t),h] <= time_ratio * MaxRampUpPerHour[t]
        );
    end

    if !(isDictAllMissing(MaxRampDownPerHour))
        @info "Maximum ramp down rates"
        @constraint(m, RampingDownLimits[(n,t)=Nodes_RampingTechs, h=Hours2],
             G_DO[(n,t),h] <= time_ratio * MaxRampDownPerHour[t]
        );
    end

    # Operating reserves
    @info "Operating reserve margin"
    @constraint(m, OperatingReserve[dr=DemandRegions, h=Hours],
        sum( CapacityDerating[n,t,h] * time_ratio * (N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)]) - G[(n,t),h]
            for (n,t) in Nodes_Dispatch if node2DemReg[n] == dr)
      + sum( PeakContribution[n,t] * time_ratio * (N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)])
            for (n,t) in Nodes_Avail_Techs if node2DemReg[n] == dr)  # or replace `Nodes_Avail_Techs` with `setdiff(Nodes_Avail_Techs,Nodes_Dispatch)`
      + sum( sqrt(Efficiency[n,sto])*STO_L[(n,sto), h]
            for (n,sto) in Nodes_Storages if node2DemReg[n] == dr)
      >= time_ratio * OperatingReserve[dr]
    );

    # Minimum stable generation levels (for inflexible technologies)

    if !(isempty(keys(MinStableGen)))
        @info "Minimum stable generation levels"
        @constraint(m, MinStableGeneration[(n,t)=keys(MinStableGen), h=Hours; !(MinStableGen[n,t] |> ismissing)],
            G[(n,t),h] >= MinStableGen[n,t] * CapacityDerating[n,t,h] * time_ratio * (N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)])
        );
        ## Average generation version:
        # @constraint(m, MinStableGeneration[(n,t)=Nodes_Techs; !(MinStableGen[n,t] |> ismissing)],
            # sum(G[(n,t),h] for h in Hours) >= MinStableGen[n,t]*length(Hours) )
    end

#=
    @info "Carbon budget upper limit (annual)"
    @constraint(m, CarbonBudgetLimit,
        sum((CarbonContent[n,t]/Efficiency[n,t]) * G[(n,t),h]
                for (n,t) in Nodes_Techs, h in Hours)
        <= CarbonBudget
    );

    next!(prog)
    println("\n")
=#
# %% * ----------------------------------------------------------------------- *
#    ***** Quotas for renewable technologies *****
#    * ----------------------------------------------------------------------- *

    #  Minimum yearly renewables requirement (con5a_minRES)
    @info "Minimum yearly renewables requirement."
    @constraint(m, MinRES[n=DemandRegions],
        sum(
           sum(G[(z,t),h] for (z,t) in Nodes_Renew if z == zone)  # Any renewable TxZone-level tech.
         + sum(G_REZ[rez,h] for (rez, z) in Nodes_Promotes if z == zone)
        for (zone, d) in Nodes_Demand if node2DemReg[d] == n
        for h in Hours
        )
        >=
        (min_res_dict[n]/100)*(
            sum(
               sum(G[(z,t),h] for (z,t) in Nodes_Techs if z == zone)
             + sum(G_REZ[rez,h] for (rez, z) in Nodes_Promotes if z == zone)
            for (zone, d) in Nodes_Demand if node2DemReg[d] == n
            for h in Hours
         )
       )
    );
    # Note: if there is NO renewable energy associated to the DemandRegion, then this constraint will
    # act to zero out ALL generation associated to the DemandRegion.

            #  + sum(STO_OUT[(z,sto),h] - STO_IN[(z,sto),h] for (z,sto) in Nodes_Storages if z == zone)

            #  + sum(H2_G2P[(z,g2p),h] for (z,g2p) in Nodes_G2P if z == zone)
            #  - sum(H2_P2G[(z,p2g),h] for (z,p2g) in Nodes_P2G if z == zone)


    #  Minimum yearly renewables requirement for whole of system.
    @info "Minimum yearly renewables requirement for whole of system."
    @constraint(m, MinRESsystem,
        sum(
            sum(
                sum(G[(z,t),h] for (z,t) in Nodes_Renew if z == zone)  # Any renewable TxZone-level tech.
                + sum(G_REZ[rez,h] for (rez, z) in Nodes_Promotes if z == zone)
                for (zone, d) in Nodes_Demand if node2DemReg[d] == n
            for h in Hours)
        for n in DemandRegions)
            >=
            (MinimumRenewShare/100)*sum(
                sum(
                    sum(G[(z,t),h] for (z,t) in Nodes_Techs if z == zone)
                    + sum(G_REZ[rez,h] for (rez, z) in Nodes_Promotes if z == zone)
                    for (zone, d) in Nodes_Demand if node2DemReg[d] == n
                for h in Hours)
            for n in DemandRegions)
    );

# %% * ----------------------------------------------------------------------- *
#    ***** Storage constraints *****
#    * ----------------------------------------------------------------------- *


    # Storage Level Dynamics Initial Condition (con4a_stolev_start)
    @info "Storage level dynamics: initial condition."
    @constraint(m, StorageLevelStart[(n,sto)=Nodes_Storages],
        STO_L[(n,sto),Hours[1]]
         ==
        StartLevel[n,sto] * (N_STO_E[(n,sto)] + CapAdd[:N_STO_E][(n,sto)])
        +   sqrt(Efficiency[n,sto])*STO_IN[(n,sto), Hours[1]]
        - 1/sqrt(Efficiency[n,sto])*STO_OUT[(n,sto), Hours[1]]
    );

    # End level equal to initial level (con4j_ending)
    @info "Storage end level equal to initial level."
    @constraint(m, StorageLevelEnd[(n,sto)=Nodes_Storages],
        STO_L[(n,sto),Hours[end]]
         ==
        StartLevel[n,sto] * (N_STO_E[(n,sto)] + CapAdd[:N_STO_E][(n,sto)])
    );

    # @constraint(m, StorageBalanceFirstHour[(n,sto) in Nodes_Storages],
    #     STO_L[(n,sto),Hours[1]]
    #     ==
    #     STO_L[(n,sto),Hours[end]]
    #     +   sqrt(Efficiency[n,sto])*STO_IN[(n,sto), Hours[1]]
    #     - 1/sqrt(Efficiency[n,sto])*STO_OUT[(n,sto), Hours[1]]
    #     );

    # "Storage Level Dynamics
    @info "Storage level dynamics."
    @constraint(m, StorageBalance[(n,sto)=Nodes_Storages,h=Hours2],
        STO_L[(n,sto), h]
        ==
        STO_L[(n,sto), h-1]
        +   sqrt(Efficiency[n,sto])*STO_IN[(n,sto), h]
        - (1/sqrt(Efficiency[n,sto]))*STO_OUT[(n,sto), h]
    );

    # Storage Power Capacity (con4c_stolev_max)
    @info "Storage power capacity."
    @constraint(m, MaxLevelStorage[(n,sto)=Nodes_Storages,h=Hours],
        STO_L[(n,sto),h] <= (N_STO_E[(n,sto)] + CapAdd[:N_STO_E][(n,sto)])
    );

    # Storage maximum inflow (con4d_maxin_sto)
    @info "Storage maximum inflow."
    @constraint(m, MaxWithdrawStorage[(n,sto)=Nodes_Storages,h=Hours],
        STO_IN[(n,sto),h] <= time_ratio * (N_STO_P[(n,sto)] + CapAdd[:N_STO_P][(n,sto)])
    );

    # Storage maximum outflow by capacity (con4e_maxout_sto)
    @info "Storage generation outflow by capacity."
    @constraint(m, MaxGenerationStorage[(n,sto)=Nodes_Storages,h=Hours],
        STO_OUT[(n,sto),h] <= time_ratio * (N_STO_P[(n,sto)] + CapAdd[:N_STO_P][(n,sto)])
    );

    if !(isDictAllMissing(dtr.parameters[:MaxEnergy]))
        @info "Storage: maximum energy allowed."
        @constraint(m, MaxEnergyStorage[(n,sto)=Nodes_Storages; !(MaxEnergy[n,sto] |> ismissing)],
            N_STO_E[(n,sto)] <= MaxEnergy[n,sto]
        );
    end

    @info "Storage: maximum power allowed."
    @constraint(m, MaxPowerStorage[(n,sto)=Nodes_Storages; !(MaxCapacity[n,sto] |> ismissing)],
        N_STO_P[(n,sto)] + CapAdd[:N_STO_P][(n,sto)] <= MaxCapacity[n,sto]
    );

    # Maximum Energy to Power ratio for certain storage technologies (con4k_PHS_EtoP)
    @info "Storage: maximum energy-to-power ratio (use time)"
    @constraint(m, EnergyToPowerRatio[(n,sto)=Nodes_Storages; !(MaxEtoP_ratio[n,sto] |> ismissing)],
        N_STO_E[(n,sto)] + CapAdd[:N_STO_E][(n,sto)] <= MaxEtoP_ratio[n,sto] * (N_STO_P[(n,sto)] + CapAdd[:N_STO_P][(n,sto)])
    );

    # Maximum storage outflow - no more than level of last period (con4h_maxout_lev)
    @info "Storage: maximum outflow - no more than level of last period."
    @constraint(m, MaxOutflowStorage[(n,sto)=Nodes_Storages,h=Hours2],
         (1/sqrt(Efficiency[n,sto]))*STO_OUT[(n,sto),h] <= STO_L[(n,sto),h-1]
    );

    # Maximum storage inflow - no more than energy capacity minus level of last period (con4i_maxin_lev)
    @info "Storage: maximum inflow - no more than energy capacity minus level of last period"
    @constraint(m, MaxInflowStorage[(n,sto)=Nodes_Storages,h=Hours2],
         sqrt(Efficiency[n,sto])*STO_IN[(n,sto),h] <= N_STO_E[(n,sto)] + CapAdd[:N_STO_E][(n,sto)] - STO_L[(n,sto),h-1]
    );

    next!(prog)
    println("\n")

    return dtr
end


function build_ev_constraints(dtr::DieterModel)
# %% * ----------------------------------------------------------------------- *
#    ***** Electric vehicle constraints *****
#    * ----------------------------------------------------------------------- *

    # Sets:
    Hours = dtr.sets[:Hours]
    Hours2 = Hours[2:end]

    EV = dtr.sets[:ElectricVehicles]

    # Parameters:
    EvDemand = dtr.parameters[:AbsoluteEvDemand]
    EvPower = dtr.parameters[:AbsoluteEvPower]
    EvCapacity = dtr.parameters[:AbsoluteEvCapacity]
    EvFuel = dtr.parameters[:EvFuel]
    EvType = dtr.parameters[:EvType]

    # Variables:
    EV_CHARGE = dtr.model.obj_dict[:EV_CHARGE]
    EV_DISCHARGE = dtr.model.obj_dict[:EV_DISCHARGE]
    EV_L = dtr.model.obj_dict[:EV_L]
    EV_PHEVFUEL = dtr.model.obj_dict[:EV_PHEVFUEL]
    EV_INF = dtr.model.obj_dict[:EV_INF]

    @info "Electric Vehicles: electric charge allowable."
    @constraint(m, MaxWithdrawEV[ev=EV,h=Hours],
        EV_CHARGE[ev,h] <= EvPower[ev][h]
    );

    @info "Electric Vehicles: electric discharge allowable."
    @constraint(m, MaxGenerationEV[ev=EV,h=Hours],
        EV_DISCHARGE[ev,h] <= EvPower[ev][h]
    );

    @info "Electric Vehicles: charge level limited by capacity."
    @constraint(m, MaxLevelEV[ev=EV,h=Hours],
        EV_L[ev,h] <= EvCapacity[ev]
    );

    @info "Electric Vehicles: energy balance."
    @constraint(m, EVBalance[ev=EV,h=Hours2],
        EV_L[ev, h]
        ==
        EV_L[ev, h-1]
        + EV_CHARGE[ev, h]*sqrt(Efficiency[ev])
        - EV_DISCHARGE[ev, h]/sqrt(Efficiency[ev])
        - EvDemand[ev][h]
        + EvType[ev]*EV_PHEVFUEL[ev, h]
        + EV_INF[ev, h]
    );

    @info "Electric Vehicles: energy balance at first and last time-steps."
    @constraint(m, EVBalanceFirstHour[ev=EV],
        EV_L[ev,Hours[end]]    # or EV_L[ev,Hours[1]]
        ==
        EV_L[ev,Hours[1]]  # or EV_L[ev,Hours[end]]
        + EV_CHARGE[ev,Hours[end]]*sqrt(Efficiency[ev])
        - EV_DISCHARGE[ev,Hours[end]]/sqrt(Efficiency[ev])
        - EvDemand[ev][Hours[end]]
        + EvType[ev]*EV_PHEVFUEL[ev, Hours[end]]
        + EV_INF[ev, Hours[end]]
    );

    # next!(prog)
    # println("\n")

    return dtr
end


function build_h2_constraints(dtr::DieterModel)
# %% * ----------------------------------------------------------------------- *
#    ***** Hydrogen constraints *****
#    * ----------------------------------------------------------------------- *

    Timestep = dtr.settings[:timestep]
    periods = round(Int,hoursInYear*(2/Timestep))

    # Sets:
    Hours = dtr.sets[:Hours]
    Hours2 = Hours[2:end]

    HoursToMonths = dtr.parameters[:HoursToMonths]

    Nodes_P2G = dtr.sets[:Nodes_P2G]
    Nodes_P2G_PEM = dtr.sets[:Nodes_P2G_PEM]
    Nodes_P2G_AE = dtr.sets[:Nodes_P2G_AE]

    Nodes_G2P = dtr.sets[:Nodes_G2P]
    Nodes_GasStorages = dtr.sets[:Nodes_GasStorages]

    # Parameters:
    time_ratio = dtr.settings[:time_ratio]

    H2Conversion = dtr.parameters[:H2Conversion]  # Units: MWh / tonne-H2 for a given P2G tech.
    # H2Demand = coalesce((dtr.settings[:h2]*1e6)/hoursInYear,0)
    H2Demand = dtr.parameters[:H2Demand] # Units: MWh / year for a given Node and P2G tech.

    CapAdd = dtr.parameters[:CapAdd]

    Capacity_Factor_AE = dtr.settings[:capacity_factor_ae] # e.g 0.8 means must produce at least 80% of hourly quota each hour.

    # if dtr.settings[:H2_G2P_Allowed_flag] == true
    RecipH2_CF_UB = dtr.settings[:capacity_factor_recipH2_ub]
    # end

    # Variables:
    H2_P2G = dtr.model.obj_dict[:H2_P2G]
    H2_G2P = dtr.model.obj_dict[:H2_G2P]
    H2_GS_L = dtr.model.obj_dict[:H2_GS_L]
    H2_GS_IN = dtr.model.obj_dict[:H2_GS_IN]
    H2_GS_OUT = dtr.model.obj_dict[:H2_GS_OUT]
    N_P2G = dtr.model.obj_dict[:N_P2G]
    N_G2P = dtr.model.obj_dict[:N_G2P]
    N_GS = dtr.model.obj_dict[:N_GS]

    @info "Hydrogen: Minimum yearly lower bound on power-to-gas."
    @constraint(dtr.model, MinYearlyP2G[(n,p2g)=Nodes_P2G],
        sum(H2_P2G[(n,p2g),h] for h in Hours) >= H2Demand[n,p2g]
    );

    # @info "Hydrogen: Minimum monthly lower bound on power-to-gas for PEM tech."
    # @constraint(dtr.model, MinMonthlyP2G_PEM[(n,p2g)=Nodes_P2G_PEM, mth=1:12],
    #     sum(H2_P2G[(n,p2g),h] for h in Hours if HoursToMonths[h] == mth)
    #         >= (H2Demand[n,p2g]/periods)*sum(1.0 for h in Hours if HoursToMonths[h] == mth)
    # );

    @info "Hydrogen: Constant hourly lower bound on power-to-gas for AE tech."
    @constraint(dtr.model, MinHourlyP2G_AE[(n,p2g)=Nodes_P2G_AE, h=Hours],
        H2_P2G[(n,p2g),h] >= Capacity_Factor_AE*(H2Demand[n,p2g]/periods)
    );

    @info "Hydrogen: Variable upper bound on power-to-gas."
    @constraint(dtr.model, MaxP2G[(n,p2g)=Nodes_P2G,h=Hours],
        H2_P2G[(n,p2g),h] <= time_ratio * (N_P2G[(n,p2g)] + CapAdd[:N_P2G][(n,p2g)])
    );

    @info "Hydrogen: Variable upper bound on gas-to-power."
    @constraint(dtr.model, MaxG2P[(n,g2p)=Nodes_G2P,h=Hours],
        H2_G2P[(n,g2p),h] <= time_ratio * (N_G2P[(n,g2p)] + CapAdd[:N_G2P][(n,g2p)])
    );

    @constraint(dtr.model, MaxCapFactorH2_G2P[(n,g2p)=Nodes_G2P],
        sum(H2_G2P[(n,g2p),h] for h in Hours) <= RecipH2_CF_UB*length(Hours)*(N_G2P[(n,g2p)] + CapAdd[:N_G2P][(n,g2p)])
    );
              
    @info "Hydrogen: Variable upper bound on gas storage."
    @constraint(dtr.model, MaxLevelGasStorage[(n,gs)=Nodes_GasStorages,h=Hours],
        H2_GS_L[(n,gs),h] <= N_GS[(n,gs)]  # + CapAdd[:N_GS][(n,p2g)] )
    );

    @constraint(dtr.model, GasStorageIn[(n,gs)=Nodes_GasStorages,h=Hours],
        H2_GS_IN[(n,gs),h] == sum(H2_P2G[(n,p2g),h]/H2Conversion[n,p2g] for (n,p2g) in Nodes_P2G)
    );

    @constraint(dtr.model, GasStorageOut[(n,gs)=Nodes_GasStorages,h=Hours],
        H2_GS_OUT[(n,gs),h] == sum(H2_G2P[(n,g2p),h]/H2Conversion[n,g2p] for (n,g2p) in Nodes_G2P)
    );

#=
    @info "Hydrogen: energy balance."
    @constraint(dtr.model, H2Balance[h=Hours],
        sum(Efficiency[p2g]*H2_P2G[p2g,h] for p2g in Nodes_P2G)  # No sqrt on Efficiency, not a roundtrip
        + sum(H2_GS_IN[gs,h] for gs in Nodes_GasStorages)   ## Changed: H2_GS_OUT -> H2_GS_IN
        # + sum(H2_GS_OUT[gs,h] for gs in GasStorages)   ## Should this be changed: H2_GS_OUT -> H2_GS_IN ?
        ==
        sum((1/Efficiency[g2p])*H2_G2P[g2p,h] for g2p in Nodes_G2P)   # No sqrt on Efficiency, not a roundtrip
        + sum(H2_GS_OUT[gs,h] for gs in Nodes_GasStorages)    ## Changed: H2_GS_IN -> H2_GS_OUT
        # + sum(H2_GS_IN[gs,h] for gs in GasStorages)    ## Should this be changed: H2_GS_IN -> H2_GS_OUT ?
        + H2Demand
    );
=#
    @info "Hydrogen: gas storage mass balance."
    @constraint(dtr.model, GasStorageBalance[(n,gs)=Nodes_GasStorages,h=Hours2],
        H2_GS_L[(n,gs), h] == H2_GS_L[(n,gs), h-1] + H2_GS_IN[(n,gs), h] - H2_GS_OUT[(n,gs), h]
    );

    @info "Hydrogen: gas storage balance at first time-steps."
    @constraint(dtr.model, GasStorageBalanceFirstHours[(n,gs)=Nodes_GasStorages],
        H2_GS_L[(n,gs), Hours[1]] == StartLevel[n,gs] * N_GS[(n,gs)] + H2_GS_IN[(n,gs),Hours[1]] - H2_GS_OUT[(n,gs),Hours[1]]  # + CapAdd[:N_GS][(n,p2g)] )
    );

    # End level equal to initial level
    @info "Hydrogen: gas storage end level equal to initial level."
    @constraint(dtr.model, GasStorageLevelEnd[(n,gs)=Nodes_GasStorages],
        H2_GS_L[(n,gs), Hours[end]] == StartLevel[n,gs] * N_GS[(n,gs)]  # + CapAdd[:N_GS][(n,p2g)] )
    );

    # next!(prog)
    # println("\n")

    return dtr
end

function build_heat_load_constraints(dtr::DieterModel)
# %% * ----------------------------------------------------------------------- *
#    ***** Heat consumption constraints *****
#    * ----------------------------------------------------------------------- *
    # Sets:
    Hours = dtr.sets[:Hours]
    Hours2 = Hours[2:end]

    BU = dtr.sets[:BuildingType]  # Building archtypes
    HP = dtr.sets[:HeatingType]   # Heating combination type

    # Parameters:
    time_ratio = dtr.settings[:time_ratio]

    HeatConsumption = dtr.parameters[:HeatConsumption]  # Units: MWh per m^2 per hour
    HeatShare = dtr.parameters[:HeatShare]
    MaxLevel = dtr.parameters[:MaxLevel]
    HeatMaxPower = dtr.parameters[:HeatMaxPower]
    StaticEfficiency = dtr.parameters[:StaticEfficiency]
    CoP = dtr.parameters[:CoP]

    # Variables:
    HEAT_STO_L = dtr.model.obj_dict[:HEAT_STO_L]
    HEAT_HP_IN = dtr.model.obj_dict[:HEAT_HP_IN]
    HEAT_INF = dtr.model.obj_dict[:HEAT_INF]

    @constraint(m, HeatBalance[bu=BU, hp=HP, h=Hours2],
        HEAT_STO_L[bu,hp,h]
        ==
        StaticEfficiency[hp] * HEAT_STO_L[bu,hp,h-1]  # h-timestep replaced with h-1
        + CoP[bu,hp][h] * HEAT_HP_IN[bu,hp,h]
        - HeatConsumption[bu,hp][h]
        + HEAT_INF[bu,hp,h]
    );

    @constraint(m, HeatBalanceFirstHour[bu=BU, hp=HP],
        HEAT_STO_L[bu,hp,Hours[1]]
        ==
        StaticEfficiency[hp] * HEAT_STO_L[bu,hp,Hours[end]]
        + CoP[bu,hp][Hours[1]] * HEAT_HP_IN[bu,hp,Hours[1]]
        - HeatConsumption[bu,hp][Hours[1]]
        + HEAT_INF[bu,hp,Hours[1]]
    );

    @constraint(m, MaxHeatLevel[bu=BU, hp=HP, h=Hours],
        HEAT_STO_L[bu,hp,h] <= MaxLevel[bu]
    );

    @constraint(m, MaxHeatPower[bu=BU, hp=HP, h=Hours],
        HEAT_HP_IN[bu,hp,h] <= time_ratio * HeatMaxPower[bu,hp]
    );

    # next!(prog)

    return dtr
end

function solve_model!(dtr::DieterModel,solver)
    @info "Setting solver to $(solver)..."
    set_optimizer(dtr.model, solver)
    @info "Starting optimization..."
    JuMP.optimize!(dtr.model)

    return dtr
end

function generate_results!(dtr::DieterModel,varNamesDict::Dict{Symbol,Array{Symbol,1}})
    @info "Storing results"

    mo = dtr.model.obj_dict

    dtr.results = Dict(v[1] => convert_jump_container_to_df(mo[v[1]], dim_names=v[2]) for v in varNamesDict)
            # value_col=v[1])
    # dtr.results = [convert_jump_container_to_df(value.(model_dict[v[1]]), dim_names=v[2]) for v in varNamesDict]

    # if abs(sum(dtr.results[:G_INF][!, :Value]))  > (1e-5) ||
    #    ( !ismissing(dtr.settings[:ev]) &&
    #    abs(sum(dtr.results[:EV_INF][!, :Value])) > (1e-5) )
    #       @warn "Problem might be infeasable"
    # end

    return nothing
end
