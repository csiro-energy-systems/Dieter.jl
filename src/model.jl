## This file is part of Dieter.jl : Model definition

# const hoursInYear = 8760
# const infeas_cost = 1000

"""
Build the JuMP model describing the optimization problem, specifying the `solver` to use.
The `nthhour` parameter should be 0 for hourly steps, 1 for 2-hourly steps and
-1 for half-hourly steps. The data must match the time-steps of the `nthhour` parameter.
"""
function build_model!(dtr::DieterModel,
    solver::OptimizerFactory;
    nthhour::Int=0)

    dtr.settings[:nthhour] = nthhour

    periods = round(Int,hoursInYear*2.0^(-nthhour))
    Hours = Base.OneTo(periods)
    Hours2 = Hours[2:end]
    corr_factor = length(Hours)/hoursInYear

    H2Demand = coalesce((dtr.settings[:h2]*1e6)/hoursInYear,0)

    # Set definitions

    Technologies = dtr.sets[:Technologies]  # Generation technologies
    Storages = dtr.sets[:Storages]          # Storage technologies

    Renewables = dtr.sets[:Renewables]           # Renewable generation technologies
    Conventional = dtr.sets[:Conventional]       # Conventional generation technologies
    Dispatchable = dtr.sets[:Dispatchable]       # Dispatchable generation technologies
    NonDispatchable = dtr.sets[:NonDispatchable] # Non-dispatchable generation technologies

    Nodes = dtr.sets[:Nodes]

    DemandRegions = dtr.sets[:DemandRegions]

    min_res = Dict{String,Float64}()
    if isa(dtr.settings[:min_res],Number)
        for n in DemandRegions
            min_res[n] = dtr.settings[:min_res]
        end
    elseif isa(dtr.settings[:min_res],Dict)
        min_res = dtr.settings[:min_res]
    else
        @error "Setting `min_res` is not of type `Number` or a `Dict`."
    end
    dtr.settings[:min_res] = min_res
    @info "The `min_res` setting is now equal to "
    @info min_res

    # Mapping set definitions
    Arcs = dtr.sets[:Arcs]
    Arcs_From = dtr.sets[:Arcs_From]
    Nodes_Techs = dtr.sets[:Nodes_Techs]
    Nodes_Storages = dtr.sets[:Nodes_Storages]

    Nodes_Renew = dtr.sets[:Nodes_Renew]
    Nodes_Conven = dtr.sets[:Nodes_Conven]
    Nodes_Dispatch = dtr.sets[:Nodes_Dispatch]
    Nodes_NonDispatch = dtr.sets[:Nodes_NonDispatch]

    Nodes_Types = dtr.sets[:Nodes_Types]
    Nodes_Promotes = dtr.sets[:Nodes_Promotes]

    Nodes_Demand = dtr.sets[:Nodes_Demand]

    TxZones = dtr.sets[:TxZones]
    REZones = dtr.sets[:REZones]
    # TxZones = [x[1] for x in filter(x -> (x[2]=="TxZone"),Nodes_Types)]
    # REZones = [x[1] for x in filter(x -> (x[2]=="REZone"),Nodes_Types)]

    ## Electric Vehicles
    EV = dtr.sets[:ElectricVehicles]

    ## Building heating
    BU = dtr.sets[:BuildingType]  # Building archtypes
    HP = dtr.sets[:HeatingType]   # Heating combination type

    ## Power-to-gas with Hydrogen
    P2G = dtr.sets[:P2G]
    G2P = dtr.sets[:G2P]
    GasStorages = dtr.sets[:GasStorages]

    # Parameter definitions

    EvDemand = dtr.parameters[:AbsoluteEvDemand]
    EvPower = dtr.parameters[:AbsoluteEvPower]
    EvCapacity = dtr.parameters[:AbsoluteEvCapacity]
    EvFuel = dtr.parameters[:EvFuel]
    EvType = dtr.parameters[:EvType]

    HeatConsumption = dtr.parameters[:HeatConsumption]  # Units: MWh per m^2 per hour
    HeatShare = dtr.parameters[:HeatShare]
    MaxLevel = dtr.parameters[:MaxLevel]
    HeatMaxPower = dtr.parameters[:HeatMaxPower]
    StaticEfficiency = dtr.parameters[:StaticEfficiency]
    CoP = dtr.parameters[:CoP]


    MarginalCost = dtr.parameters[:MarginalCost]
    InvestmentCost = dtr.parameters[:InvestmentCost]
    InvestmentCostPower = dtr.parameters[:InvestmentCostPower]
    InvestmentCostEnergy = dtr.parameters[:InvestmentCostEnergy]
    FixedCost = dtr.parameters[:FixedCost]
    Availability = dtr.parameters[:Availability]
    MaxCapacity = dtr.parameters[:MaxCapacity]
    MaxEnergy = dtr.parameters[:MaxEnergy]
    MaxEtoPratio = dtr.parameters[:MaxEnergyToPowerRatio]
    Efficiency = dtr.parameters[:Efficiency]
    StartLevel = dtr.parameters[:StartLevel]
    Load = dtr.parameters[:Load]

    CurtailmentCost = dtr.settings[:cu_cost]

    LoadIncreaseCost = dtr.parameters[:LoadIncreaseCost]
    LoadDecreaseCost = dtr.parameters[:LoadDecreaseCost]

    @info "Start of model building:"

    prog = Progress(7, dt=0.01, desc="Building Model...         \n", barlen=30)

    m = Model(solver)

    @info "Variable definitions."
    @variables(m, begin
        Z, (base_name="Total_cost_objective", lower_bound=0)
        G[Nodes_Techs, Hours], (base_name="Generation_level", lower_bound=0) # Units: MWh; Generation level - all generation tech.
        G_UP[Nodes_Dispatch, Hours] , (base_name="Generation_upshift", lower_bound=0) # Units: MWh; Generation level change up
        G_DO[Nodes_Dispatch, Hours], (base_name="Generation_downshift", lower_bound=0) # Units: MWh Generation level change down
        G_INF[Nodes, Hours], (base_name="Generation_infeasible", lower_bound=0) # Units: MWh; Infeasibility term for Energy Balance
        G_REZ[REZones,Hours], (base_name="Generation_renewable_zones", lower_bound=0) # Units: MWh; Generation level - renewable energy zone tech. & stor.
        G_TxZ[TxZones,Hours], (base_name="Generation_transmission_zones", lower_bound=0) # Units: MWh; Generation level - transmission zone tech. & stor.
        # G_RES[Nodes_Renew, h in HOURS], (base_name="Generation_renewable", lower_bound=0) # Units: MWh; Generation level - renewable gen. tech.
        CU[Nodes_NonDispatch, Hours], (base_name="Curtailment_renewables", lower_bound=0) # Units: MWh; Non-dispatchable curtailment
        STO_IN[Nodes_Storages, Hours], (base_name="Storage_inflow", lower_bound=0) # Units: MWh per h; Storage energy inflow
        STO_OUT[Nodes_Storages, Hours], (base_name="Storage_outflow", lower_bound=0) # Units: MWh per h; Storage energy outflow
        STO_L[Nodes_Storages, Hours], (base_name="Storage_level", lower_bound=0) # Units: MWh; Storage energy level
        N_TECH[Nodes_Techs], (base_name="Technology_capacity", lower_bound=0) # Units: MW; Technology capacity built
        # N_RES[Nodes_Renew], (base_name="Renewable_capacity", lower_bound=0) # Units: MW; Renewable technology capacity built
        N_STO_E[Nodes_Storages], (base_name="Storage_build_energy", lower_bound=0) # Units: MWh; Storage technology built - Energy
        N_STO_P[Nodes_Storages], (base_name="Storage_capacity", lower_bound=0) # Units: MW; Storage loading and discharging capacity built
        FLOW[Arcs,Hours], (base_name="Internodal_flow") # Units: MWh; Power flow between nodes in topology
        EV_CHARGE[EV, Hours], (base_name="EV_charging", lower_bound=0) # Units: MWh; Electric vehicle charge for vehicle profile in set EV
        EV_DISCHARGE[EV, Hours], (base_name="EV_discharging", lower_bound=0) # Units: MWh; Electric vehicle dischargw for vehicle profile in set EV
        EV_L[EV, Hours], (base_name="EV_charge_level", lower_bound=0) # Units: MWh; Electric vehicle charging level for vehicle profile in set EV
        EV_PHEVFUEL[EV, Hours], (base_name="EV_PHEV_fuel_use", lower_bound=0) #  Plug in hybrid electric vehicle conventional fuel use
        EV_INF[EV, Hours], (base_name="EV_infeasible", lower_bound=0) # Units: MWh; Infeasibility term for Electric vehicle energy balance
        H2_P2G[P2G, Hours], (base_name="H2_power_to_gas", lower_bound=0) # Units: MWh; Power-to-gas conversion
        H2_G2P[G2P, Hours], (base_name="H2_gas_to_power", lower_bound=0) # Units: MWh; Gas-to-power conversion
        H2_GS_L[GasStorages, Hours], (base_name="H2_storage_level", lower_bound=0) # Units: Current gas storage level
        H2_GS_IN[GasStorages, Hours], (base_name="H2_storage_inflow", lower_bound=0)
        H2_GS_OUT[GasStorages, Hours], (base_name="H2_storage_outflow", lower_bound=0)
        N_P2G[P2G], (base_name="H2_P2G_capacity", lower_bound=0)
        N_G2P[G2P], (base_name="H2_G2P_capacity", lower_bound=0)
        N_GS[GasStorages], (base_name="H2_storage_capacity", lower_bound=0)
        HEAT_STO_L[BU,HP,Hours], (base_name="Heat_storage_level", lower_bound=0) # Units: MWh; Heating: storage level
        HEAT_HP_IN[BU,HP, Hours], (base_name="Heat_heat_pump_", lower_bound=0)   # Units: MWh; Heating: electricity demand from heat pump
        HEAT_INF[BU,HP, Hours], (base_name="Heat_infeasible", lower_bound=0)  # Units: MWh; Heating: Infeasibility term for Electric vehicle energy balance
    end)

    # @variable(m, G[Nodes_Techs, Hours] >= 0)
    # @variable(m, G_UP[Nodes_Dispatch, Hours] >= 0)
    # @variable(m, G_DO[Nodes_Dispatch, Hours] >= 0)
    #
    # @variable(m, CU[Nodes_NonDispatch, Hours] >= 0)
    # @variable(m, STO_IN[Nodes_Storages, Hours] >= 0)
    # @variable(m, STO_OUT[Nodes_Storages, Hours] >= 0)
    # @variable(m, STO_L[Nodes_Storages, Hours] >= 0)
    # @variable(m, N_TECH[Nodes_Techs] >= 0)
    # @variable(m, N_STO_E[Nodes_Storages] >= 0)
    # @variable(m, N_STO_P[Nodes_Storages] >= 0)
    # @variable(m, G_INF[Nodes,Hours] >= 0)
    #
    # @variable(m, FLOW[Arcs,Hours])
    #
    # @variable(m, EV_CHARGE[EV, Hours] >= 0)
    # @variable(m, EV_DISCHARGE[EV, Hours] >= 0)
    # @variable(m, EV_L[EV, Hours] >= 0)
    # @variable(m, EV_PHEVFUEL[EV, Hours] >= 0)
    # @variable(m, EV_INF[EV, Hours] >= 0)
    #
    # @variable(m, H2_P2G[P2G, Hours] >= 0)
    # @variable(m, H2_G2P[G2P, Hours] >= 0)
    # @variable(m, H2_GS_L[GasStorages, Hours] >= 0)
    # @variable(m, H2_GS_IN[GasStorages, Hours] >= 0)
    # @variable(m, H2_GS_OUT[GasStorages, Hours] >= 0)
    # @variable(m, N_P2G[P2G] >= 0)
    # @variable(m, N_G2P[G2P] >= 0)
    # @variable(m, N_GS[GasStorages] >= 0)
    #
    # @variable(m, HEAT_STO_L[BU,HP,Hours] >= 0)
    # @variable(m, HEAT_HP_IN[BU,HP, Hours] >= 0)
    # @variable(m, HEAT_INF[BU,HP, Hours] >= 0)


# %%   * --------------------------------------------------------------------- *
#    ***** Objective function *****
#    * ----------------------------------------------------------------------- *

    @objective(m, Min, Z)

    @info "Objective function."
    @constraint(m, ObjectiveFunction,
        Z ==
        sum(MarginalCost[n,t] * G[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours)

        + sum(LoadIncreaseCost[n,t] * G_UP[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours2)
        + sum(LoadDecreaseCost[n,t] * G_DO[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours)

        + sum(CurtailmentCost * CU[(n,t),h] for (n,t) in Nodes_NonDispatch, h in Hours)

        + sum(infeas_cost * G_INF[n,h] for n in Nodes, h in Hours)

        + sum(MarginalCost[n,sto] * (STO_OUT[(n,sto),h] + STO_IN[(n,sto),h])
            for (n,sto) in Nodes_Storages, h in Hours)

        + sum(MarginalCost[ev] * EV_DISCHARGE[ev, h]
              + EvFuel[ev] * EV_PHEVFUEL[ev, h]
              + infeas_cost * EV_INF[ev, h] for ev in EV, h in Hours)

        + sum(MarginalCost[g2p] * H2_G2P[g2p,h] for g2p in G2P, h in Hours)

        + sum(infeas_cost * HEAT_INF[bu,hp,h] for bu in BU, hp in HP, h in Hours)


        + corr_factor *
            (sum(InvestmentCost[n,t] * N_TECH[(n,t)] for (n,t) in Nodes_Techs)
            + sum(InvestmentCostPower[n,sto] * N_STO_P[(n,sto)] for (n,sto) in Nodes_Storages)
            + sum(InvestmentCostEnergy[n,sto] * N_STO_E[(n,sto)] for (n,sto) in Nodes_Storages)

            + sum(FixedCost[n,t] * N_TECH[(n,t)] for (n,t) in Nodes_Techs)
            + sum(FixedCost[n,sto] * 0.5*(N_STO_P[(n,sto)] + N_STO_E[(n,sto)]) for (n,sto) in Nodes_Storages)

            + sum(InvestmentCost[p2g] * N_P2G[p2g] for p2g in P2G)
            + sum(InvestmentCost[g2p] * N_G2P[g2p] for g2p in G2P)
            + sum(InvestmentCost[gs] * N_GS[gs] for gs in GasStorages)

            + sum(FixedCost[p2g] * N_P2G[p2g] for p2g in P2G)
            + sum(FixedCost[g2p] * N_G2P[g2p] for g2p in G2P)
            + sum(FixedCost[gs] * N_GS[gs] for gs in GasStorages)
           )
    );

    next!(prog)
    @info "\n"

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
    # In particular, G for NonDispatchable aggregates generation from Renewables via other constraints
    # while G for Dispatchable aggregates from Transmission regions.

    @info "Definition of REZone generation book-keeping variables"
    @constraint(m, REZoneGen[rez=REZones,h=Hours],
        G_REZ[rez,h] ==
            sum(G[(z,t),h] for (z,t) in Nodes_Techs if z == rez)
         +  sum(STO_OUT[(q,sto),h] for (q,sto) in Nodes_Storages if q == rez)
         -  sum( STO_IN[(q,sto),h] for (q,sto) in Nodes_Storages if q == rez)
    );

    @info "Definition of TxZone generation book-keeping variables"
    @constraint(m, TxZoneGen[zone=TxZones,h=Hours],
        G_TxZ[zone,h] ==
           sum(G[(z,t),h] for (z,t) in Nodes_Techs if z == zone)
        + sum(G_REZ[rez,h] for (rez, z) in Nodes_Promotes if z == zone)
        + sum(STO_OUT[(z,sto),h] for (z,sto) in Nodes_Storages if z == zone)
        - sum( STO_IN[(z,sto),h] for (z,sto) in Nodes_Storages if z == zone)
        - sum(FLOW[(from,to),h] for (from,to) in Arcs if from == zone)
    );


    # Energy balance at each demand node:
    @info "Energy balance at each demand node."
    @constraint(m, EnergyBalance[n=DemandRegions,h=Hours],
      # sum(G[(p,t),h] for (p,t) in Nodes_Techs if p == n)

      sum(G_TxZ[zone,h] for (zone, d) in Nodes_Demand if d == n)
        # + sum(STO_OUT[(q,sto),h] for (q,sto) in Nodes_Storages if q == n)
        + sum(EV_DISCHARGE[ev,h] for ev in EV)
        + sum(H2_G2P[g2p,h] for g2p in G2P)
        ==
      Load[n,h]
        # = sum(Load[a,h] for (a,b) in Nodes_Demand if b == n)
        # + sum(STO_IN[(q,sto),h] for (q,sto) in Nodes_Storages if q == n)
        + sum(EV_CHARGE[ev,h] for ev in EV)
        + sum(H2_P2G[p2g,h] for p2g in P2G)
        + sum(HEAT_HP_IN[bu,hp,h] for bu in BU, hp in HP)

    );

    # Energy flow reflexive constraint:
    @info "Energy flow reflexive constraint."
    @constraint(m, FlowEnergyReflex[(from,to)=Arcs,h=Hours; from in Arcs_From],
        FLOW[(from,to),h] + FLOW[(to,from),h] == 0
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

    # Variable upper bound on dispatchable generation by capacity
    @info "Variable upper bound on non-dispatchable generation by capacity."
    @constraint(m, MaxGenerationDisp[(n,t)=Nodes_Dispatch,h=Hours],
        G[(n,t),h] <= N_TECH[(n,t)]
    );

    @info "Variable upper bound on non-dispatchable generation by capacity."
    @constraint(m, MaxGenerationNonDisp[(n,t)=Nodes_NonDispatch,h=Hours],
        G[(n,t),h] + CU[(n,t),h] == Availability[n,t,h] * N_TECH[(n,t)]
    );

    # Maximum capacity allowed
    @info "Maximum capacity allowed."
    @constraint(m, MaxCapacityBound[(n,t)=Nodes_Techs; !(MaxCapacity[n,t] |> ismissing)],
        N_TECH[(n,t)] <= MaxCapacity[n,t]
    );

    # Maximum generated energy allowed.
    @info "Maximum generated energy allowed."
    @constraint(m, MaxEnergyGenerated[(n,t)=Nodes_Techs; !(MaxEnergy[n,t] |> ismissing)],
        sum(G[(n,t),h] for h in Hours) <= corr_factor * MaxEnergy[n,t]
    );

    next!(prog)

# %% * ----------------------------------------------------------------------- *
#    ***** Quotas for renewable technologies *****
#    * ----------------------------------------------------------------------- *

    #  Minimum yearly renewables requirement (con5a_minRES)
    @info "Minimum yearly renewables requirement."
    @constraint(m, MinRES[n=DemandRegions],
        sum(G[(p,res),h] for (p,res) in Nodes_Renew if p == n for h in Hours)
        + sum(STO_OUT[(q,sto),h] - STO_IN[(q,sto),h] for (q,sto) in Nodes_Storages if q == n for h in Hours)
        + sum(H2_G2P[g2p,h] - H2_P2G[p2g,h] for p2g in P2G, g2p in G2P, h in Hours)
        >=
        (min_res[n]/100)*(
        sum(G[(p,t),h] for (p,t) in Nodes_Techs if p == n for h in Hours)
        + sum(STO_OUT[(q,sto),h] - STO_IN[(q,sto),h] for (q,sto) in Nodes_Storages if q == n for h in Hours)
        + sum(H2_G2P[g2p,h] - H2_P2G[p2g,h] for p2g in P2G, g2p in G2P, h in Hours)
        )
    );

# %% * ----------------------------------------------------------------------- *
#    ***** Storage constraints *****
#    * ----------------------------------------------------------------------- *


    # Storage Level Dynamics Initial Condition (con4a_stolev_start)
    @info "Storage level dynamics: initial condition."
    @constraint(m, StorageLevelStart[(n,sto)=Nodes_Storages],
        STO_L[(n,sto),Hours[1]]
         ==
        StartLevel[n,sto] * N_STO_E[(n,sto)]
        +   sqrt(Efficiency[n,sto])*STO_IN[(n,sto), Hours[1]]
        - 1/sqrt(Efficiency[n,sto])*STO_OUT[(n,sto), Hours[1]]
    );

    # End level equal to initial level (con4j_ending)
    @info "Storage end level equal to initial level."
    @constraint(m, StorageLevelEnd[(n,sto)=Nodes_Storages],
        STO_L[(n,sto),Hours[end]]
         ==
        StartLevel[n,sto] * N_STO_E[(n,sto)]
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
        STO_L[(n,sto),h] <= N_STO_E[(n,sto)]
    );

    # Storage maximum inflow (con4d_maxin_sto)
    @info "Storage maximum inflow."
    @constraint(m, MaxWithdrawStorage[(n,sto)=Nodes_Storages,h=Hours],
        STO_IN[(n,sto),h] <= N_STO_P[(n,sto)]
    );

    # Storage maximum outflow by capacity (con4e_maxout_sto)
    @info "Storage generation outflow by capacity."
    @constraint(m, MaxGenerationStorage[(n,sto)=Nodes_Storages,h=Hours],
        STO_OUT[(n,sto),h] <= N_STO_P[(n,sto)]
    );

    @info "Storage: maximum energy allowed."
    @constraint(m, MaxEnergyStorage[(n,sto)=Nodes_Storages; !(MaxEnergy[n,sto] |> ismissing)],
        N_STO_E[(n,sto)] <= MaxEnergy[n,sto]
    );

    @info "Storage: maximum power allowed."
    @constraint(m, MaxPowerStorage[(n,sto)=Nodes_Storages; !(MaxCapacity[n,sto] |> ismissing)],
        N_STO_P[(n,sto)] <= MaxCapacity[n,sto]
    );

    # Maximum Energy to Power ratio for certain storage technologies (con4k_PHS_EtoP)
    @info "Storage: maximum energy-to-power ratio (use time)"
    @constraint(m, EnergyToPowerRatio[(n,sto)=Nodes_Storages; !(MaxEtoPratio[n,sto] |> ismissing)],
        N_STO_E[(n,sto)] <= MaxEtoPratio[n,sto]*N_STO_P[(n,sto)]
    );

    # Maximum storage outflow - no more than level of last period (con4h_maxout_lev)
    @info "Storage: maximum outflow - no more than level of last period."
    @constraint(m, MaxOutflowStorage[(n,sto)=Nodes_Storages,h=Hours2],
         (1/sqrt(Efficiency[n,sto]))*STO_OUT[(n,sto),h] <= STO_L[(n,sto),h-1]
    );

    # Maximum storage inflow - no more than energy capacity minus level of last period (con4i_maxin_lev)
    @info "Storage: maximum inflow - no more than energy capacity minus level of last period"
    @constraint(m, MaxInflowStorage[(n,sto)=Nodes_Storages,h=Hours2],
         sqrt(Efficiency[n,sto])*STO_IN[(n,sto),h] <= N_STO_E[(n,sto)] - STO_L[(n,sto),h-1]
    );

    next!(prog)

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

    next!(prog)

    @constraint(m, MaxP2G[p2g=P2G,h=Hours],
        H2_P2G[p2g,h] <= N_P2G[p2g]
    );

    @constraint(m, MaxG2P[g2p=G2P,h=Hours],
        H2_G2P[g2p,h] <= N_G2P[g2p]
    );

    @constraint(m, MaxLevelGasstorage[gs=GasStorages,h=Hours],
        H2_GS_L[gs,h] <= N_GS[gs]
    );

    @constraint(m, H2Balance[h=Hours],
        sum(sqrt(Efficiency[p2g])*H2_P2G[p2g,h] for p2g in P2G)  # Added sqrt to Efficiency
        + sum(H2_GS_IN[gs,h] for gs in GasStorages)   ## Changed: H2_GS_OUT -> H2_GS_IN
        # + sum(H2_GS_OUT[gs,h] for gs in GasStorages)   ## Should this be changed: H2_GS_OUT -> H2_GS_IN ?
        ==
        sum((1/sqrt(Efficiency[g2p]))*H2_G2P[g2p,h] for g2p in G2P)  # Added sqrt to Efficiency
        + sum(H2_GS_OUT[gs,h] for gs in GasStorages)    ## Changed: H2_GS_IN -> H2_GS_OUT
        # + sum(H2_GS_IN[gs,h] for gs in GasStorages)    ## Should this be changed: H2_GS_IN -> H2_GS_OUT ?
        + H2Demand
    );

    @constraint(m, GasstorageBalance[gs=GasStorages,h=Hours2],
        H2_GS_L[gs, h]
        ==
        H2_GS_L[gs, h-1]
        + H2_GS_IN[gs, h]
        - H2_GS_OUT[gs, h]
    );

    @constraint(m, GasstorageBalanceFirstHours[gs=GasStorages],
        H2_GS_L[gs, Hours[1]]
        ==
        H2_GS_L[gs, Hours[end]]
        + H2_GS_IN[gs,Hours[1]]
        - H2_GS_OUT[gs,Hours[1]]
    );

    next!(prog)

    @constraint(m, HeatBalance[bu=BU, hp=HP, h=Hours2],
        HEAT_STO_L[bu,hp,h]
        ==
        StaticEfficiency[hp] * HEAT_STO_L[bu,hp,h-1]  # h-nthhour replaced with h-1
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
        HEAT_HP_IN[bu,hp,h] <= HeatMaxPower[bu,hp]
    );

    next!(prog)

    dtr.model = m

    return dtr
end

function solve_model!(dtr::DieterModel)
    @info "Starting optimization..."
    optimize!(dtr.model)

    return dtr
end

function generate_results!(dtr::DieterModel)
    @info "Storing results"

    vars = [
        :Z => [],
        :G => [:Nodes_Techs, :Hours],
        :G_UP => [:Nodes_Dispatch, :Hours],
        :G_DO => [:Nodes_Dispatch, :Hours],
        :G_INF => [:Nodes,:Hours],
        :CU => [:Nodes_NonDispatch, :Hours],
        :STO_IN => [:Nodes_Storages, :Hours],
        :STO_OUT => [:Nodes_Storages, :Hours],
        :STO_L => [:Nodes_Storages, :Hours],
        :N_TECH => [:Nodes_Techs],
        :N_STO_E => [:Nodes_Storages],
        :N_STO_P => [:Nodes_Storages],
        :FLOW => [:Arcs,:Hours],

        :EV_CHARGE => [:EV, :Hours],
        :EV_DISCHARGE => [:EV, :Hours],
        :EV_L => [:EV, :Hours],
        :EV_PHEVFUEL => [:EV, :Hours],
        :EV_INF => [:EV, :Hours],

        :H2_P2G => [:P2G, :Hours],
        :H2_G2P => [:G2P, :Hours],
        :H2_GS_L => [:GasStorages, :Hours],
        :H2_GS_IN => [:GasStorages, :Hours],
        :H2_GS_OUT => [:GasStorages, :Hours],
        :N_P2G => [:P2G],
        :N_G2P => [:G2P],
        :N_GS => [:GasStorages],

        :HEAT_STO_L => [:BU, :HP, :Hours],
        :HEAT_HP_IN => [:BU, :HP, :Hours],
        :HEAT_INF => [:BU, :HP, :Hours]
    ]

    m = dtr.model

    model_dict = m.obj_dict

    dtr.results = Dict(v[1] => convert_jump_container_to_df(model_dict[v[1]], dim_names=v[2]) for v in vars)
    # dtr.results = [convert_jump_container_to_df(value.(model_dict[v[1]]), dim_names=v[2]) for v in vars]

    if abs(sum(dtr.results[:G_INF][!, :Hours]))  > (1e-5) ||
       abs(sum(dtr.results[:EV_INF][!, :Hours])) > (1e-5)
          @warn "Problem might be infeasable"
    end

    return nothing
end
