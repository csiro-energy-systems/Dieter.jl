## This file is part of Dieter.jl : Model definition

const hoursInYear = 8760

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
    Hours2 = collect(2:length(Hours))
    corr_factor = length(Hours)/hoursInYear

    H2Demand = coalesce((dtr.settings[:h2]*1e6)/hoursInYear,0)
    min_res = dtr.settings[:min_res]

    # Set definitions

    Technologies = dtr.sets[:Technologies]  # Generation technologies
    Storages = dtr.sets[:Storages]          # Storage technologies

    Renewables = dtr.sets[:Renewables]           # Renewable generation technologies
    Conventional = dtr.sets[:Conventional]       # Conventional generation technologies
    Dispatchable = dtr.sets[:Dispatchable]       # Dispatchable generation technologies
    NonDispatchable = dtr.sets[:NonDispatchable] # Non-dispatchable generation technologies

    Nodes = dtr.sets[:Nodes]

    DemandNodes = ["DE"]  # Fix this.

    # Mapping set definitions
    Arcs = dtr.sets[:Arcs]
    Nodes_Techs = dtr.sets[:Nodes_Techs]
    Nodes_Storages = dtr.sets[:Nodes_Storages]

    Nodes_Renew = dtr.sets[:Nodes_Renew]
    Nodes_Conven = dtr.sets[:Nodes_Conven]
    Nodes_Dispatch = dtr.sets[:Nodes_Dispatch]
    Nodes_NonDispatch = dtr.sets[:Nodes_NonDispatch]

    Nodes_Levels = dtr.sets[:Nodes_Levels]
    Nodes_Promotes = dtr.sets[:Nodes_Promotes]

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

    HeatConsumption = dtr.parameters[:HeatConsumption]
    HeatShare = dtr.parameters[:HeatShare]
    MaxLevel = dtr.parameters[:MaxLevel]
    HeatMaxPower = dtr.parameters[:HeatMaxPower]
    StaticEfficiency = dtr.parameters[:StaticEfficiency]
    COP = dtr.parameters[:COP]


    MarginalCost = dtr.parameters[:MarginalCost]
    InvestmentCost = dtr.parameters[:InvestmentCost]
    InvestmentCostPower = dtr.parameters[:InvestmentCostPower]
    InvestmentCostEnergy = dtr.parameters[:InvestmentCostEnergy]
    FixedCost = dtr.parameters[:FixedCost]
    Availability = dtr.parameters[:Availability]
    MaxInstallable = dtr.parameters[:MaxInstallable]
    MaxEnergy = dtr.parameters[:MaxEnergy]
    MaxPower = dtr.parameters[:MaxPower]
    Efficiency = dtr.parameters[:Efficiency]
    Load = dtr.parameters[:Load]

    CurtailmentCost = dtr.settings[:cu_cost]

    LoadIncreaseCost = dtr.parameters[:LoadIncreaseCost]
    LoadDecreaseCost = dtr.parameters[:LoadDecreaseCost]

    @info "Building Model"

    m = Model(solver)

    @variable(m, G[Nodes_Techs, Hours] >= 0)
    @variable(m, G_UP[Nodes_Dispatch, Hours] >= 0)
    @variable(m, G_DO[Nodes_Dispatch, Hours] >= 0)

    @variable(m, CU[Nodes_NonDispatch, Hours] >= 0)
    @variable(m, STO_IN[Nodes_Storages, Hours] >= 0)
    @variable(m, STO_OUT[Nodes_Storages, Hours] >= 0)
    @variable(m, STO_L[Nodes_Storages, Hours] >= 0)
    @variable(m, N[Nodes_Techs] >= 0)
    @variable(m, N_STO_E[Nodes_Storages] >= 0)
    @variable(m, N_STO_P[Nodes_Storages] >= 0)
    @variable(m, G_INF[Nodes,Hours] >= 0)

    @variable(m, EV_CHARGE[EV, Hours] >= 0)
    @variable(m, EV_DISCHARGE[EV, Hours] >= 0)
    @variable(m, EV_L[EV, Hours] >= 0)
    @variable(m, EV_PHEVFUEL[EV, Hours] >= 0)
    @variable(m, EV_INF[EV, Hours] >= 0)

    @variable(m, H2_P2G[P2G, Hours] >= 0)
    @variable(m, H2_G2P[G2P, Hours] >= 0)
    @variable(m, H2_GS_L[GasStorages, Hours] >= 0)
    @variable(m, H2_GS_IN[GasStorages, Hours] >= 0)
    @variable(m, H2_GS_OUT[GasStorages, Hours] >= 0)
    @variable(m, N_P2G[P2G] >= 0)
    @variable(m, N_G2P[G2P] >= 0)
    @variable(m, N_GS[GasStorages] >= 0)

    @variable(m, HEAT_STO_L[BU,HP,Hours] >= 0)
    @variable(m, HEAT_HP[BU,HP, Hours] >= 0)
    @variable(m, HEAT_INF[BU,HP, Hours] >= 0)

    prog = Progress(7, 0.01, "Building Model...         ", 40)

    @objective(m, Min,

        sum(MarginalCost[n,t] * G[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours)

        + sum(LoadIncreaseCost[n,t] * G_UP[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours)
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
            (sum(InvestmentCost[n,t] * N[(n,t)] for (n,t) in Nodes_Techs)
            + sum(InvestmentCostPower[n,sto] * N_STO_P[(n,sto)] for (n,sto) in Nodes_Storages)
            + sum(InvestmentCostEnergy[n,sto] * N_STO_E[(n,sto)] for (n,sto) in Nodes_Storages)

            + sum(FixedCost[n,t] * N[(n,t)] for (n,t) in Nodes_Techs)
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

    # @constraint(m, test, 1 >= 0)

    @constraint(m, EnergyBalance[n=DemandNodes,h=Hours],
        sum(G[(n,t),h] for (n,t) in Nodes_Dispatch)
        ## TODO : Expand to NonDispatchable / Renewables
        + sum(STO_OUT[(n,sto),h] for (n,sto) in Nodes_Storages)
        + sum(EV_DISCHARGE[ev,h] for ev in EV)
        + sum(H2_G2P[g2p,h] for g2p in G2P)
        ==
        sum(STO_IN[(n,sto),h] for (n,sto) in Nodes_Storages)
        + sum(EV_CHARGE[ev,h] for ev in EV)
        + sum(H2_P2G[p2g,h] for p2g in P2G)
        + sum(HEAT_HP[bu,hp,h] for bu in BU, hp in HP)
        + Load[n,h]
    );

    next!(prog)

    @constraint(m, MaxGenerationDisp[(n,t)=Nodes_Dispatch,h=Hours],
        G[(n,t),h] <= N[(n,t)]
        );

    @constraint(m, MaxGenerationNondisp[(n,t)=Nodes_NonDispatch,h=Hours],
        G[(n,t),h] + CU[(n,t),h] == Availability[n,t,h] * N[(n,t)]
        );

    @constraint(m, MaxInstallableBound[(n,t)=Nodes_Techs; !(MaxInstallable[n,t] |> ismissing)],
        N[(n,t)] <= MaxInstallable[n,t]
        );

    @constraint(m, MaxEnergyGenerated[(n,t)=Nodes_Techs; !(MaxEnergy[n,t] |> ismissing)],
        sum(G[(n,t),h] for h in Hours) <= corr_factor * MaxEnergy[n,t]
        );

    next!(prog)


    @constraint(m, MinRES,
        sum(G[(n,r),h] for (n,r) in Nodes_Renew, h in Hours)
        + sum(STO_OUT[(n,sto),h] - STO_IN[(n,sto),h] for (n,sto) in Nodes_Storages, h in Hours)
        + sum(H2_G2P[g2p,h] - H2_P2G[p2g,h] for p2g in P2G, g2p in G2P, h in Hours)

        >=
        (min_res/100) *(
            sum(G[(n,t),h] for (n,t)=Nodes_Techs, h in Hours)
            + sum(STO_OUT[(n,sto),h] - STO_IN[(n,sto),h] for (n,sto) in Nodes_Storages, h in Hours)
            + sum(H2_G2P[g2p,h] - H2_P2G[p2g,h] for p2g in P2G, g2p in G2P, h in Hours)
            )
        );


    @constraint(m, MaxWithdrawStorage[(n,sto) in Nodes_Storages,h=Hours],
        STO_IN[(n,sto),h] <= N_STO_P[(n,sto)]
        );

    @constraint(m, MaxGenerationStorage[(n,sto) in Nodes_Storages,h=Hours],
        STO_OUT[(n,sto),h] <= N_STO_P[(n,sto)]
        );

    @constraint(m, MaxLevelStorage[(n,sto) in Nodes_Storages,h=Hours],
        STO_L[(n,sto),h] <= N_STO_E[(n,sto)]
        );

    @constraint(m, MaxEnergyStorage[(n,sto) in Nodes_Storages; MaxEnergy[(n,sto)] >= 0],
        N_STO_E[(n,sto)] <= MaxEnergy[(n,sto)]
        );

    @constraint(m, MaxPowerStorage[(n,sto) in Nodes_Storages; MaxPower[(n,sto)] >= 0],
        N_STO_P[(n,sto)] <= MaxPower[(n,sto)]
        );

    @constraint(m, StorageBalance[(n,sto) in Nodes_Storages,h=Hours2],
        STO_L[(n,sto), Hours[h]]
        ==
        STO_L[(n,sto), Hours[h-1]]
        + STO_IN[(n,sto), Hours[h]]*sqrt(Efficiency[(n,sto)])
        - STO_OUT[(n,sto), Hours[h]]/sqrt(Efficiency[(n,sto)])
        );

    @constraint(m, StorageBalanceFirstHour[(n,sto) in Nodes_Storages],
        STO_L[(n,sto),Hours[1]]
        ==
        STO_L[(n,sto),Hours[end]]
        + STO_IN[(n,sto),Hours[1]]*sqrt(Efficiency[(n,sto)])
        - STO_OUT[(n,sto),Hours[1]]/sqrt(Efficiency[(n,sto)])
        );

    next!(prog)

    @constraint(m, MaxWithdrawEV[ev=EV,h=Hours],
        EV_CHARGE[ev,h] <= EvPower[ev][h]
        );

    @constraint(m, MaxGenerationEV[ev=EV,h=Hours],
        EV_DISCHARGE[ev,h] <= EvPower[ev][h]
        );

    @constraint(m, MaxLevelEV[ev=EV,h=Hours],
        EV_L[ev,h] <= EvCapacity[ev]
        );

    @constraint(m, EVBalance[ev=EV,h=Hours2],
        EV_L[ev, Hours[h]]
        ==
        EV_L[ev, Hours[h-1]]
        + EV_CHARGE[ev, Hours[h]]*sqrt(Efficiency[ev])
        - EV_DISCHARGE[ev, Hours[h]]/sqrt(Efficiency[ev])
        - EvDemand[ev][Hours[h]]
        + EV_PHEVFUEL[ev, Hours[h]] * EvType[ev]
        + EV_INF[ev, Hours[h]]
        );

    @constraint(m, EVBalanceFirstHour[ev=EV],
        EV_L[ev,Hours[1]]
        ==
        EV_L[ev,Hours[end]]
        + EV_CHARGE[ev,Hours[1]]*sqrt(Efficiency[ev])
        - EV_DISCHARGE[ev,Hours[1]]/sqrt(Efficiency[ev])
        - EvDemand[ev][Hours[1]]
        + EV_PHEVFUEL[ev, Hours[1]] * EvType[ev]
        + EV_INF[ev, Hours[1]]
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
        sum(Efficiency[p2g]*H2_P2G[p2g,h] for p2g in P2G)
        + sum(H2_GS_OUT[gs,h] for gs in GasStorages)
        ==
        sum((1/Efficiency[g2p])*H2_G2P[g2p,h] for g2p in G2P)
        + sum(H2_GS_IN[gs,h] for gs in GasStorages)
        + H2Demand
        );

    @constraint(m, GasstorageBalance[gs=GasStorages,h=Hours2],
        H2_GS_L[gs, Hours[h]]
        ==
        H2_GS_L[gs, Hours[h-1]]
        + H2_GS_IN[gs,Hours[h]]
        - H2_GS_OUT[gs,Hours[h]]
        );

    @constraint(m, GasstorageBalanceFirstHours[gs=GasStorages],
        H2_GS_L[gs, Hours[1]]
        ==
        H2_GS_L[gs, Hours[end]]
        + H2_GS_IN[gs,Hours[1]]
        - H2_GS_OUT[gs,Hours[1]]
        );

    next!(prog)

    @constraint(m, HeatBalance[bu=BU, hp=HP, h=Hours[2:length(Hours)]],
        HEAT_STO_L[bu,hp,h]
        ==
        StaticEfficiency[hp] * HEAT_STO_L[bu,hp,h-nthhour]
        + COP[bu,hp][h] * HEAT_HP[bu,hp,h]
        - HeatConsumption[bu,hp][h]
        + HEAT_INF[bu,hp,h]
        );


    @constraint(m, HeatBalanceFirstHour[bu=BU, hp=HP],
        HEAT_STO_L[bu,hp,Hours[1]]
        ==
        StaticEfficiency[hp] * HEAT_STO_L[bu,hp,Hours[end]]
        + COP[bu,hp][Hours[1]] * HEAT_HP[bu,hp,Hours[1]]
        - HeatConsumption[bu,hp][Hours[1]]
        + HEAT_INF[bu,hp,Hours[1]]
        );

    @constraint(m, MaxHeatLevel[bu=BU, hp=HP, h=Hours],
        HEAT_STO_L[bu,hp,h] <= MaxLevel[bu]
        );

    @constraint(m, MaxHeatPower[bu=BU, hp=HP, h=Hours],
        HEAT_HP[bu,hp,h] <= HeatMaxPower[bu,hp]
        );

    next!(prog)

    dtr.model = m

    return dtr
end

function solve_model!(dtr::DieterModel)
    @info "Start optimization"
    optimize!(dtr.model)

    return dtr
end

function generate_results!(dtr::DieterModel)
    @info "Storing results"

    vars = [
        :G => [:Nodes_Techs, :Hours],
        :G_UP => [:Nodes_Dispatch, :Hours],
        :G_DO => [:Nodes_Dispatch, :Hours],
        :CU => [:Nodes_NonDispatch, :Hours],
        :STO_IN => [:Nodes_Storages, :Hours],
        :STO_OUT => [:Nodes_Storages, :Hours],
        :STO_L => [:Nodes_Storages, :Hours],
        :N => [:Nodes_Techs],
        :N_STO_E => [:Nodes_Storages],
        :N_STO_P => [:Nodes_Storages],
        :G_INF => [:Nodes,:Hours],

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
        :HEAT_HP => [:BU, :HP, :Hours],
        :HEAT_INF => [:BU, :HP, :Hours]
    ]

    m = dtr.model

    var_dict = m.obj_dict

    dtr.results = Dict(v[1] => convert_jump_container_to_df(var_dict[v[1]], dim_names=v[2]) for v in vars)
    # dtr.results = [convert_jump_container_to_df(value.(var_dict[v[1]]), dim_names=v[2]) for v in vars]

    sum(dtr.results[:G_INF][!, :Hours]) > 0 && @warn "Problem might be infeasable"

    return nothing
end
