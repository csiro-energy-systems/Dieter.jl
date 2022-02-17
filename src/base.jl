
"""
Build the data for the model, returning intermediate `DataFrame`s containing parsed data.
The function will not parse data where the corresponding model setting equals `missing`.
"""
function parse_data_to_model!(dtr::AbstractDieterModel; dataname::AbstractString="")

    datapath = dtr.settings[:datapath]

    fileDict = dtr.data["files"]

    dtr.data["dataframes"] = Dict{String,DataFrame}()

    dfDict = dtr.data["dataframes"]

    ## Base data

    dfDict["nodes"] = parse_file(fileDict["nodes"]; dataname=dataname)
    parse_nodes!(dtr)

    # e.g. fileDict["tech"] = joinpath(datapath,"base","technologies.csv")
    dfDict["tech"] = parse_file(fileDict["tech"]; dataname=dataname)
    parse_base_technologies!(dtr)

    # e.g. fileDict["storage"] = joinpath(datapath,"base","storages.csv")
    dfDict["storage"] = parse_file(fileDict["storage"]; dataname=dataname)
    parse_storages!(dtr)

    # e.g. fileDict["load"] = joinpath(datapath,"base","load.csv")
    dfDict["load"] = parse_file(fileDict["load"]; dataname=dataname)
    parse_load!(dtr, dfDict["load"])

    # e.g. fileDict["avail"] = joinpath(datapath,"base","availability.csv")
    dfDict["avail"] = parse_file(fileDict["avail"]; dataname=dataname)
    parse_availibility!(dtr)

    # Relations (i.e set-to-set correspondences)
    dfDict["map_node_demand"] = parse_file(fileDict["map_node_demand"]; dataname=dataname)
    dfDict["map_node_tech"] = parse_file(fileDict["map_node_tech"]; dataname=dataname)
    dfDict["map_node_storages"] = parse_file(fileDict["map_node_storages"]; dataname=dataname)
    dfDict["arcs"] = parse_file(fileDict["arcs"]; dataname=dataname)
    initialise_set_relation_data!(dtr)
    parse_set_relations!(dtr)
    STABLE.parse_arcs!(dtr)

    # Calculated base parameters
    calc_base_parameters!(dtr)

    return dfDict

end

function parse_extensions!(dtr::AbstractDieterModel; dataname::AbstractString="")

    datapath = dtr.settings[:datapath]
    fileDict = dtr.data["files"]
    dfDict = dtr.data["dataframes"]

    ## EV
    if !(dtr.settings[:ev] |> ismissing)
        # e.g. fileDict["ev"] = joinpath(datapath,"ev","ev.csv")
        dfDict["ev"] = parse_file(fileDict["ev"]; dataname=dataname)
        parse_ev_technologies!(dtr, dfDict["ev"])

        # e.g. fileDict["ev_demand"] = joinpath(datapath,"ev","ev_demand.csv")
        dfDict["ev_demand"] = parse_file(fileDict["ev_demand"]; dataname=dataname)
        parse_ev_demand!(dtr, dfDict["ev_demand"])

        # e.g. fileDict["ev_power"] = joinpath(datapath,"ev","ev_power.csv")
        dfDict["ev_power"] = parse_file(fileDict["ev_power"]; dataname=dataname)
        parse_ev_power!(dtr, dfDict["ev_power"])

        calc_ev_quantity!(dtr)
    else
        @info "No Electric Vehicle aspects present in model."
        dtr.sets[:ElectricVehicles] = Array([])
        for s in [:AbsoluteEvDemand,:AbsoluteEvPower,:AbsoluteEvCapacity,:EvFuel,:EvType]
            dtr.parameters[s] = Dict()
        end
    end

    ## Heat
    if !(dtr.settings[:heat] |> ismissing)
        # e.g. fileDict["heat"] = joinpath(datapath,"heat","heat.csv")
        dfDict["heat"] = parse_file(fileDict["heat"]; dataname=dataname)
        parse_heat!(dtr, dfDict["heat"])

        # e.g. fileDict["heat_technologies"] = joinpath(datapath,"heat","heat_technologies.csv")
        dfDict["heat_technologies"] = parse_file(fileDict["heat_technologies"]; dataname=dataname)
        parse_heat_technologies!(dtr, dfDict["heat_technologies"])

        # e.g. fileDict["buildings"] = joinpath(datapath,"heat","buildings.csv")
        dfDict["buildings"] = parse_file(fileDict["buildings"]; dataname=dataname)
        parse_buildings!(dtr, dfDict["buildings"])

        # e.g. fileDict["temperature"] = joinpath(datapath,"heat","temperature.csv")
        dfDict["temperature"] = parse_file(fileDict["temperature"]; dataname=dataname)
        parse_temperature!(dtr, dfDict["temperature"])

        # e.g. fileDict["heat_demand"] = joinpath(datapath,"heat","heat_demand.csv")
        dfDict["heat_demand"] = parse_file(fileDict["heat_demand"]; dataname=dataname)
        parse_heat_demand!(dtr, dfDict["heat_demand"])

        # e.g. fileDict["dhw_demand"] = joinpath(datapath,"heat","dhw_demand.csv")
        dfDict["dhw_demand"] = parse_file(fileDict["dhw_demand"]; dataname=dataname)
        parse_dhw_demand!(dtr, dfDict["dhw_demand"])

        calc_hp_cop!(dtr)
        calc_heat_demand!(dtr)
    else
        @info "No Heat aspects present in model."
        dtr.sets[:BuildingType] = Array([])
        dtr.sets[:HeatingType]  = Array([])
        for s in [:HeatConsumption, :HeatShare, :MaxLevel, :HeatMaxPower, :StaticEfficiency, :CoP]
            dtr.parameters[s] = Dict()
        end
    end

    ## Hydrogen (h2)
    if !(dtr.settings[:h2] |> ismissing)
        # e.g. fileDict["h2"] = joinpath(datapath,"h2","h2_technologies.csv")
        dfDict["h2_technologies"] = parse_file(fileDict["h2_technologies"]; dataname=dataname)

        # Overwrite values from "h2_technologies" table by Scenario:
        ScenCostPower = dtr.parameters[:ScenCostPower]

        tech_types = dtr.settings[:tech_types]

        H2ElectrolyserTypes = tech_types[:H2ElectrolyserTypes]
        H2RecipEngTypes = tech_types[:H2RecipEngTypes]

        ## Map OvernightCost data into table
        @info "H2 Technology cost overwritten with scenario values."
        dfDict["h2_technologies"] = @eachrow dfDict["h2_technologies"] begin
            if :H2Technologies in keys(ScenCostPower)
                :H2OvernightCost = ScenCostPower[:H2Technologies]
                h2_electrolyser = :H2Technologies
                region = :Region
                # @info "H2 Technology ($(h2_electrolyser)) cost was overwritten in $(region) with scenario value $(ScenCostPower[h2_electrolyser])."
            # else
            #     @warn "H2 Electrolyser ($(h2_electrolyser)) cost was NOT overwritten by scenario value"
            end
        end

        parse_h2_technologies!(dtr, dfDict["h2_technologies"])

        # # H2 Reciprocating Engine Fuel prices
        InstanceYear = dtr.settings[:inst_year]
        InstanceYear_Sym = Symbol("FYE$(InstanceYear)")
        node2DemReg = dtr.parameters[:node_demreg_dict]

        fuel_prices_sql_db_path = dtr.settings[:datapaths]["Fuel_Price_SQL_db"]
        H2_Cost_Tech_Map = dtr.settings[:h2_cost_tech_map]
        
        if dtr.settings[:H2_G2P_Allowed_flag] == true
            fileDict["fuel_prices_h2"] = joinpath(datapath,"base","fuel_prices_h2.sql")
            dfDict["fuel_prices_h2"] = parse_file(fileDict["fuel_prices_h2"]; dataname=fuel_prices_sql_db_path)
            
            df_fuel_h2 = @subset(dfDict["fuel_prices_h2"],:TechID .== H2_Cost_Tech_Map["N_RecipH2"])
            # Rename the InstanceYear column as `:FuelCost`:
            select!(df_fuel_h2,:Region, (InstanceYear_Sym => :FuelCost))
            h2_fuel_dict = Dict(eachrow(df_fuel_h2))

            for ((n,t),v) in dtr.parameters[:FuelCost] 
                    if t == H2_Cost_Tech_Map["N_RecipH2"]
                        dtr.parameters[:FuelCost][(n,t)] = h2_fuel_dict[node2DemReg[n]]
                    end
            end
        end

        calc_inv_gas!(dtr)
    else
        @info "No Power-to-Gas/Hydrogen aspects present in model."
        dtr.sets[:H2Technologies] = Array([])
        dtr.sets[:P2G] = Array([])
        dtr.sets[:G2P] = Array([])
        dtr.sets[:GasStorages] = Array([])
        dtr.sets[:Nodes_H2Tech] = Array([])
        dtr.sets[:Nodes_P2G] = Array([])
        dtr.sets[:Nodes_G2P] = Array([])
        dtr.sets[:Nodes_GasStorages] = Array([])

        dtr.sets[:Nodes_P2G_PEM] = Array([])
        dtr.sets[:Nodes_P2G_AE]  = Array([])

        dtr.parameters[:H2Conversion] = Dict{Any,Any}()
        # for s in []
        #     dtr.parameters[s] = Dict()
        # end
    end

    return nothing
end

# Expected input e.g.  df = dtr.data["dataframes"]["nodes"]
# Columns: | Nodes │ NodeType │ NodePromote │ DemandRegion | IncludeLevel |
# Expected that NodeType has (at least) entries including ["TxZone", "REZone"]
function parse_nodes!(dtr::DieterModel)

    df = dtr.data["dataframes"]["nodes"]

    dtr.sets[:Nodes] = disallowmissing(unique(df[!, :Nodes]))
    dtr.sets[:NodeType] = disallowmissing(unique(df[!, :NodeType]))
    dtr.sets[:DemandRegions] = collect(skipmissing(unique(df[!, :DemandRegion])))
    # dtr.sets[:Regions] = disallowmissing(unique(df[!, :Region]))

    # Obtain the Transmission Zones
    dtr.sets[:TxZones] = (@subset(df,:NodeType .== "TxZone"))[!,:Nodes]
    # Obtain the Renewable Energy Zones
    dtr.sets[:REZones] = (@subset(df,:NodeType .== "REZone"))[!,:Nodes]

    # Equate Demand Zones to Transmission Zones
    dtr.sets[:DemandZones] = dtr.sets[:TxZones]

    # Provide a dict (lookup map) between nodes and the associated demand region
    dtr.parameters[:node_demreg_dict] = Dict(zip(df[!,:Nodes],df[!,:DemandRegion]))

    # Provide a dict (lookup map) between nodes and the associated "IncludeLevel", a relative measure of cost increase:
    dtr.parameters[:include_level_dict] = Dict(zip(df[!,:Nodes],df[!,:IncludeLevel]))

    return nothing
end

"""
Expected Dieter model parameters read into `dtr.parameters`:
CapacityFactor, FuelCost, FixedCost, VariableCost, OvernightCost,
CurtailmentCost, LoadIncreaseCost, LoadDecreaseCost, Lifetime,
MaxCapacity, MaxEnergy,
TechTypeCategory, Renewable, Dispatchable,
FuelSource, Efficiency, CarbonContent, CO2_price
"""
function parse_base_technologies!(dtr::DieterModel)

    df = dtr.data["dataframes"]["tech"]

    dtr.sets[:Technologies] = disallowmissing(unique(df[!, :Technologies]))

    dtr.sets[:Renewables] = disallowmissing(unique([row[:Technologies] for row in eachrow(df)
        if row[:Renewable] == 1]))

    dtr.sets[:Conventional] = disallowmissing(unique([row[:Technologies] for row in eachrow(df)
        if row[:Renewable] == 0]))

    dtr.sets[:Dispatchable] = disallowmissing(unique([row[:Technologies]
        for row in eachrow(df) if row[:Dispatchable] == 1]))

    dtr.sets[:NonDispatchable] = disallowmissing(unique([row[:Technologies]
        for row in eachrow(df) if row[:Dispatchable] == 0]))

    # params = map_idcol(df, skip_cols=[:Region, :Renewable, :Dispatchable])
    params = map_idcol(df, [:Region, :Technologies], skip_cols=Symbol[])
    merge!(dtr.parameters, params)

    # If a column of parameters is entirely missing, convert container to a Union type.
    # Otherwise, we get a parameter Dict of type Dict{Tuple{R,T},Missing} :
    DefaultType = Float64 # default type to contain general parameters
    R = typeof(df[1,:Region])
    T = typeof(df[1,:Technologies])

    for (k,v) in dtr.parameters
        if typeof(v) <: Dict{Tuple{R,T},Missing}
            dtr.parameters[k] = convert(Dict{Tuple{R,T},Union{DefaultType,Missing}},v)
        end
    end

    return nothing
end

function parse_storages!(dtr::DieterModel)

    df = dtr.data["dataframes"]["storage"]

    dtr.sets[:Storages] = disallowmissing(unique(df[!,:Storages]))

    params = map_idcol(df, [:Region, :Storages], skip_cols=Symbol[])
    for (k,v) in params
        try
            update_dict!(dtr.parameters, k, v)
        catch e
            print(e)
            error(
                "Update failed: dtr.parameters[$k] is type
                    $(typeof(dtr.parameters[k])) \n while new parameter values are of type
                    $(typeof(v))."
            )
        end
    end

    return nothing
end

function parse_load!(dtr::DieterModel, df::DataFrame)

    params = map_idcol(df, [:Nodes, :TimeIndex], skip_cols=Symbol[])
    # params = map_idcol(df, [:DemandRegion, :TimeIndex], skip_cols=Symbol[])
    for (k,v) in params update_dict!(dtr.parameters, k, v) end

    return nothing
end

function parse_availibility!(dtr::DieterModel; check_normalisation::Bool=true)

    df = dtr.data["dataframes"]["avail"]

    ## Checks
    if !isempty(setdiff([:RenewRegionID, :TechTypeID, :TimeIndex, :Availability], propertynames(df)))
        @error "Required columns were not found in the passed dataframe."
    end

    if any(df[:,:Availability] .< 0)
        error("Negative availability factor values found in column :Availability from dtr.data[\"dataframes\"][\"avail\"].")
    elseif check_normalisation && any(df[:,:Availability] .> 1)
        @info "Checking for whether availability factor values are no greater than 1."
        error("Availability factor values greater than 1 found in column :Availability from dtr.data[\"dataframes\"][\"avail\"].")
    elseif !(check_normalisation)
        @info "Skipping check for whether availability factor values are no greater than 1."
    end

    params = map_idcol(df, [:RenewRegionID, :TechTypeID, :TimeIndex], skip_cols=Symbol[])
    for (k,v) in params update_dict!(dtr.parameters, k, v) end

    # Construct the set of node/technology pairs that have availability traces.
    # This is used to validate the data before constraint construction.
    availpairs = unique(df[!,[:RenewRegionID,:TechTypeID]])
    dtr.sets[:Nodes_Avail_Techs] = [Tuple(x) for x in eachrow(availpairs)]

    if !isempty( setdiff(dtr.sets[:Nodes_Avail_Techs],dtr.sets[:Nodes_Techs]) )
          println(setdiff(dtr.sets[:Nodes_Avail_Techs],dtr.sets[:Nodes_Techs]))
          error("Available trace data includes a technology an excluded region.")
    end

    return nothing
end


"Create the mathematical set relations for doing set-correspondences in models."
function initialise_set_relation_data!(dtr)
# See util.jl/create_relation(...) for the function definition
# See struct.jl/parse_data_to_model!(...) for construction of the Dict dtr.data["maps"]

# map_<...> is a a reference to a DataFrame containing side-by-side columns of set correspondences
# rel_<...> is a filter function created from the DataFrame map_<...>
# rel_<...> applied is a Julia array object with Tuples containing the mathematical relation for later iteration.

    # Input for this scope (the inputs need to have been parsed previously using initialise_data_file_dict! & parse_data_to_model! )
    dfDict = dtr.data["dataframes"]

    # Output stored in model's data field:
    dtr.data["relations"] = Dict{String,Any}()

    # Relation between a node and associated technologies
    rel_node_tech = create_relation(dfDict["map_node_tech"],:Nodes,:Technologies,:IncludeFlag)
    # rel_node_tech = filter(x->f_node_tech(x[1],x[2]), [(n,t) for n in N for t in T] )

    # Relation between a node and associated storage technologies
    rel_node_storages = create_relation(dfDict["map_node_storages"],:Nodes,:Technologies,:IncludeFlag)
    # rel_node_storages = filter(x->f_node_storages(x[1],x[2]),[(n,s) for n in N for s in S])

    # Relation between a node and associated type
    rel_node_type = create_relation(dfDict["nodes"],:Nodes,:NodeType,:IncludeLevel)

    # Relation between a node and associated type above (its promotion)
    rel_node_promote = create_relation(dfDict["nodes"],:Nodes,:NodePromote,:IncludeLevel)

    # Relation between a node and demand
    rel_node_demand = create_relation(dfDict["map_node_demand"],:Nodes,:DemandZone,:IncludeFlag)
    # rel_node_demand = create_relation(dfDict["map_node_demand"],:Nodes,:DemandRegion,:IncludeFlag)

    # Relation between a node and incident/connected nodes (inter-connectors)
    rel_node_incidence = create_relation(dfDict["arcs"],:FromNode,:ToNode,:TransferCapacity)

    # Store the results:
    dtr.data["relations"]["rel_node_tech"] = rel_node_tech
    dtr.data["relations"]["rel_node_storages"] = rel_node_storages
    dtr.data["relations"]["rel_node_type"] = rel_node_type
    dtr.data["relations"]["rel_node_promote"] = rel_node_promote
    dtr.data["relations"]["rel_node_demand"] = rel_node_demand
    dtr.data["relations"]["rel_node_incidence"] = rel_node_incidence

    return nothing

end

function parse_set_relations!(dtr)

    Nodes = dtr.sets[:Nodes]
    NodeType = dtr.sets[:NodeType]
    DemandRegions = dtr.sets[:DemandRegions]
    DemandZones = dtr.sets[:DemandZones]

    Technologies = dtr.sets[:Technologies]
    Storages = dtr.sets[:Storages]

    Renewables = dtr.sets[:Renewables]
    Conventional = dtr.sets[:Conventional]
    Dispatchable = dtr.sets[:Dispatchable]
    NonDispatchable = dtr.sets[:NonDispatchable]

    # Relation between a node and associated technologies
    rel_node_tech = dtr.data["relations"]["rel_node_tech"]
    dtr.sets[:Nodes_Techs] = tuple2_filter(rel_node_tech, Nodes, Technologies)
    dtr.sets[:Nodes_Renew] = tuple2_filter(rel_node_tech, Nodes, Renewables)
    dtr.sets[:Nodes_Conven] = tuple2_filter(rel_node_tech, Nodes, Conventional)
    dtr.sets[:Nodes_Dispatch] = tuple2_filter(rel_node_tech, Nodes, Dispatchable)
    dtr.sets[:Nodes_NonDispatch] = tuple2_filter(rel_node_tech, Nodes, NonDispatchable)


    # Relation between a node and associated storage technologies
    rel_node_storages = dtr.data["relations"]["rel_node_storages"]
    dtr.sets[:Nodes_Storages] = tuple2_filter(rel_node_storages, Nodes, Storages)

    # Relation between a node and associated type
    rel_node_type = dtr.data["relations"]["rel_node_type"]
    dtr.sets[:Nodes_Types] = tuple2_filter(rel_node_type, Nodes, NodeType)

    # Obtain the Transmission Zones
    # dtr.sets[:TxZones] = [x[1] for x in filter(x -> (x[2]=="TxZone"),dtr.sets[:Nodes_Types])]
    # Obtain the Renewable Energy Zones
    # dtr.sets[:REZones] = [x[1] for x in filter(x -> (x[2]=="REZone"),dtr.sets[:Nodes_Types])]

    # Relation between a node and associated node above it (its "promotion") if that is also a Node
    # Typically, creates a relation of REZone Nodes to the connected TxZone.
    rel_node_promote = dtr.data["relations"]["rel_node_promote"]
    dtr.sets[:Nodes_Promotes] = tuple2_filter(rel_node_promote, Nodes, Nodes)

    # Relation between a node and associated demand
    rel_node_demand = dtr.data["relations"]["rel_node_demand"]
    dtr.sets[:Nodes_Demand] = tuple2_filter(rel_node_demand, Nodes, DemandZones)
    # dtr.sets[:Nodes_Demand] = tuple2_filter(rel_node_demand, Nodes, DemandRegions)

    return nothing

end

function parse_arcs!(dtr::DieterModel)

    df = dtr.data["dataframes"]["arcs"]

    Nodes = dtr.sets[:Nodes]

    # Relation between a node and incident/connected nodes (inter-connectors)
    rel_node_incidence = dtr.data["relations"]["rel_node_incidence"]
    dtr.sets[:Arcs] = tuple2_filter(rel_node_incidence, Nodes, Nodes)

    dtr.sets[:Arcs_From] = disallowmissing(unique(df[!,:FromNode]))

    dtr.sets[:Arcs_REZones] = [(rez,txz) for (rez,txz) in dtr.sets[:Nodes_Promotes]
                                if rez in dtr.sets[:REZones]]

    params = map_idcol(df, [:FromNode, :ToNode], skip_cols=Symbol[])
    for (k,v) in params update_dict!(dtr.parameters, k, v) end

    return nothing
end

function calc_inv_trans_exp!(dtr::DieterModel)

    tec = dtr.parameters[:TransExpansionCost]
    lt_Tx = dtr.settings[:lifetime_Tx]
    i = dtr.settings[:interest]

    dict = Dict((from,to) => tec[from,to]*annuity(i, lt_Tx) for (from,to) in keys(tec))
    update_dict!(dtr.parameters, :InvestmentCostTransExp, dict)

    return nothing
end

annuity(i,lifetime) = i*((1+i)^lifetime) / (((1+i)^lifetime)-1)
# Equivalent alternative: annuity(i,lifetime) = i +  i/( ((1+i)^lifetime)-1 )

"Calculate the amortised investment cost of technologies from capital costs."
function calc_inv_tech!(dtr::DieterModel)
    Nodes_Techs = dtr.sets[:Nodes_Techs]
    # Nodes = dtr.sets[:Nodes]
    # Technologies = dtr.sets[:Technologies]

    oc = dtr.parameters[:OvernightCostPower]
    lt = dtr.parameters[:Lifetime]
    i = dtr.settings[:interest]

    dict = Dict((n,t) => oc[n,t]*annuity(i, lt[n,t]) for (n,t) in Nodes_Techs)
    update_dict!(dtr.parameters, :InvestmentCost, dict)

    return nothing
end

"Calculate the amortised investment cost of technologies from capital costs."
function calc_inv_storages!(dtr::DieterModel)
    Nodes_Storages = dtr.sets[:Nodes_Storages]
    # Nodes = dtr.sets[:Nodes]
    # Storages = dtr.sets[:Storages]

    oce = dtr.parameters[:OvernightCostEnergy]
    ocp = dtr.parameters[:OvernightCostPower]
    lt = dtr.parameters[:Lifetime]
    i = dtr.settings[:interest]

    cost_energy = Dict((n,s) => oce[n,s]*annuity(i, lt[n,s]) for (n,s) in Nodes_Storages)
    cost_power  = Dict((n,s) => ocp[n,s]*annuity(i, lt[n,s]) for (n,s) in Nodes_Storages)

    update_dict!(dtr.parameters, :InvestmentCostEnergy, cost_energy)
    update_dict!(dtr.parameters, :InvestmentCostPower, cost_power)

    return nothing
end

"Calculated marginal costs for plants including variable costs"
function calc_mc!(dtr::DieterModel)
    Nodes_Techs = dtr.sets[:Nodes_Techs]
    Nodes_Storages = dtr.sets[:Nodes_Storages]
    # Nodes = dtr.sets[:Nodes]
    # Technologies = dtr.sets[:Technologies]
    # Storages = dtr.sets[:Storages]

    fc = dtr.parameters[:FuelCost]      # FuelCost in currency/MWh-thermal-input
    eff = dtr.parameters[:Efficiency]   # Efficiency is unitless (MWh-output/MWh-thermal-input)
    cc = dtr.parameters[:CarbonContent] # CarbonContent in units of t-CO2/MWh-thermal-input
    co2 = dtr.settings[:co2]            # Price on carbon in currency/t-CO2
    vc = dtr.parameters[:VariableCost]  # VariableCost in currency/MWh-output

    marginalcost = Dict((n,t) => fc[n,t]/eff[n,t] + (cc[n,t]*co2)/eff[n,t] + vc[n,t] for (n,t) in Nodes_Techs)
    update_dict!(dtr.parameters, :MarginalCost, marginalcost)

    marginalcost_sto = Dict((n,s) => fc[n,s]/eff[n,s] + vc[n,s] for (n,s) in Nodes_Storages)
    update_dict!(dtr.parameters, :MarginalCost, marginalcost_sto)

    return nothing
end


function calc_base_parameters!(dtr::DieterModel)
    calc_inv_tech!(dtr)
    calc_inv_storages!(dtr)
    calc_mc!(dtr)

    return nothing
end
