function parse_h2_technologies!(dtr::DieterModel, df::DataFrame)
# function parse_h2_technologies!(dtr::DieterModel, path::AbstractString)
    # df = df::DataFrame.read(path)
    dtr.sets[:H2Technologies] = disallowmissing(unique(df[!,:H2Technologies]))

    dtr.sets[:P2G] = disallowmissing(unique([row[:H2Technologies] for row in eachrow(df)
        if row[:H2Type] == "P2G"]))

    dtr.sets[:G2P] = disallowmissing(unique([row[:H2Technologies] for row in eachrow(df)
        if row[:H2Type] == "G2P"]))

    dtr.sets[:GasStorages] = disallowmissing(unique([row[:H2Technologies] for row in eachrow(df)
        if row[:H2Type] == "GasStorage"]))

    params = map_idcol(df, [:Region, :H2Technologies], skip_cols=Symbol[])
    for (k,v) in params update_dict!(dtr.parameters, k, v) end

    rel_node_h2tech = create_relation(df,:Region,:H2Technologies,:Efficiency)
    dtr.data["relations"]["rel_node_h2tech"] = rel_node_h2tech

    Nodes = dtr.sets[:Nodes]
    H2Technologies = dtr.sets[:H2Technologies]
    dtr.sets[:Nodes_H2Tech] = Dieter.tuple2_filter(rel_node_h2tech, Nodes, dtr.sets[:H2Technologies])
    dtr.sets[:Nodes_P2G] = Dieter.tuple2_filter(rel_node_h2tech, Nodes, dtr.sets[:P2G])
    dtr.sets[:Nodes_G2P] = Dieter.tuple2_filter(rel_node_h2tech, Nodes, dtr.sets[:G2P])
    dtr.sets[:Nodes_GasStorages] = Dieter.tuple2_filter(rel_node_h2tech, Nodes, dtr.sets[:GasStorages])

    return nothing
end

function calc_inv_gas!(dtr::DieterModel)
    Nodes_H2Tech = dtr.sets[:Nodes_H2Tech]

    T = dtr.sets[:H2Technologies]
    oc = dtr.parameters[:H2OvernightCost]
    lt = dtr.parameters[:Lifetime]
    i = dtr.settings[:interest]

    dict = Dict((n,t) => oc[n,t]*annuity(i, lt[n,t]) for (n,t) in Nodes_H2Tech if t in T)
    update_dict!(dtr.parameters, :InvestmentCost, dict)

    fc = dtr.parameters[:FuelCost]      # FuelCost in currency/MWh-thermal-input
    eff = dtr.parameters[:Efficiency]   # Efficiency is unitless (MWh-output/MWh-thermal-input)
    vc = dtr.parameters[:VariableCost]

    marginalcost = Dict((n,t) => fc[n,t]/eff[n,t] + vc[n,t] for (n,t) in Nodes_H2Tech)

    update_dict!(dtr.parameters, :MarginalCost, marginalcost)

    return nothing
end

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
