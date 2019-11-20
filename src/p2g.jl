function parse_h2_technologies!(dtr::DieterModel, df::DataFrame)
# function parse_h2_technologies!(dtr::DieterModel, path::AbstractString)
    # df = df::DataFrame.read(path)
    dtr.sets[:H2Technologies] = disallowmissing(unique(df[!,:H2Technologies]))

    dtr.sets[:P2G] = disallowmissing([row[:H2Technologies] for row in eachrow(df)
        if row[:H2Type] == "P2G"])

    dtr.sets[:G2P] = disallowmissing([row[:H2Technologies] for row in eachrow(df)
        if row[:H2Type] == "G2P"])

    dtr.sets[:GasStorages] = disallowmissing([row[:H2Technologies] for row in eachrow(df)
        if row[:H2Type] == "GasStorage"])

    params = map_idcol(df)
    for (k,v) in params update_dict!(dtr.parameters, k, v) end

    return nothing
end

function calc_inv_gas!(dtr::DieterModel)
    T = dtr.sets[:H2Technologies]
    oc = dtr.parameters[:OvernightCost]
    lt = dtr.parameters[:Lifetime]
    i = dtr.settings[:interest]

    dict = Dict(t => oc[t]*annuity(i, lt[t]) for t in T)
    update_dict!(dtr.parameters, :InvestmentCost, dict)
    # for (k,v) in dict update_dict!(dtr.parameters[:InvestmentCost], k, v) end

    return nothing
end
