
function parse_base_technologies!(dtr::DieterModel, path::AbstractString)
    df = CSV.read(path)
    dtr.sets[:Technologies] = disallowmissing(unique(df[:Technologies]))
    dtr.sets[:Regions] = disallowmissing(unique(df[:Region]))

    dtr.sets[:Renewables] = disallowmissing([row[:Technologies] for row in eachrow(df)
        if row[:Renewable] == 1])

    dtr.sets[:Conventional] = disallowmissing([row[:Technologies] for row in eachrow(df)
        if row[:Renewable] == 0])

    dtr.sets[:Dispatchable] = disallowmissing([row[:Technologies]
        for row in eachrow(df) if row[:Dispatchable] == 1])

    dtr.sets[:NonDispatchable] = disallowmissing([row[:Technologies]
        for row in eachrow(df) if row[:Dispatchable] == 0])

    params = map_idcol(df, skip=[:Region, :Renewable, :Dispatchable])
    merge!(dtr.parameters, params)

    return nothing
end

function parse_storages!(dtr::DieterModel, path::AbstractString)
    df = CSV.read(path)
    dtr.sets[:Storages] = disallowmissing(unique(df[:Storages]))

    params = map_idcol(df, skip=[:Region])
    for (k,v) in params update_dict!(dtr.parameters, k, v) end

    return nothing
end

function parse_load!(dtr::DieterModel, path::AbstractString)
    df = CSV.read(path)
    dtr.parameters[:Load] = disallowmissing(df[:Load])

    return nothing
end

function parse_availibility!(dtr::DieterModel, path::AbstractString)
    params = CSV.read(path) |> map_dfheader_to_col
    update_dict!(dtr.parameters, :Availability, params)

    return nothing
end

annuity(i,lifetime) = i*((1+i)^lifetime) / (((1+i)^lifetime)-1)

function calc_inv_tech!(dtr::DieterModel)
    T = dtr.sets[:Technologies]
    oc = dtr.parameters[:OvernightCost]
    lt = dtr.parameters[:Lifetime]
    i = dtr.settings[:interest]

    dict = Dict(t => oc[t]*annuity(i, lt[t]) for t in T)
    update_dict!(dtr.parameters, :InvestmentCost, dict)
end


function calc_inv_storages!(dtr::DieterModel)
    S = dtr.sets[:Storages]
    oce = dtr.parameters[:OvernightCostEnergy]
    ocp = dtr.parameters[:OvernightCostPower]
    lt = dtr.parameters[:Lifetime]
    i = dtr.settings[:interest]

    energy = Dict(s => oce[s]*annuity(i, lt[s]) for s in S)
    power = Dict(s => ocp[s]*annuity(i, lt[s]) for s in S)

    update_dict!(dtr.parameters, :InvestmentCostEnergy, energy)
    update_dict!(dtr.parameters, :InvestmentCostPower, power)

    return nothing
end

function calc_mc!(dtr::DieterModel)
    T = dtr.sets[:Technologies]
    fc = dtr.parameters[:FuelCost]
    eff = dtr.parameters[:Efficiency]
    vc = dtr.parameters[:VariableCost]
    cc = dtr.parameters[:CarbonContent]
    co2 = dtr.settings[:co2]

    dict = Dict(t => fc[t]/eff[t] + (cc[t]*co2)/eff[t] + vc[t] for t in T)
    update_dict!(dtr.parameters, :MarginalCost, dict)
end


function calc_base_parameters!(dtr::DieterModel)
    calc_inv_tech!(dtr)
    calc_inv_storages!(dtr)
    calc_mc!(dtr)
    return nothing
end
