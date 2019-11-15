function parse_ev_technologies!(dtr::DieterModel, df::DataFrame)
# function parse_ev_technologies!(dtr::DieterModel, path::AbstractString)
    # df = CSV.read(path)
    dtr.sets[:ElectricVehicles] = disallowmissing(unique(df[!,:EV]))

    params = map_idcol(df)
    for (k,v) in params update_dict!(dtr.parameters, k, v) end

    return nothing
end


function parse_ev_demand!(dtr::DieterModel, df::DataFrame)
# function parse_ev_demand!(dtr::DieterModel, path::AbstractString)
    # params = CSV.read(path) |> map_dfheader_to_col
    params = map_dfheader_to_col(df)
    update_dict!(dtr.parameters, :EvDemand, params)

    return nothing
end

function parse_ev_power!(dtr::DieterModel, df::DataFrame)
# function parse_ev_power!(dtr::DieterModel, path::AbstractString)
    # params = CSV.read(path) |> map_dfheader_to_col
    params = map_dfheader_to_col(df)
    update_dict!(dtr.parameters, :EvPower, params)

    return nothing
end

function calc_ev_quantity!(dtr::DieterModel)
    p = dtr.parameters
    EV = dtr.sets[:ElectricVehicles]
    s = p[:EvShare]
    q = dtr.settings[:ev] * 1e6

    p[:AbsoluteEvDemand] = Dict(ev => s[ev] * q .* p[:EvDemand][ev] for ev in EV)
    p[:AbsoluteEvPower] = Dict(ev => s[ev] * q .* p[:EvPower][ev] for ev in EV)
    p[:AbsoluteEvCapacity] = Dict(ev => s[ev] * q * p[:EvCapacity][ev] for ev in EV)

    return nothing
end
