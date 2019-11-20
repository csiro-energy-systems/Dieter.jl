
function parse_heat!(dtr::DieterModel, df::DataFrame)
# function parse_heat!(dtr::DieterModel, path::AbstractString)
    # df = CSV.read(path)
    dtr.sets[:BuildingType] = disallowmissing(unique(df[!,:BuildingType]))
    dtr.sets[:HeatingType] = disallowmissing(unique(df[!,:HeatingType]))

    params = map_idcol(df, id_cols=2)
    merge!(dtr.parameters, params)

    return nothing
end

function parse_heat_technologies!(dtr::DieterModel, df::DataFrame)
# function parse_heat_technologies!(dtr::DieterModel, path::AbstractString)
    # df = CSV.read(path)
    params = map_idcol(df)
    merge!(dtr.parameters, params)

    return nothing
end

function parse_buildings!(dtr::DieterModel, df::DataFrame)
# function parse_buildings!(dtr::DieterModel, path::AbstractString)
    # df = CSV.read(path)

    params = map_idcol(df)
    for (k,v) in params update_dict!(dtr.parameters, k, v) end

    return nothing
end

function parse_temperature!(dtr::DieterModel, df::DataFrame)
# function parse_temperature!(dtr::DieterModel, path::AbstractString)
    # params = CSV.read(path) |> map_dfheader_to_col
    params = map_dfheader_to_col(df)
    update_dict!(dtr.parameters, :Temperature, params)

    return nothing
end

function parse_heat_demand!(dtr::DieterModel, df::DataFrame)
# function parse_heat_demand!(dtr::DieterModel, path::AbstractString)
    # params = CSV.read(path) |> map_dfheader_to_col
    params = map_dfheader_to_col(df)
    update_dict!(dtr.parameters, :HeatDemand, params)

    return nothing
end

function parse_dhw_demand!(dtr::DieterModel, df::DataFrame)
# function parse_dhw_demand!(dtr::DieterModel, path::AbstractString)
    # params = CSV.read(path) |> map_dfheader_to_col
    params = map_dfheader_to_col(df)
    update_dict!(dtr.parameters, :HotwaterDemand, params)

    return nothing
end

function calc_hp_cop!(dtr::DieterModel)

    cop = Dict{Tuple{String, String}, Array{Float64,1}}()

    for bu in dtr.sets[:BuildingType], hp in dtr.sets[:HeatingType]

        if dtr.parameters[:Outside][hp] == 1
            temp_diff = dtr.parameters[:TemperatureSink][bu] .- dtr.parameters[:Temperature][hp]
        else
            temp_diff = fill(dtr.parameters[:TemperatureSink][bu] - dtr.parameters[:TemperatureSource][bu], 8760)
        end

        delta = (dtr.parameters[:TemperatureSink][bu] + 273.15) ./ temp_diff
        cop[bu,hp] = dtr.parameters[:DynamicEfficiency][hp] .* delta
    end

    dtr.parameters[:COP] = cop

    return nothing
end


function calc_heat_demand!(dtr::DieterModel)
# function calc_heat_demand!(dtr::DieterModel, heat::Number)
    heat = dtr.settings[:heat]
    heat_demand = Dict{Tuple{String, String}, Array{Float64,1}}()
    total_share = heat/100

    for bu in dtr.sets[:BuildingType], hp in dtr.sets[:HeatingType]
        hs = dtr.parameters[:HeatShare][bu,hp]
        af = dtr.parameters[:AreaFloor][bu]
        d = dtr.parameters[:HeatDemand][bu] + dtr.parameters[:HotwaterDemand][bu]

        heat_demand[bu,hp] = total_share * (hs * (af .* d))
    end

    dtr.parameters[:HeatConsumption] = heat_demand

    return nothing
end
