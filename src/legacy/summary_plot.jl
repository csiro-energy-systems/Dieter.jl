
function get_result(rdir::AbstractString, var::Union{String, Symbol}; filter_by::Dict=Dict())

    df = joinpath(rdir,String(var)*".feather") |> Feather.read

    for n in names(df)
        if typeof(df[!,n]) <: Arrow.BitPrimitive
            df[!,n] = df[!,n] |> collect
        end
    end

    for (k,v) in filter_by
        if typeof(v) <: AbstractArray
            df = filter(x-> x[k] in v, df)
        else
            df = filter(x-> x[k] == v, df)
        end
    end

    return df
end

function c_gradient(g::Array{RGB{Float64}, 1}, n::Int)
    if n > 1
        return RGB[get(ColorScheme(g), z) for z=range(1, stop=0, length=n)]
    elseif n==1
        return RGB[get(ColorScheme(g), 1.0)]
    end
end

function c_gradient(g::ColorScheme, n::Int)
    if n > 1
        return RGB[get(g, z) for z=range(1, stop=0, length=n)]
    elseif n==1
        return RGB[get(g, 1.0)]
    end
end

function create_filter_dict(sector::Symbol)
    filter_dict = Dict(:ev => 0, :heat => 0, :h2 => 0)
    delete!(filter_dict, sector)
    return filter_dict
end


function plot_generation_investments(rdir::AbstractString,
    color_dict::Dict;
    techs=["Solar PV", "Wind Onshore", "Wind Offshore", "Biomasse"],
    sector::Symbol,
    diff=false,
    kwargs...)

    capacity = get_result(rdir, :N, filter_by=create_filter_dict(sector))

    renewable = filter(row-> row[:fuel] in techs, capacity)
    sort!(renewable, [:min_res, sector])

    renewable = by(renewable, [:min_res, sector, :fuel], :Value => sum)
    # Note: The above line creates a column-name "Value_sum" in the results.


    plt = plot(legend=:topleft)
    fuels = unique(renewable[!,:fuel])
    len_sector = unique(capacity[!,sector]) |> length

    for f in fuels
        df = renewable[renewable.fuel .== f, :]
        c = c_gradient(color_dict[f], len_sector) |> permutedims
        x = df[!,:min_res]
        y =  df[!,:Value_sum]
        g = df[!,sector]
        l = vcat(f, fill("", len_sector-1)) |> permutedims
        t = "Generation Capacity"
        if sum(df[!,:Value_sum]) > 0
            plot!(x, y; group=g, c=c, label=l, width=2, title=t, grid=false, kwargs...)
        end
    end

    return plt
end


function plot_storage_investments(rdir::AbstractString,
    color_dict::Dict;
    sector::Symbol,
    kwargs...)


    fd = create_filter_dict(sector)
    storage_power = get_result(rdir, :N_STO_P, filter_by=fd)
    storage_energy = get_result(rdir, :N_STO_E, filter_by=fd)
    p2g = get_result(rdir, :N_P2G, filter_by=fd)
    g2p = get_result(rdir, :N_G2P, filter_by=fd)
    gs = get_result(rdir, :N_GS, filter_by=fd)

    for df in [storage_power, storage_energy, p2g, g2p, gs]
        if !isempty(df)
            sort!(df, [:min_res, sector])
        end
    end

    plt = plot(legend=:topleft)
    len_sector = unique(storage_power[!,sector]) |> length

    storage_power = by(storage_power, [:min_res, sector], :Value => sum)
    x = storage_power[!,:min_res]
    y = storage_power[!,:Value_sum]
    g = storage_power[!,sector]
    c = c_gradient(color_dict["Storages"], len_sector)
    plt = plot(x, y; group=g, color=permutedims(c), width=2,
        title="Storage Power", legend=false, grid=false, kwargs...)

    if !isempty(p2g)
        p2g = by(p2g, [:min_res, sector], :Value => sum)
        x = p2g[!,:min_res]
        y = p2g[!,:Value_sum]
        g = p2g[!,sector]
        c = c_gradient(color_dict["Hydrogen"], len_sector)
        plot!(plt, x, y; group=g, color=permutedims(c), width=2,
            label="Power to gas", legend=false, grid=false, linestyle=:dash)
    end

    if !isempty(g2p)
        g2p = by(g2p, [:min_res, sector], :Value => sum)
        x = g2p[!,:min_res]
        y = g2p[!,:Value_sum]
        g = g2p[!,sector]
        c = c_gradient(color_dict["Hydrogen"], len_sector)
        plot!(plt, x,y; group=g, color=permutedims(c), width=2,
            label="Gas to Power", legend=false, grid=false, linestyle=:dot)
    end

    storage_energy = by(storage_energy, [:min_res, sector], :Value => sum)
    x = storage_energy[!,:min_res]
    y = storage_energy[!,:Value_sum]
    g = storage_energy[!,sector]
    c = c_gradient(color_dict["Storages"], len_sector)
    plt2 = plot(x, y; group=g, color=permutedims(c), width=2, legend=:topleft,
        title="Storage Energy", grid=false, kwargs...)

    # if !isempty(gs)
    #     gs = by(gs, [:min_res, sector], :Value => sum)
    #     x = gs[:min_res]
    #     y = gs[:Value_sum]
    #     g = gs[sector]
    #     c = c_gradient(color_dict["Hydrogen"], len_sector)
    #     plot!(plt2, x,y; group=g, color=permutedims(c), width=2,
    #         label="h2 Storage", legend=false, grid=false, linestyle=:dash)
    # end

    return plt, plt2
end


function plot_curtailment(rdir::AbstractString,
    color_dict::Dict;
    sector::Symbol,
    kwargs...)

    fd = create_filter_dict(sector)
    cu = get_result(rdir, "summary_generation", filter_by=fd)
    sort!(cu, [:min_res, sector])

    plt = plot(legend=false)
    len_sector = unique(cu[!,sector]) |> length

    x = cu[!,:min_res]
    y = cu[!,:relative_curtailment]
    g = cu[!,sector]
    c = c_gradient(color_dict["Curtailment"], len_sector)
    plt = plot(x, y; group=g, color=permutedims(c), width=2, legend=:topleft,
        title="Curtailment [%]", grid=false, kwargs...)

    return plt
end

function plot_all(rdir::AbstractString,
    c::Dict;
    sector::Symbol,
    kwargs...)

    cap = plot_generation_investments(rdir, c, sector=sector)
    pow, ener = plot_storage_investments(rdir, c, sector=sector)
    cur = plot_curtailment(rdir, c, sector=sector)

    plt = plot(cap, cur, pow, ener; grid=(2,2), kwargs...)
    return plt
end
