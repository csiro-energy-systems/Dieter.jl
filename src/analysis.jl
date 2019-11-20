


function get_installed_capacities(dtr::DieterModel)
    df = copy(dtr.results[:N])
    df[!, :Category] = [dtr.parameters[:Fuel][x] for x in df[!,:Technologies]]

    df = aggregate(df[!,[:Category, :Value]], :Category, sum)
    rename!(df, :Value_sum => :InstalledPower)

    return df
end

function get_generation(dtr::DieterModel)

    df = copy(dtr.results[:G])
    hours = unique(df[!,:Hours])
    df[!,:Category] = [dtr.parameters[:Fuel][x] for x in df[!,:Technologies]]
    df = aggregate(df[!,[:Category,:Value]], :Category, sum)
    rename!(df, :Value_sum => :Generation)

    return df
end

function summarize_result(dtr::DieterModel; del_zeros=true)

    cap = get_installed_capacities(dtr)
    gen = get_generation(dtr)
    results = dtr.results
    df = join(cap, gen, on=:Category)

    push!(df,
        (Category = "Generation Technologies",
        InstalledPower = sum(df[!,:InstalledPower]),
        Generation = sum(df[!,:Generation])))

    push!(df,
        (Category = "Storage",
        InstalledPower = sum(results[:N_STO_P][!,:Value]),
        Generation = sum(results[:STO_OUT][!,:Value] .- results[:STO_IN][!,:Value])))

    push!(df,
        (Category = "Storage Energy",
        InstalledPower = sum(results[:N_STO_E][!,:Value]),
        Generation = maximum(results[:STO_L][!,:Value])))

    if !(dtr.settings[:h2] |> ismissing)
        push!(df,
            (Category = "Gas2Power",
            InstalledPower = sum(results[:N_G2P][!,:Value]),
            Generation = sum(results[:H2_G2P][!,:Value])))

        push!(df,
            (Category = "Power2Gas",
            InstalledPower = -sum(results[:N_P2G][!,:Value]),
            Generation = -sum(results[:H2_P2G][!,:Value])))
    end

    if !(dtr.settings[:ev] |> ismissing)
        push!(df,
            (Category = "EV Demand",
            InstalledPower = -maximum(results[:EV_CHARGE][!,:Value]),
            Generation = -sum(results[:EV_CHARGE][!,:Value])))

        push!(df,
            (Category = "EV Discharge",
            InstalledPower = maximum(results[:EV_DISCHARGE][!,:Value]),
            Generation = sum(results[:EV_DISCHARGE][!,:Value])))
    end

    if !(dtr.settings[:heat] |> ismissing)
        push!(df,
            (Category = "Heat Demand",
            InstalledPower = -maximum(results[:HEAT_HP][!,:Value]),
            Generation = -sum(results[:HEAT_HP][!,:Value])))
    end

    hours = unique(dtr.results[:G][!,:Hours])
    load = dtr.parameters[:Load][hours]

    push!(df,
        (Category = "Load",
        InstalledPower = -maximum(load),
        Generation = -sum(load)))

    push!(df,
        (Category = "Curtailment",
        InstalledPower = maximum(results[:CU][!,:Value]),
        Generation = sum(results[:CU][!,:Value])))


    df[!,:InstalledPower] = round.(df[!,:InstalledPower] ./ 1000, digits=1)
    df[!,:Generation] = round.(df[!,:Generation] ./ 1e6, digits=1)

    if del_zeros
        df = df[df[!,:InstalledPower] .!= 0, :]
    end

    return df
end
