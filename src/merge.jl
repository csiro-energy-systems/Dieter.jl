
const FOLDER_LIST =
    ["CU",
    "EV_CHARGE",
    "EV_DISCHARGE",
    "EV_INF",
    "EV_L",
    "EV_PHEVFUEL",
    "G",
    "G_INF",
    "H2_G2P",
    "H2_GS_IN",
    "H2_GS_L",
    "H2_GS_OUT",
    "H2_P2G",
    "HEAT_HP",
    "HEAT_INF",
    "HEAT_STO_L",
    "N",
    "N_G2P",
    "N_GS",
    "N_P2G",
    "N_STO_E",
    "N_STO_P",
    "STO_IN",
    "STO_L",
    "STO_OUT"]

"""
Remove subfolders of `rdir`, the results directory.
Note that the default is for individual variable data files to be saved to the directory
".../<save_name>/single_results/<name>"
when the function save_results(dtr::DieterModel,"<save_name>") has been called.
"""
function rm_folders(rdir::AbstractString; singledir::AbstractString="single_results")
    for folder in FOLDER_LIST
        path = joinpath(rdir,singledir,folder)
        isdir(path) && rm(path, recursive=true)
    end
end


function merge_results(rdir::AbstractString; singledir::AbstractString="single_results")

    prog = Progress(length(FOLDER_LIST), 0.01, "Merging results...   ", 40)

    for folder in FOLDER_LIST
        f = joinpath(rdir, singledir, folder)
        if isdir(f)
            dfs = [Feather.read(joinpath(f,file))
                for file in readdir(f) if rsplit(file,".";limit=2)[end] == "feather"]
            merged_df = vcat(dfs...)
            Feather.write(joinpath(rdir,folder*".feather"), merged_df)
        end
        next!(prog)
    end

    return nothing
end

function join_dfs(arr::Vector{DataFrame}; on::Array{Symbol, 1})

    merge_df = join(arr[1], arr[2], on=on)
    for df in arr[3:end]
        merge_df = join(merge_df, df, on=on, kind = :outer)
    end

    return merge_df
end

function need_to_merge(rdir::AbstractString; singledir::AbstractString="single_results")

    for f in FOLDER_LIST
        p = joinpath(rdir, singledir, f)
        pp = readdir(p)
        m = maximum([mtime(joinpath(p,x)) for x in pp])
        merged_file = joinpath(rdir, f*".feather")

        if isfile(merged_file)
            m2 = mtime(merged_file)
        else
            println("No merged file yet.")
            return true
        end

        if m2 < m
            println("Newer version of $f found.")
            return true
        end
    end

    return false
end

const COLS = [:nthhour, :scen, :min_res, :ev, :heat, :h2]

function create_generation_summary(rdir::AbstractString)

    g = Feather.read(joinpath(rdir,"G.feather"))
    cu = Feather.read(joinpath(rdir,"CU.feather"))
    sto_in = Feather.read(joinpath(rdir,"STO_IN.feather"))
    sto_out = Feather.read(joinpath(rdir,"STO_OUT.feather"))

    ev_charge = Feather.read(joinpath(rdir,"EV_CHARGE.feather"))
    ev_discharge = Feather.read(joinpath(rdir,"EV_DISCHARGE.feather"))

    heat_demand = Feather.read(joinpath(rdir,"HEAT_HP.feather"))

    p2g = Feather.read(joinpath(rdir,"H2_P2G.feather"))
    g2p = Feather.read(joinpath(rdir,"H2_G2P.feather"))

    gen_summary = DataFrame[]
    push!(gen_summary, by(g, COLS, Generation_total = :Value => sum))

    gen_by_fuel = by(g, vcat(COLS,:fuel), Generation_by_fuel = :Value => sum)
    push!(gen_summary, unstack(gen_by_fuel,  :fuel, :Generation_by_fuel))

    gen_by_res = by(g, vcat(COLS,:renewable), Generation_by_res = :Value => sum)
    gen_by_res= unstack(gen_by_res,  :renewable, :Generation_by_res)
    rename!(gen_by_res, Dict(Symbol(:false) => :conventional, Symbol(:true) => :renewable))
    push!(gen_summary, gen_by_res)

    push!(gen_summary, by(cu, COLS, Curtailment = :Value => sum))

    push!(gen_summary, by(ev_charge, COLS, Charged_by_ev = :Value => sum))
    push!(gen_summary, by(ev_discharge, COLS, Discharged_by_ev = :Value => sum))

    push!(gen_summary, by(heat_demand, COLS, Demand_from_hp = :Value => sum))

    push!(gen_summary, by(sto_in, COLS, Storage_in = :Value => sum))
    push!(gen_summary, by(sto_out, COLS, Storage_out = :Value => sum))

    push!(gen_summary, by(p2g, COLS, Power2Gas = :Value => sum))
    push!(gen_summary, by(g2p, COLS, Gas2power = :Value => sum))

    summary_generation = join_dfs(gen_summary, on=COLS)

    summary_generation[!,:res_share] =
        summary_generation[!,:renewable].*100 ./ summary_generation[!,:Generation_total]

    summary_generation[!,:relative_curtailment] =
        summary_generation[!,:Curtailment].*100 ./ summary_generation[!,:Generation_total]

    Feather.write(joinpath(rdir,"summary_generation.feather"), summary_generation)
    CSV.write(joinpath(rdir,"summary_generation.csv"), summary_generation)

end


function create_investment_summary(rdir::AbstractString)

    inv = Feather.read(joinpath(rdir,"N.feather"))
    inv_power = Feather.read(joinpath(rdir,"N_STO_P.feather"))
    inv_storage = Feather.read(joinpath(rdir,"N_STO_E.feather"))

    n_p2g = Feather.read(joinpath(rdir,"N_P2G.feather"))
    n_g2p = Feather.read(joinpath(rdir,"N_G2P.feather"))
    n_gs = Feather.read(joinpath(rdir,"N_GS.feather"))

    inv_summary = DataFrame[]

    inv = by(inv, vcat(COLS,:fuel), Investment_by_fuel = :Value => sum)
    push!(inv_summary, unstack(inv,  :fuel, :Investment_by_fuel))
    push!(inv_summary, by(inv_power, COLS, Investment_power = :Value => sum))
    push!(inv_summary, by(inv_storage, COLS, Investment_storage = :Value => sum))


    push!(inv_summary, by(n_p2g, COLS, Investment_p2g = :Value => sum))
    push!(inv_summary, by(n_g2p, COLS, Investment_g2p = :Value => sum))
    push!(inv_summary, by(n_gs, COLS, Investment_gs = :Value => sum))

    summary_investment = join_dfs(inv_summary, on=COLS)

    Feather.write(joinpath(rdir,"summary_investment.feather"), summary_investment)
    CSV.write(joinpath(rdir,"summary_investment.csv"), summary_investment)

end


function post_process_results(rdir::AbstractString)

    if need_to_merge(rdir)
        merge_results(rdir)
        create_generation_summary(rdir)
        create_investment_summary(rdir)
    else
        println("Did not find any new files.")
    end

    return nothing
end
