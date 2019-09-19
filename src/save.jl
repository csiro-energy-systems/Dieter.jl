

function create_folder(p::AbstractString)
    if isdir(p)
        return joinpath(pwd(),p)
    else
        mkdir(p)
        return joinpath(pwd(),p)
    end
end





function add_info!(dtr::DieterModel)
    r = dtr.results

    for (k,df) in r
        for s in [:min_res, :ev, :heat, :h2, :nthhour, :scen]
            df[s] = dtr.settings[s]
        end
        df[:Value] = df[:Value] ./ 1000
    end

    for k  in [:N, :G]
        r[k][:fuel] = [dtr.parameters[:Fuel][x] for x in r[k][:Technologies]]
        r[k][:renewable] = [dtr.parameters[:Renewable][x] == 1 for x in r[k][:Technologies]]
        r[k][:dispatchable] = [dtr.parameters[:Dispatchable][x] == 1 for x in r[k][:Technologies]]
    end


    return dtr
end


function save_results(dtr::DieterModel, rdir::AbstractString)

    add_info!(dtr)

    min_res = dtr.settings[:min_res]
    ev = dtr.settings[:ev]
    heat = dtr.settings[:heat]
    h2 = dtr.settings[:h2]
    nthhour = dtr.settings[:nthhour]
    scen = dtr.settings[:scen]

    fn = "-scen_$(scen)-res_$(min_res)%-ev_$(ev)-heat_$(heat)%-h2_$(h2)twh-hour_$(nthhour).csv"
    fn2 = "-scen_$(scen)-res_$(min_res)%-ev_$(ev)-heat_$(heat)%-h2_$(h2)twh-hour_$(nthhour).feather"

    create_folder(rdir)
    srdir = create_folder(joinpath(rdir,"single_results"))

    for (k,v) in dtr.results
        if !isempty(v)
            path = create_folder(joinpath(srdir,string(k)))

            # fp = joinpath(path, string(k)*fn)
            # CSV.write(fp, v)

            fp2 = joinpath(path, string(k)*fn2)
            Feather.write(fp2, v)
        end
    end

    # CSV.write(joinpath(rdir, "load.csv"), DataFrame(Load=dtr.parameters[:Load]))
    fpath = joinpath(rdir, "load.feather")
    isfile(fpath) || Feather.write(fpath, DataFrame(Load=dtr.parameters[:Load]))

    nothing
end
