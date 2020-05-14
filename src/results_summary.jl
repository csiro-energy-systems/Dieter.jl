# This file: Produce summaries of various aspects of the results:

import XLSX
using OrderedCollections
import OrderedCollections: OrderedDict

# cf. write.jl
# We annotate with the label "Grouped" after summarising the results into technology groupings:
xlsx_output_file = joinpath(resultsdir,"STABLE_summary-Grouped-$(scenario_timestamp).xlsx")
# xlsx_output_file = joinpath(resultsdir,"STABLE_summary-Grouped-Output-Template-Testing.xlsx")

# Utility functions:

prepare_df_xlsx(df) = ( collect(DataFrames.eachcol(df)), DataFrames.names(df) )


# "This has the obvious caveat that the input Dict needs to be one-to-one."
# function invert_Dict(d::Dict)
#     if length(unique(values(d))) !== length(values(d))
#         @warn "The values of the input `Dict` are not unique, and the output is ambiguous."
#     end
#     return Dict([(v,k) for (k,v) in d])
# end

# We need to summarise all technology types into certain smaller categories:

summary_tech_map = OrderedDict{Symbol,String}(
    :BlackCoal    => "Black coal",
    :BrownCoal    => "Brown coal",
    :GasOther     => "Gas recip.",
    :CCGT         => "CCGT",
    :OCGT         => "OCGT",
    :Distillate   => "Diesel",
    :Biomass      => "Biomass",
    :Hydro        => "Hydro",
    :Wind         => "Wind",
    :SolarPV      => "Solar PV"
    # :SolarThermal => "Solar thermal",
)

summary_storage_map = OrderedDict{Symbol,String}(
    :PumpedHydro  => "Pumped hydro",
    :Battery      => "Large-scale battery",
    :BatteryVPP   => "VPP battery"
)
#     :Interflow => "Imports",
#     :Load => "Demand"
# )

summary_map = (union(summary_tech_map, summary_storage_map))
summary_tech = [k for (k,v) in summary_map]

summary_dict = Dict(summary_map)

# summary_reverse_map = invert_Dict(summary_tech_map)

# The map_output table has the `values` of summary_map as column headings.

fileDict["map_output"] = joinpath(datapath,"base","map_output.sql")
dfDict["map_output"] = parse_file(fileDict["map_output"]; dataname=dataname)


techGroups = Dict()
for k in summary_tech
    df = dfDict["map_output"]
    techGroups[k] = [df[row,:TechID] for row in 1:size(df,1) if df[row,k] == 1]
    # Alternative ?:
    # techGroups[k] = [row[:TechID] for row in eachrow(df) if row[k] == 1]
end

# Run checks that the map_output data corresponds to our data model:
function collateAll(dc)
    stac = []
    for k in keys(dc)
        push!(stac,dc[k]...)
    end
    return stac
end

techColl = collateAll(techGroups)
unknown_tech = setdiff(techColl, dtr.sets[:Technologies], dtr.sets[:Storages])
@assert isempty(unknown_tech)

missed_tech = setdiff(union(dtr.sets[:Technologies], dtr.sets[:Storages]),techColl)
# @assert isempty(missed_tech)

# Find Capacity for each tech. group and region

function map_group(dc, lookup)
    found_groups = []
    for k in keys(dc)
        if lookup in dc[k]
            push!(found_groups, k)
        end
    end
    if isempty(found_groups)
        @warn "Lookup value not found"
    elseif length(found_groups) > 1
        @warn "More than one group found for lookup value."
    end

    return found_groups
end

# e.g. map_group(techGroups,"Hydro_New")

function create_output_frame(tech_group_dict, df_input,
        values_col::Symbol, column_dimension::Symbol; transpose_flag=false, map_names::Bool=false, map_dict=Dict())
    df = @byrow! df_input begin
                    @newcol TechGroup::Array{String}
                    group = map_group(tech_group_dict,:Technologies)
                    if length(group) == 1
                        if map_names == true
                            :TechGroup = map_dict[group[1]]
                        else
                            :TechGroup = String(group[1])
                        end
                    end
                end

    df_sum = by(df, [column_dimension, :TechGroup], values_col => sum)

    values_rename = Symbol(String(values_col)*"_sum")

    if transpose_flag == false
        df_pivot = unstack(df_sum, :TechGroup, column_dimension, values_rename)
    else
        df_pivot = unstack(df_sum, column_dimension, :TechGroup, values_rename)
    end

    df_pivot = coalesce.(df_pivot, 0)

    return df_pivot
end

df_cap_tech_pivot = create_output_frame(techGroups, resSplit[:CAPACITY_GEN],
                            :N_TECH, :DemandRegion, map_names=true, map_dict=summary_dict)
df_cap_sto_power = create_output_frame(techGroups, resSplit[:CAPACITY_STO],
                            :N_STO_P, :DemandRegion, map_names=true, map_dict=summary_dict)
df_cap_sto_energy = create_output_frame(techGroups, resSplit[:CAPACITY_STO],
                            :N_STO_E, :DemandRegion, map_names=true, map_dict=summary_dict)

df_cap_rez = @where(resSplit[:CAPACITY_GEN], rel_node_tech_rez.(:Nodes,:Technologies) .== true)
df_cap_rez_pivot = create_output_frame(techGroups, df_cap_rez,
                            :N_TECH, :Nodes, transpose_flag=true, map_names=true, map_dict=summary_dict)

XLSX.openxlsx(xlsx_output_file, mode="w") do xf
    XLSX.rename!(xf[1],"DG_Cap_Tech")
    # XLSX.addsheet!(xf,"DG_Cap_Tech")
        XLSX.writetable!(xf["DG_Cap_Tech"], prepare_df_xlsx(df_cap_tech_pivot)...)
    XLSX.addsheet!(xf,"DS_Cap_Power")
        XLSX.writetable!(xf["DS_Cap_Power"], prepare_df_xlsx(df_cap_sto_power)...)
    XLSX.addsheet!(xf,"DS_Cap_Energy")
        XLSX.writetable!(xf["DS_Cap_Energy"], prepare_df_xlsx(df_cap_sto_energy)...)
    XLSX.addsheet!(xf,"DG_Cap_REZ")
        XLSX.writetable!(xf["DG_Cap_REZ"], prepare_df_xlsx(df_cap_rez_pivot)...)
    XLSX.addsheet!(xf,"DG_Gen")
end

# %%  Annual Generation for each tech. group and region

function scale_dataframe(df::DataFrame,exclude_col::Symbol,scalar::Number)
    fixed_col = df[!,exclude_col]
    scaled_df = mapcols(x -> scalar.*x, select(df,Not(exclude_col)))
    insertcols!(scaled_df,1, exclude_col => fixed_col)
    return scaled_df
end

df_txz_gen = resSplit[:TxZ_GEN]
df_rez_gen = resSplit[:REZ_GEN]

df_txz_gen_pivot = create_output_frame(techGroups, resSplit[:TxZ_GEN],
                            :G, :DemandRegion, map_names=true, map_dict=summary_dict)
df_rez_gen_pivot = create_output_frame(techGroups, resSplit[:REZ_GEN],
                            :G, :DemandRegion, map_names=true, map_dict=summary_dict))

df_gen_all = vcat(df_txz_gen_pivot, df_rez_gen_pivot)

df_gen_GWh = scale_dataframe(df_gen_all, :TechGroup, 1e-3)

df_sto_in_pivot = create_output_frame(techGroups, resSplit[:STORAGE],
                            :STO_OUT, :DemandRegion, map_names=true, map_dict=summary_dict)
df_sto_in_GWh = scale_dataframe(df_sto_in_pivot, :TechGroup, 1e-3)

append!(df_gen_GWh,df_sto_in_GWh)

XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
        XLSX.writetable!(xf["DG_Gen"], prepare_df_xlsx(df_gen_GWh)...)
end

# %%  Find min. and max. renewable days in each region,
#  with corresponding traces for all tech. groups, imports and demand.

# Find a one-day period that satisfies the given criteria:
# e.g. # function find_period_by_criteria(resultsDict)

df_rez = resSplit[:REZ_GEN_CU]

# unique(df_re,[:Nodes,:Technologies])

df_var = @linq df_rez |>
            where(.!occursin.(r"Hydro",:Technologies))
             # |>
            # select(:DemandRegion, :Technologies, :Hours, :AvailCap)

df_rez_days = @byrow! df_rez begin
                @newcol Day::Array{Int64}
                :Day = div(:Hours-1,48/timestep) + 1
            end

# hrs = select(df_var,:Hours)
df_var_days = @byrow! df_var begin
                @newcol Day::Array{Int64}
                :Day = div(:Hours-1,48/timestep) + 1
            end
df_var_sum = by(df_var_days, [:DemandRegion, :Day], :AvailCap => sum)
df_var_sum_pivot = unstack(df_var_sum, :Day, :DemandRegion, :AvailCap_sum)

x_days = Dict{Symbol,Dict{String,Int64}}()
df_x_days = Dict{Symbol,Dict{String,DataFrame}}()

df_x_results = Dict{Symbol,Dict{String,DataFrame}}()

make_window(n,window_size) = (n-window_size):1:(n+window_size)
# make_window(4,2) |> collect

window_size = 1
test_day = 320
test_window = make_window(test_day,window_size)

# Use inter-region flow DataFrame and filter for each region to obtain _net_ flow _from_ region:
interflow_dict = Dict{String, DataFrame}()
for DR in dtr.sets[:DemandRegions]
    interflow_dict[DR] = @linq resSplit[:INTERFLOW] |>
              where(:FromRegion .== DR) |>
              by(:Hours, Level = sum(:FLOW))
end

DemandRegions = dtr.sets[:DemandRegions]
DR_Symbols = Symbol.(dtr.sets[:DemandRegions])
for DemReg in DR_Symbols
    # DemReg = :TAS1
    # day_type = "max_day"

    df_dr = select(df_var_sum_pivot, [:Day, DemReg])
    # names(df_dr) |> display
    sort!(df_dr, DemReg)
    x_days[DemReg] = Dict{String,Int64}()
    x_days[DemReg]["min_day"] = df_dr[!,:Day][1]
    x_days[DemReg]["max_day"] = df_dr[!,:Day][end]
    # TODO: check the end days are not first and last periods of year.

    # min_day = x_days[DemReg]["min_day"]
    # max_day = x_days[DemReg]["max_day"]

    df_x_days[DemReg] = Dict{String,DataFrame}()
    df_x_results[DemReg] = Dict{String,DataFrame}()

    for day_type in ["min_day", "max_day"]
        df_x_days[DemReg][day_type] = @linq df_rez_days |>
            where(:DemandRegion .== String(DemReg)) |>
            where(in.(:Day,[make_window(x_days[DemReg][day_type],window_size)]))

        df_x_rez = select(df_x_days[DemReg][day_type],names(df_txz_gen))

        hours_list = sort(unique(df_x_rez[!,:Hours]))

        df_x_txz = @where(df_txz_gen,
                         :DemandRegion .== String(DemReg),
                         in.(:Hours,[hours_list])
                     )

        df_x = vcat(df_x_rez,df_x_txz)

        df_x_pivot = create_output_frame(techGroups, df_x, :G, :Hours, transpose_flag=true)

        df_demand = @where(dfDict["load"],
                          :DemandRegion .== String(DemReg),
                          in.(:TimeIndex,[hours_list])
                    )

        df_flow = @where(interflow_dict[String(DemReg)],
                        in.(:Hours,[hours_list])
                    )

        insertcols!(df_x_pivot, ncol(df_x_pivot)+1, :Imports => df_flow[!,:Level])
        insertcols!(df_x_pivot, ncol(df_x_pivot)+1, :Demand => df_demand[!,:Load])

        df_sto_dr = @where(resSplit[:STORAGE],
                          :DemandRegion .== String(DemReg),
                          in.(:Hours,[hours_list])
                    )

        df_sto_dr = @byrow! df_sto_dr begin
                        @newcol STO_NET::Array{Float64}
                        :STO_NET = :STO_OUT - :STO_IN
                    end

        df_s = create_output_frame(techGroups, df_sto_dr, :STO_NET, :Hours, transpose_flag=true)

        df_x_final = join(df_x_pivot, df_s, on=[:Hours])

        for tech in summary_tech
            if !(tech in names(df_x_final))
                insertcols!(df_x_final, ncol(df_x_final)+1, tech => zeros(size(df_x_final)[1]))
            end
        end

        select!(df_x_final,:Hours,summary_tech...,:Imports,:Demand)

        rename!(df_x_final, summary_map)

        xlsx_sheet_name = "DG_$(String(DemReg))_$(day_type)"
        XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
            XLSX.addsheet!(xf,xlsx_sheet_name)
                XLSX.writetable!(xf[xlsx_sheet_name], prepare_df_xlsx(df_x_final)...)
        end

        df_x_results[DemReg][day_type] = df_x_final

    end
end

# Eye-balling plots:
# x_days
# @df df_var_sum_pivot plot(:Day, cols(DR_Symbols), legend=:outertopright)
# @df df_var_sum_pivot plot(:Day, [:SA1], legend=:none)
