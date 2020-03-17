
using DataFrames
import Serialization
import CSV
import XLSX

prepare_df_xlsx(df) = ( collect(DataFrames.eachcol(df)), DataFrames.names(df) )
# %% Get the results:

run_timestamp = "2020-02-12-results"
instance_name = "2020-02-12-results.dat"

if Base.Sys.isapple()
      resultsdir = joinpath(ENV["HOME"],"Documents/Projects/ESM/","results_STABLE")
elseif Base.Sys.iswindows()
      resultsdir = joinpath("F:\\STABLE\\","results_STABLE")
end

# resultsIndex = generate_results!(dtr)
# res = dtr.results
results_read_filename = joinpath(resultsdir,instance_name)
res = Serialization.deserialize(results_read_filename)

xlsx_output_file = joinpath(resultsdir,"STABLE_summary-$(run_timestamp).xlsx")

# %% Start transforming for output:
resSplit = deepcopy(res)

delete!(resSplit,:Z)
Dieter.remove_all_empty!(resSplit)

# The first Symbol is the column to split,
# the second are the new split column names,
# the third are the remaining columns to keep:
resultsSplitIndex = [
    # :Z => [],
    :G => [:Nodes_Techs, [:Nodes, :Technologies],[:Hours,:Value]],
    :G_UP => [:Nodes_Dispatch, [:Nodes, :Technologies],[:Hours,:Value]],
    :G_DO => [:Nodes_Dispatch, [:Nodes, :Technologies],[:Hours,:Value]],
    :G_REZ => [],
    :G_TxZ => [],
    :STO_IN => [:Nodes_Storages, [:Nodes, :Technologies],[:Hours,:Value]],
    :STO_OUT => [:Nodes_Storages, [:Nodes, :Technologies],[:Hours,:Value]],
    :STO_L => [:Nodes_Storages, [:Nodes, :Technologies],[:Hours,:Value]],
    :N_TECH => [:Nodes_Techs, [:Nodes, :Technologies],[:Value]],
    :N_STO_E => [:Nodes_Storages, [:Nodes, :Technologies],[:Value]],
    :N_STO_P => [:Nodes_Storages, [:Nodes, :Technologies],[:Value]],
    :N_RES_EXP => [],
    :FLOW => [:Arcs,[:From,:To],[:Hours,:Value]]
]

std_sort = [:Technologies, :Nodes,:Hours]
sortIndex = Dict(
    :G => std_sort,
    :G_UP => std_sort,
    :G_DO => std_sort,
    :G_REZ => [:REZones],
    :G_TxZ => [:TxZones],
    :STO_IN => std_sort,
    :STO_OUT => std_sort,
    :STO_L => std_sort,
    :N_TECH => [:Technologies, :Nodes],
    :N_STO_E => [:Technologies, :Nodes],
    :N_STO_P => [:Technologies, :Nodes],
    :N_RES_EXP => [:REZones],
    :FLOW => [:From,:To]
)
# resultsSplitDict = Dict{Symbol,Vector{<:Any}}(resultsSplitIndex)

# %% Write results to a regular format using Excel .xlsx files

# Splitting Tuples, e.g.
# df_flow = split_df_tuple(res[:FLOW],:Arcs,[:From,:To])

for (sym, V) in resultsSplitIndex
    if length(V) > 1 && typeof(V[1]) <: Symbol
        Dieter.split_df_tuple!(resSplit[sym],V[1],V[2])
        DataFrames.select!(resSplit[sym],[V[2]...,V[3]...])
        # println(sym," ",[V[2]...,V[3]...])
    end
end

for (sym, df) in resSplit
    DataFrames.rename!(df,Dict(:Value => sym))
    DataFrames.sort!(df,sortIndex[sym])
end

# %% Technology capacity

N_TECH_Include = @where(resSplit[:N_TECH],:N_TECH .> 1e-3)
N_STO_Include = @where(resSplit[:N_STO_P],:N_STO_P .> 1e-3)

resSplit[:N_STO_all] = join(resSplit[:N_STO_P],resSplit[:N_STO_E], on =[:Nodes, :Technologies])

# XLSX.writetable(xlsx_output_file,
#     CAPACITY_GEN=prepare_df_xlsx(N_TECH_Include),
#     CAPACITY_STO=prepare_df_xlsx(resSplit[:N_STO_all])
#     )

# Write capacity data and create place-holders for larger datasets:
XLSX.openxlsx(xlsx_output_file, mode="w") do xf
    XLSX.rename!(xf[1],"CAPACITY_GEN")
        XLSX.writetable!(xf["CAPACITY_GEN"], prepare_df_xlsx(N_TECH_Include)...)
    XLSX.addsheet!(xf,"CAPACITY_STO")
        XLSX.writetable!(xf["CAPACITY_STO"], prepare_df_xlsx(resSplit[:N_STO_all])...)
    XLSX.addsheet!(xf,"CAPACITY_REZ_EXP")
        XLSX.writetable!(xf["CAPACITY_REZ_EXP"], prepare_df_xlsx(resSplit[:N_RES_EXP])...)
    XLSX.addsheet!(xf,"DISPATCH")
        # XLSX.writetable!(xf["DISPATCH"], prepare_df_xlsx(resSplit[:DISPATCH])...)
    XLSX.addsheet!(xf,"STORAGE")
        # XLSX.writetable!(xf["STORAGE"], prepare_df_xlsx(resSplit[:STORAGE])...)
    XLSX.addsheet!(xf,"FLOW")
        # XLSX.writetable!(xf["FLOW"], prepare_df_xlsx(resSplit[:FLOW])...)
end

CSV.write(joinpath(resultsdir,"$(run_timestamp)-Flow.csv"),resSplit[:FLOW])


# %% Dispatch:
# Filter for absent technologies
rel_node_tech_built = create_relation(N_TECH_Include,:Nodes,:Technologies,:N_TECH)
resSplit[:G] = @where(resSplit[:G], rel_node_tech_built.(:Nodes,:Technologies) .== true)
resSplit[:G_UP] = @where(resSplit[:G_UP], rel_node_tech_built.(:Nodes,:Technologies) .== true)
resSplit[:G_DO] = @where(resSplit[:G_DO], rel_node_tech_built.(:Nodes,:Technologies) .== true)
# Joins
tmp1 = join(resSplit[:G_UP],resSplit[:G_DO], on =[:Nodes, :Technologies,:Hours])
resSplit[:DISPATCH] = join(resSplit[:G], tmp1, on =[:Nodes, :Technologies,:Hours])

CSV.write(joinpath(resultsdir,"$(run_timestamp)-Dispatch.csv"),resSplit[:DISPATCH])

# XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
#     XLSX.writetable!(xf["DISPATCH"], prepare_df_xlsx(resSplit[:DISPATCH])...)
# end


# %% Storage:
# Filter for absent technologies
rel_node_sto_built = create_relation(N_STO_Include,:Nodes,:Technologies,:N_STO_P)
# Note: STO_L is energy (MWh), STO_IN / STO_OUT is energy in time (power)
resSplit[:STO_L] = @where(resSplit[:STO_L], rel_node_sto_built.(:Nodes,:Technologies) .== true)
resSplit[:STO_IN] = @where(resSplit[:STO_IN], rel_node_sto_built.(:Nodes,:Technologies) .== true)
resSplit[:STO_OUT] = @where(resSplit[:STO_OUT], rel_node_sto_built.(:Nodes,:Technologies) .== true)
# Joins
tmp2 = join(resSplit[:STO_IN],resSplit[:STO_OUT], on =[:Nodes, :Technologies,:Hours])
resSplit[:STORAGE] = join(resSplit[:STO_L], tmp2, on =[:Nodes, :Technologies,:Hours])

CSV.write(joinpath(resultsdir,"$(run_timestamp)-Storage.csv"),resSplit[:STORAGE])

# XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
#     XLSX.writetable!(xf["STORAGE"], prepare_df_xlsx(resSplit[:STORAGE])...)
# end
