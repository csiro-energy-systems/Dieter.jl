
using DataFrames
using Dieter
import Serialization
import CSV
import XLSX

prepare_df_xlsx(df) = ( collect(DataFrames.eachcol(df)), DataFrames.names(df) )

function remove_all_empty!(dict::Dict)
    for (k,v) in dict
        if isempty(v)
            delete!(dict,k)
        end
    end
end
# %% Get the results:

if Base.Sys.isapple()
      resultsdir = joinpath(ENV["HOME"],"Documents/Projects/ESM/","results_STABLE")
elseif Base.Sys.iswindows()
      resultsdir = joinpath("F:\\STABLE\\","results_STABLE")
end

run_timestamp = scen_settings[:scen]
results_filename = run_timestamp*"-results-Julia_Serial.dat"
results_read_filename = joinpath(resultsdir,results_filename)
#

#
# resultsIndex = generate_results!(dtr)
# res = dtr.results

# res = Serialization.deserialize(results_read_filename)

xlsx_output_file = joinpath(resultsdir,"STABLE_summary-$(run_timestamp).xlsx")

# %% Start transforming for output:
resSplit = deepcopy(res)

delete!(resSplit,:Z)
remove_all_empty!(resSplit)

# %% Define how to split Tuples in the results.

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
    :N_SYNC => [],
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
    :N_SYNC => [:DemandRegion],
    :FLOW => [:From,:To]
)

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

# %% Filters
df_nodes = dtr.data["dataframes"]["nodes"]
node2type = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:NodeTypes]))

# rel_nodes_rez = Dieter.create_relation(@where(df_nodes,:NodeTypes .== "REZone"),:Nodes,:NodeTypes,:IncludeFlag)
# rel_nodes_txz = Dieter.create_relation(@where(df_nodes,:NodeTypes .== "TxZone"),:Nodes,:NodeTypes,:IncludeFlag)

REZones = dtr.sets[:REZones]
TxZones = dtr.sets[:TxZones]

# %% Technology capacity
resSplit[:CAPACITY_REZ_EXP] = resSplit[:N_RES_EXP]

df_N_STO_Include = @where(resSplit[:N_STO_P],:N_STO_P .> 1e-3)
resSplit[:CAPACITY_STO] = join(df_N_STO_Include,resSplit[:N_STO_E], on =[:Nodes, :Technologies])

df_N_TECH_Include = @where(resSplit[:N_TECH],:N_TECH .> 1e-3)
resSplit[:CAPACITY_GEN] = df_N_TECH_Include

# Find tech. in renewable energy zones
df_N_TECH_rez = filter(row -> row[:Nodes] in REZones, df_N_TECH_Include)
# Find tech. in transmission zones
df_N_TECH_txz = filter(row -> row[:Nodes] in TxZones, df_N_TECH_Include)

# Functional relations
rel_node_tech_built = create_relation(df_N_TECH_Include,:Nodes,:Technologies,:N_TECH)
rel_node_sto_built = create_relation(df_N_STO_Include,:Nodes,:Technologies,:N_STO_P)

rel_node_tech_rez = create_relation(df_N_TECH_rez,:Nodes,:Technologies,:N_TECH)
rel_node_tech_txz = create_relation(df_N_TECH_txz,:Nodes,:Technologies,:N_TECH)

# %% Generation of non-renewables:
resSplit[:G] = @where(resSplit[:G], rel_node_tech_built.(:Nodes,:Technologies) .== true)

resSplit[:TxZ_GEN] = @where(resSplit[:G], rel_node_tech_txz.(:Nodes,:Technologies) .== true)
# check: unique(resSplit[:GEN_TxZ][!,:Technologies])

# # Dispatch information
# Filter for absent technologies

# resSplit[:G_UP] = @where(resSplit[:G_UP], rel_node_tech_built.(:Nodes,:Technologies) .== true)
# resSplit[:G_DO] = @where(resSplit[:G_DO], rel_node_tech_built.(:Nodes,:Technologies) .== true)
# # Joins
# tmp1 = join(resSplit[:G_UP],resSplit[:G_DO], on =[:Nodes, :Technologies,:Hours])
# resSplit[:DISPATCH] = join(resSplit[:G], tmp1, on =[:Nodes, :Technologies,:Hours])

# %% Generation of renewables:
resSplit[:REZ_GEN] = @where(resSplit[:G], rel_node_tech_rez.(:Nodes,:Technologies) .== true)

# %% Storage:
# Filter for absent technologies
# Note: STO_L is energy (MWh), STO_IN / STO_OUT is energy in time (power)
resSplit[:STO_L] = @where(resSplit[:STO_L], rel_node_sto_built.(:Nodes,:Technologies) .== true)
resSplit[:STO_IN] = @where(resSplit[:STO_IN], rel_node_sto_built.(:Nodes,:Technologies) .== true)
resSplit[:STO_OUT] = @where(resSplit[:STO_OUT], rel_node_sto_built.(:Nodes,:Technologies) .== true)
# Joins
tmp2 = join(resSplit[:STO_IN],resSplit[:STO_OUT], on =[:Nodes, :Technologies,:Hours])
resSplit[:STORAGE] = join(resSplit[:STO_L], tmp2, on =[:Nodes, :Technologies,:Hours])

# %% Augment with DemandRegion (State) information:

DR_augmentIndex = Dict(
  :TxZ_GEN => :Nodes,
  :REZ_GEN => :Nodes,
  :STORAGE => :Nodes,
  :CAPACITY_GEN => :Nodes,
  :CAPACITY_STO => :Nodes,
  :CAPACITY_REZ_EXP => :REZones
  # :FLOW => :From
)

# Potential - annotate REZones with corresponding Transmission Zone
TxZ_augmentIndex = Dict(
  :REZ_GEN => :Nodes,
  :CAPACITY_REZ_EXP => :REZones
  # :FLOW => :From
)
#  :G => :Nodes,  :G_UP => :Nodes,  :G_DO => :Nodes,  :G_REZ => :REZones,  :G_TxZ => :TxZones,
#  :STO_IN => :Nodes,  :STO_OUT => :Nodes,  :STO_L => :Nodes,  :N_TECH => :Nodes,
#  :N_STO_E => :Nodes,  :N_STO_P => :Nodes,  :N_RES_EXP => :REZones,
#  :FLOW => :From, :FLOW => [:From,:To]
#  :G_REZ => [:REZones],  :G_TxZ => [:TxZones],  :N_RES_EXP => [:REZones],

node2DemReg = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:DemandRegion]))

function DR_map!(df::DataFrame,colname::Symbol)
    demreg = map(x -> node2DemReg[x], df[!,colname])
    df[!,:DemandRegion] = demreg
    return nothing
end

for (ks, vs) in DR_augmentIndex
    # Add the DemandRegion column
    DR_map!(resSplit[ks],vs)
    # Make the new column the first
    L = ncol(resSplit[ks])
    resSplit[ks] = select(resSplit[ks],[L,1:(L-1)...])
end

# Augment with a column showing the Demand Regions of the flow's nodes:
df_flow = resSplit[:FLOW]
df_flow[!,:FromRegion] = map(x -> node2DemReg[x], df_flow[!,:From])
df_flow[!,:ToRegion] = map(x -> node2DemReg[x], df_flow[!,:To])

L = ncol(df_flow)
resSplit[:FLOW] = select(df_flow,[L-1,L,1:(L-2)...])

# Inter-state flow:
df_interflow = @where(df_flow, :FromRegion .!== :ToRegion)
DemandRegions = dtr.sets[:DemandRegions]
df_regflow = Dict{String,DataFrame}([r => DataFrame() for r in DemandRegions])

for DR in DemandRegions
    df_regflow[DR] = @linq df_interflow |>
              where(:FromRegion .== DR) |>
              by(:Hours, Level = sum(:FLOW))
end

# XLSX.openxlsx("test.xlsx", mode="w") do xf
#     for DR in DemandRegions
#         XLSX.addsheet!(xf,"INTERFLOW_"*DR)
#            XLSX.writetable!(xf["INTERFLOW_"*DR], prepare_df_xlsx(df_regflow[DR])...)
#     end
# end
# %% Write out:

# Write capacity data and create place-holders for larger datasets:
XLSX.openxlsx(xlsx_output_file, mode="w") do xf
    XLSX.rename!(xf[1],"OPDEMAND")
    XLSX.addsheet!(xf,"CAPACITY_GEN")
        XLSX.writetable!(xf["CAPACITY_GEN"], prepare_df_xlsx(resSplit[:CAPACITY_GEN])...)
    XLSX.addsheet!(xf,"CAPACITY_STO")
        XLSX.writetable!(xf["CAPACITY_STO"], prepare_df_xlsx(resSplit[:CAPACITY_STO])...)
    XLSX.addsheet!(xf,"CAPACITY_REZ_EXP")
        XLSX.writetable!(xf["CAPACITY_REZ_EXP"], prepare_df_xlsx(resSplit[:CAPACITY_REZ_EXP])...)
    XLSX.addsheet!(xf,"CAPACITY_SYNC")
        XLSX.writetable!(xf["CAPACITY_SYNC"], prepare_df_xlsx(resSplit[:N_SYNC])...)
    XLSX.addsheet!(xf,"GEN_TxZ")
        # XLSX.writetable!(xf["GEN_TxZ"], prepare_df_xlsx(resSplit[:TxZ_GEN])...)
    XLSX.addsheet!(xf,"GEN_REZ")
        # XLSX.writetable!(xf["GEN_REZ"], prepare_df_xlsx(resSplit[:REZ_GEN])...)
    XLSX.addsheet!(xf,"STORAGE")
        # XLSX.writetable!(xf["STORAGE"], prepare_df_xlsx(resSplit[:STORAGE])...)
    XLSX.addsheet!(xf,"FLOW")
        # XLSX.writetable!(xf["FLOW"], prepare_df_xlsx(resSplit[:FLOW])...)
    for DR in DemandRegions
        XLSX.addsheet!(xf,"INTERFLOW_"*DR)
           XLSX.writetable!(xf["INTERFLOW_"*DR], prepare_df_xlsx(df_regflow[DR])...)
    end
end

CSV.write(joinpath(resultsdir,"$(run_timestamp)-Gen-TxZ.csv"),resSplit[:TxZ_GEN])
CSV.write(joinpath(resultsdir,"$(run_timestamp)-Gen-REZ.csv"),resSplit[:REZ_GEN])

CSV.write(joinpath(resultsdir,"$(run_timestamp)-Storage.csv"),resSplit[:STORAGE])
CSV.write(joinpath(resultsdir,"$(run_timestamp)-Flow.csv"),resSplit[:FLOW])

# CSV.write(joinpath(resultsdir,"$(run_timestamp)-Dispatch.csv"),resSplit[:DISPATCH])
# XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
#     XLSX.writetable!(xf["DISPATCH"], prepare_df_xlsx(resSplit[:DISPATCH])...)
# end

# XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
#     XLSX.writetable!(xf["STORAGE"], prepare_df_xlsx(resSplit[:STORAGE])...)
# end

# XLSX.writetable(xlsx_output_file,
#     CAPACITY_GEN=prepare_df_xlsx(df_N_TECH_Include),
#     CAPACITY_STO=prepare_df_xlsx(resSplit[:N_STO_all])
#     )
