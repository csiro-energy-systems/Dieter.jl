# using DataFrames
# using Dieter
# import Serialization

# run_timestamp = scen_settings[:scen]
# run_timestamp = "2020-04-04-H18-Scen1_BAU-ScYr2030-Testing"
#
# results_filename = run_timestamp*"-results-Julia_Serial.dat"
# results_read_filename = joinpath(resultsdir,results_filename)
# #
# # resultsIndex = generate_results!(dtr)
# # res = dtr.results
# # # or:
# res = Serialization.deserialize(results_read_filename)

function remove_all_empty!(dict::Dict)
    for (k,v) in dict
        if isempty(v)
            delete!(dict,k)
        end
    end
end

# %% Start transforming for output:
resSplit = deepcopy(res)

delete!(resSplit,:Z)
remove_all_empty!(resSplit)

# %% Add Curtailment
Nodes_Avail_Techs = dtr.sets[:Nodes_Avail_Techs]
Hours = dtr.sets[:Hours]
Availability = dtr.parameters[:Availability]
N_TECH_Dict = Dict(eachrow(res[:N_TECH]))

df_Gen_Avail = filter(row -> row[:Nodes_Techs] in Nodes_Avail_Techs, res[:G])
# G_Dict = Dict{Any,Float64}()
# for row in eachrow(df_Gen_Avail)
#     G_Dict[row[:Nodes_Techs],row[:Hours]] = row[:Value]
# end

df_CU = @byrow! df_Gen_Avail begin
            @newcol CU::Array{Float64}
            (n,t) = :Nodes_Techs
            h = :Hours
            G = :Value
            :CU = Availability[n,t,h] * time_ratio * N_TECH_Dict[(n,t)] - G
        end

resSplit[:CU] = df_CU[!,[:Nodes_Techs, :Hours, :CU]]
DataFrames.rename!(resSplit[:CU],Dict(:CU => :Value))

# %% Define how to split Tuples in the results.

# The first Symbol is the column to split,
# the second are the new split column names,
# the third are the remaining columns to keep:
resultsSplitIndex = [
    # :Z => [],
    :G => [:Nodes_Techs, [:Nodes, :Technologies],[:Hours,:Value]],
    :CU => [:Nodes_Techs, [:Nodes, :Technologies],[:Hours,:Value]],
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
    :FLOW => [:Arcs,[:From,:To],[:Hours,:Value]],
    :H2_P2G => [:Nodes_P2G, [:Nodes, :Technologies],[:Hours,:Value]],
    :H2_G2P => [:Nodes_G2P, [:Nodes, :Technologies],[:Hours,:Value]],
    :H2_GS_L => [:Nodes_GasStorages, [:Nodes, :Technologies],[:Hours,:Value]],
    :H2_GS_IN => [:Nodes_GasStorages, [:Nodes, :Technologies],[:Hours,:Value]],
    :H2_GS_OUT => [:Nodes_GasStorages, [:Nodes, :Technologies],[:Hours,:Value]],
    :N_P2G => [:Nodes_P2G, [:Nodes, :Technologies],[:Value]],
    :N_G2P => [:Nodes_G2P, [:Nodes, :Technologies],[:Value]],
    :N_GS => [:Nodes_GasStorages, [:Nodes, :Technologies],[:Value]]
]

std_sort = [:Technologies, :Nodes]
std_sort_hours = [:Technologies, :Nodes, :Hours]
sortIndex = Dict(
    :G => std_sort_hours,
    :CU => std_sort_hours,
    :G_UP => std_sort_hours,
    :G_DO => std_sort_hours,
    :G_REZ => [:REZones],
    :G_TxZ => [:TxZones],
    :STO_IN => std_sort_hours,
    :STO_OUT => std_sort_hours,
    :STO_L => std_sort_hours,
    :N_TECH => std_sort,
    :N_STO_E => std_sort,
    :N_STO_P => std_sort,
    :N_RES_EXP => [:REZones],
    :N_SYNC => [:DemandRegion],
    :FLOW => [:From,:To],
    :H2_P2G => std_sort_hours,
    :H2_G2P => std_sort_hours,
    :H2_GS_L => std_sort_hours,
    :H2_GS_IN => std_sort_hours,
    :H2_GS_OUT => std_sort_hours,
    :N_P2G => std_sort,
    :N_G2P => std_sort,
    :N_GS => std_sort
)

# %%  Splitting Tuples, e.g.
# df_flow = split_df_tuple(res[:FLOW],:Arcs,[:From,:To])

for (sym, V) in resultsSplitIndex
    if sym in keys(resSplit)
        if length(V) > 1 && typeof(V[1]) <: Symbol
            Dieter.split_df_tuple!(resSplit[sym],V[1],V[2])
            DataFrames.select!(resSplit[sym],[V[2]...,V[3]...])
            # println(sym," ",[V[2]...,V[3]...])
        end
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
resSplit[:CU_GEN] = @where(resSplit[:CU], rel_node_tech_rez.(:Nodes,:Technologies) .== true)

resSplit[:REZ_GEN_CU] = join(resSplit[:REZ_GEN], resSplit[:CU_GEN], on=[:Nodes, :Technologies,:Hours], kind = :left)

# resSplit[:REZ_GEN_CU] = @byrow! resSplit[:REZ_GEN_CU] if ismissing(:CU); :CU = 0 end
resSplit[:REZ_GEN_CU] = @byrow! resSplit[:REZ_GEN_CU] begin
            @newcol AvailCap::Array{Float64}
            n = :Nodes
            t = :Technologies
            h = :Hours
            if ismissing(:CU)
                :CU = 0
                :AvailCap = time_ratio * N_TECH_Dict[(n,t)]
            else #if !ismissing(:CU)
                :AvailCap = Availability[n,t,h] * time_ratio * N_TECH_Dict[(n,t)]
            end
        end

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
  :REZ_GEN_CU => :Nodes,
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
