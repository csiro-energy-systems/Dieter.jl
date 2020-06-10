# using DataFrames
# using Dieter
# import Serialization

# timestep = 2
timestep = dtr.settings[:timestep]
periods = round(Int,Dieter.hoursInYear*(2/timestep))
time_ratio = Dieter.hoursInYear//periods

Hours = Base.OneTo(periods)
dtr.sets[:Hours] = Hours

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
    :N_IC_EXP => [:Arcs,[:From,:To],[:Value]],
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
    :N_IC_EXP => [:From, :To],
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
node2type = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:NodeType]))

# rel_nodes_rez = Dieter.create_relation(@where(df_nodes,:NodeType .== "REZone"),:Nodes,:NodeType,:IncludeFlag)
# rel_nodes_txz = Dieter.create_relation(@where(df_nodes,:NodeType .== "TxZone"),:Nodes,:NodeType,:IncludeFlag)

REZones = dtr.sets[:REZones]
TxZones = dtr.sets[:TxZones]

# %% Technology capacity
resSplit[:CAPACITY_REZ_EXP] = resSplit[:N_RES_EXP]

ExpansionLimit = dtr.parameters[:ExpansionLimit]
resSplit[:CAPACITY_REZ_EXP] = @byrow! resSplit[:CAPACITY_REZ_EXP] begin
            @newcol ExpUpperBound::Array{Float64}
            :ExpUpperBound = ExpansionLimit[:REZones]
        end

ExpansionLimit_Tx = dtr.parameters[:ExpansionLimit_Tx]
resSplit[:CAPACITY_IC_EXP] = resSplit[:N_IC_EXP]
resSplit[:CAPACITY_IC_EXP] = @byrow! resSplit[:CAPACITY_IC_EXP] begin
            @newcol ExpUpperBound::Array{Float64}
            :ExpUpperBound = ExpansionLimit_Tx[:From,:To]
        end

df_N_STO_Include = @where(resSplit[:N_STO_P],:N_STO_P .> 1e-3)
resSplit[:CAPACITY_STO] = innerjoin(df_N_STO_Include,resSplit[:N_STO_E], on =[:Nodes, :Technologies])

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
# check: unique(resSplit[:TxZ_GEN][!,:Technologies])
# check: unique(resSplit[:TxZ_GEN][!,[:Nodes,:Technologies]])

# # Dispatch information
# Filter for absent technologies

# resSplit[:G_UP] = @where(resSplit[:G_UP], rel_node_tech_built.(:Nodes,:Technologies) .== true)
# resSplit[:G_DO] = @where(resSplit[:G_DO], rel_node_tech_built.(:Nodes,:Technologies) .== true)
# # Joins
# tmp1 = innerjoin(resSplit[:G_UP],resSplit[:G_DO], on =[:Nodes, :Technologies,:Hours])
# resSplit[:DISPATCH] = innerjoin(resSplit[:G], tmp1, on =[:Nodes, :Technologies,:Hours])

# %% Generation by renewables:

# Filter out absent technologies:
resSplit[:REZ_GEN] = @where(resSplit[:G], rel_node_tech_rez.(:Nodes,:Technologies) .== true)
resSplit[:CU_GEN] = @where(resSplit[:CU], rel_node_tech_rez.(:Nodes,:Technologies) .== true)

resSplit[:REZ_GEN_CU] = leftjoin(resSplit[:REZ_GEN], resSplit[:CU_GEN], on=[:Nodes, :Technologies,:Hours])

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
# Filter out absent technologies:
# Note: STO_L is energy (MWh), STO_IN / STO_OUT is energy in time (power)
resSplit[:STO_L] = @where(resSplit[:STO_L], rel_node_sto_built.(:Nodes,:Technologies) .== true)
resSplit[:STO_IN] = @where(resSplit[:STO_IN], rel_node_sto_built.(:Nodes,:Technologies) .== true)
resSplit[:STO_OUT] = @where(resSplit[:STO_OUT], rel_node_sto_built.(:Nodes,:Technologies) .== true)
# Joins
tmp2 = innerjoin(resSplit[:STO_IN],resSplit[:STO_OUT], on =[:Nodes, :Technologies,:Hours])
resSplit[:STORAGE] = innerjoin(resSplit[:STO_L], tmp2, on =[:Nodes, :Technologies,:Hours])

# %% Hydrogen:
if !ismissing(dtr.settings[:h2])
    # H2 Power-to-gas:
    df_p2g = @byrow! resSplit[:N_P2G] begin
                @newcol TechType::Array{String}
                :TechType = "P2G"
            end
    df_p2g = rename(df_p2g, Dict(:N_P2G => :Capacity))

    # H2 Gas-to-power
    df_g2p = @byrow! resSplit[:N_G2P] begin
                @newcol TechType::Array{String}
                :TechType = "G2P"
            end
    df_g2p = rename(df_g2p, Dict(:N_G2P => :Capacity))

    # H2 gas storage
    df_gs = @byrow! resSplit[:N_GS] begin
                @newcol TechType::Array{String}
                :TechType = "GS"
            end
    df_gs = rename(df_gs, Dict(:N_GS => :Capacity))

    df_h2 = vcat(df_p2g, df_g2p, df_gs)
    resSplit[:CAPACITY_H2] = select(df_h2, [:TechType, :Nodes, :Technologies, :Capacity])
end

# %% Augment with DemandRegion (State) information:

DemandRegion_augmentIndex = Dict(
  :TxZ_GEN => :Nodes,
  :REZ_GEN => :Nodes,
  :REZ_GEN_CU => :Nodes,
  :STORAGE => :Nodes,
  :CAPACITY_GEN => :Nodes,
  :CAPACITY_STO => :Nodes,
  :CAPACITY_REZ_EXP => :REZones,
  :CAPACITY_H2 => :Nodes
)

# Potential - annotate REZones with corresponding Transmission Zone
TxZone_augmentIndex = Dict(
  :REZ_GEN_CU => :Nodes,
  :CAPACITY_REZ_EXP => :REZones
  # :FLOW => :From
)
#  :G => :Nodes,  :G_UP => :Nodes,  :G_DO => :Nodes,  :G_REZ => :REZones,  :G_TxZ => :TxZones,
#  :STO_IN => :Nodes,  :STO_OUT => :Nodes,  :STO_L => :Nodes,  :N_TECH => :Nodes,
#  :N_STO_E => :Nodes,  :N_STO_P => :Nodes,  :N_RES_EXP => :REZones,
#  :FLOW => :From, :FLOW => [:From,:To]
#  :G_REZ => [:REZones],  :G_TxZ => [:TxZones],  :N_RES_EXP => [:REZones],

df_nodes_rez = filter(row -> row[:NodeType] == "REZone", df_nodes)
rez2TxZone = Dict(zip(df_nodes_rez[!,:Nodes],df_nodes_rez[!,:NodePromote]))

function REZ_TxZone_map!(df::DataFrame,input_colname::Symbol, output_colname::Symbol)
    df[!,output_colname] = map(x -> rez2TxZone[x], df[!,input_colname])
    return nothing
end

for (ks, vs) in TxZone_augmentIndex
    # Add the DemandRegion column
    REZ_TxZone_map!(resSplit[ks],vs,:TxZone)
    # Make the newly created (last) column the first
    L = ncol(resSplit[ks])
    resSplit[ks] = select(resSplit[ks],[L,1:(L-1)...])
end

node2DemReg = Dict(zip(df_nodes[!,:Nodes],df_nodes[!,:DemandRegion]))

function DemandRegion_map!(df::DataFrame,input_colname::Symbol, output_colname::Symbol)
    df[!,output_colname] = map(x -> node2DemReg[x], df[!,input_colname])
    return nothing
end

for (ks, vs) in DemandRegion_augmentIndex
    if ks in keys(resSplit)
        # Add the DemandRegion column
        DemandRegion_map!(resSplit[ks],vs,:DemandRegion)
        # Make the newly created (last) column the first
        L = ncol(resSplit[ks])
        resSplit[ks] = select(resSplit[ks],[L,1:(L-1)...])
    end
end

# Augment with a column showing the Demand Regions of the flow's nodes:
df_flow = resSplit[:FLOW]
DemandRegion_map!(df_flow, :From, :FromRegion)
DemandRegion_map!(df_flow, :To, :ToRegion)
# # Equivalent to:
# df_flow[!,:FromRegion] = map(x -> node2DemReg[x], df_flow[!,:From])
# df_flow[!,:ToRegion] = map(x -> node2DemReg[x], df_flow[!,:To])
L = ncol(df_flow)
resSplit[:FLOW] = select(df_flow,[L-1,L,1:(L-2)...])

# Construct inter-state / inter-region flow DataFrame:
resSplit[:INTERFLOW] = @where(resSplit[:FLOW], :FromRegion .!== :ToRegion)

DemandRegion_map!(resSplit[:N_IC_EXP], :From, :FromRegion)
DemandRegion_map!(resSplit[:N_IC_EXP], :To, :ToRegion)
L = ncol(resSplit[:N_IC_EXP])
resSplit[:N_IC_EXP] = select(resSplit[:N_IC_EXP],[L-1,L,1:(L-2)...])
