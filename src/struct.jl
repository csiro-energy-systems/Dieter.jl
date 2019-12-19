
# using JuMP

"""
The `def` macro is used to build other macros that can insert the same block of
Julia code into different parts of a program.
This is macro is used to generate a standard set of fields inside a model type
hierarchy. (Original code from InfrastructureModels.jl)
"""
macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end


abstract type AbstractDieterModel end

"""
A macro for adding the basic fields to an AbstractDieterModel type definition.
"""
@def dm_fields begin
    model::JuMP.AbstractModel

    data::Dict{String,<:Any}

    sets::Dict{Symbol,Array{<:Any,1}}
    parameters::Dict{Symbol,<:Any}

    settings::Dict{Symbol,<:Any}

    results::Dict{Symbol,<:Any}

    # ref::Dict{Symbol,<:Any} # reference data
    # var::Dict{Symbol,<:Any} # JuMP variables
    # con::Dict{Symbol,<:Any} # JuMP constraint references
end

mutable struct DieterModel <: AbstractDieterModel @dm_fields end

#=
 Note: We avoid defining a custom inner constructor method for subtypes of AbstractDieterModel
 while defining the basic fields using a custom macro @dm_fields.
 If we want to enforce constraints or invariants on the fields data, we may need an inner constructor.
 This gives the option to define outer constructors as needed for AbstractDieterModel subtype.
=#

const current_path = @__DIR__

"""
Default constructor for a Dieter model, initialised with the Model DataType and
a user-specified `data` Dict associated with the model in the constructor.
"""
function InitialiseDieterModel(ModelType::Type, data::Dict{String,T} where T;
            datapath::String="",
            jump_model::JuMP.AbstractModel=JuMP.Model(),
            settings=Dict{Symbol,Any}()
         )

        @assert ModelType <: AbstractDieterModel
        # initialise data structures for the model

        sets = Dict{Symbol, Array{Any,1}}()
        parameters = Dict{Symbol, Any}()

        # initialise settings

        default_settings = Dict{Symbol, Union{String,Number}}(
            :scen => "Default",
            :datapath => "",
            :interest => 0.04,
            :co2 => 0,     # Carbon Price
            :min_res => 0,  # Minimum Renewable Share (%)
            :cu_cost => 0,   # Curtailment Cost
            :ev => 0,
            :heat => 0,
            :h2 => 0
        )
        default_settings[:datapath] = datapath
        # default_settings[:interest] = interest

        initial_settings = merge(default_settings, settings)

        # initialise results
        # results = Dict{String,Any}()
        results = Dict{Symbol,DataFrame}()

        dtr = ModelType(
            jump_model,
            data,
            sets,
            parameters,
            initial_settings,
            results
        )

    return dtr
end

# %%
"""
Initialise a basic file structure for accessing the model data.
The default directory structure is:
<datapath>/
├── base
│   ├── availability.csv
│   ├── load.csv
│   ├── storages.csv
│   └── technologies.csv
├── ev
│   ├── ev.csv
│   ├── ev__demand.csv
│   └── ev__power.csv
├── h2
│   └── h2__technologies.csv
└── heat
    ├── buildings.csv
    ├── dhw__demand.csv
    ├── heat.csv
    ├── heat__demand.csv
    ├── heat__technologies.csv
    └── temperature.csv
"""
function initialise_data_dir_structure(datapath;remove=false)
    if !(isdir(datapath))
        @info "The given input directory \"$datapath\" does not exist. It will be created."
        # return false
    end

    datadirs = ["base","ev","h2","heat"]

    for dir in datadirs
        thispath=joinpath(datapath,dir)
        if isdir(thispath)
            @info "Directory \"$dir\" in \"$datapath\"/ already exists, skipping..."
            continue
        else
            @info "Creating directory \"$dir\" in \"$datapath\"/ "
            mkpath(thispath)
        end
    end

    return true
end

"""
Initialise the strings giving the file-path references to the model data.
The argument `sfx` is the file suffix of each model data file.
Currently, `sfx` can be either "csv" or "sql"
"""
function initialise_data_file_dict!(dtr::AbstractDieterModel,sfx::String)

    datapath = dtr.settings[:datapath]

    dtr.data["files"] = Dict{String,String}()

    fileDict = dtr.data["files"]

    # %% Base
    fileDict["nodes"] = joinpath(datapath,"base","nodes."*sfx)
    fileDict["arcs"] = joinpath(datapath,"base","arcs."*sfx)
    fileDict["tech"] = joinpath(datapath,"base","technologies."*sfx)
    fileDict["storage"] = joinpath(datapath,"base","storages."*sfx)
    fileDict["load"] = joinpath(datapath,"base","load."*sfx)
    fileDict["avail"] = joinpath(datapath,"base","availability."*sfx)

    fileDict["map_node_tech"] = joinpath(datapath,"base","map_node_tech."*sfx)
    fileDict["map_node_storages"] = joinpath(datapath,"base","map_node_storages."*sfx)

    # %% EV
    if !(dtr.settings[:ev] |> ismissing)
        fileDict["ev"] = joinpath(datapath,"ev","ev."*sfx)
        fileDict["ev_demand"] = joinpath(datapath,"ev","ev_demand."*sfx)
        fileDict["ev_power"] = joinpath(datapath,"ev","ev_power."*sfx)
    end

    # %% Heat
    if !(dtr.settings[:heat] |> ismissing)
        fileDict["heat"] = joinpath(datapath,"heat","heat."*sfx)
        fileDict["heat_technologies"] = joinpath(datapath,"heat","heat_technologies."*sfx)
        fileDict["buildings"] = joinpath(datapath,"heat","buildings."*sfx)
        fileDict["temperature"] = joinpath(datapath,"heat","temperature."*sfx)
        fileDict["heat_demand"] = joinpath(datapath,"heat","heat_demand."*sfx)
        fileDict["dhw_demand"] = joinpath(datapath,"heat","dhw_demand."*sfx)
    end

    # %% Hydrogen (h2)
    if !(dtr.settings[:h2] |> ismissing)
        fileDict["h2"] = joinpath(datapath,"h2","h2_technologies."*sfx)
    end

    return dtr
end
