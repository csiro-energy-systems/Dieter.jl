
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

    sets::Dict{Symbol,Array{String,1}}
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

        sets = Dict{Symbol, Array{String,1}}()
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
    fileDict["tech"] = joinpath(datapath,"base","technologies."*sfx)
    fileDict["storage"] = joinpath(datapath,"base","storages."*sfx)
    fileDict["load"] = joinpath(datapath,"base","load."*sfx)
    fileDict["avail"] = joinpath(datapath,"base","availability."*sfx)

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

function check_data_files_exist(fileDict::Dict{String,String})
    for (k,v) in fileDict
        if !isfile(fileDict[k])
            @error "File \"$(fileDict[k])\" not found."
            return false
        end
    end

    return true
end

"Parse a given data file and return the desired data in a `DataFrame`."
function parse_file(file::String; dataname::String="")
    filetype=split(lowercase(file), '.')[end]
    if filetype == "csv"
        df = CSV.read(file)
    elseif filetype == "sql" ## TODO
        # We expect dataname to be the path to a SQLite database file
        if !(isfile(dataname))
            @error "The given dataname \"$dataname\" is not a valid file."
        end
        try
            db = SQLite.DB(dataname)
        catch
            @error "Opening of SQLite database file \"$dataname\" unsuccessful."
        end
        db = SQLite.DB(dataname)
        queryString = read(file, String)
        df = SQLite.Query(db,queryString) |> DataFrame
    else
        error("Unrecognised filetype (as parsed from the file extension).")
    end

    return df
end

# %%
## TODO ? : Split the parsing to the DataFrame Dict from the model transformations.
"""
Build the data for the model, returning intermediate `DataFrame`s containing parsed data.
The function will not parse data where the corresponding model setting equals `missing`.
"""
function parse_data_to_model!(dtr::AbstractDieterModel; dataname::AbstractString="")

    datapath = dtr.settings[:datapath]

    fileDict = dtr.data["files"]

    dfDict = Dict{String,DataFrame}()

    ## Base data
    # e.g. fileDict["tech"] = joinpath(datapath,"base","technologies.csv")
    dfDict["tech"] = parse_file(fileDict["tech"]; dataname=dataname)
    parse_base_technologies!(dtr, dfDict["tech"])

    # e.g. fileDict["storage"] = joinpath(datapath,"base","storages.csv")
    dfDict["storage"] = parse_file(fileDict["storage"]; dataname=dataname)
    parse_storages!(dtr, dfDict["storage"])

    # e.g. fileDict["load"] = joinpath(datapath,"base","load.csv")
    dfDict["load"] = parse_file(fileDict["load"]; dataname=dataname)
    parse_load!(dtr, dfDict["load"])

    # e.g. fileDict["avail"] = joinpath(datapath,"base","availability.csv")
    dfDict["avail"] = parse_file(fileDict["avail"]; dataname=dataname)
    parse_availibility!(dtr,dfDict["avail"])

    calc_base_parameters!(dtr)

    ## EV
    if !(dtr.settings[:ev] |> ismissing)
        # e.g. fileDict["ev"] = joinpath(datapath,"ev","ev.csv")
        dfDict["ev"] = parse_file(fileDict["ev"]; dataname=dataname)
        parse_ev_technologies!(dtr, dfDict["ev"])

        # e.g. fileDict["ev_demand"] = joinpath(datapath,"ev","ev_demand.csv")
        dfDict["ev_demand"] = parse_file(fileDict["ev_demand"]; dataname=dataname)
        parse_ev_demand!(dtr, dfDict["ev_demand"])

        # e.g. fileDict["ev_power"] = joinpath(datapath,"ev","ev_power.csv")
        dfDict["ev_power"] = parse_file(fileDict["ev_power"]; dataname=dataname)
        parse_ev_power!(dtr, dfDict["ev_power"])

        calc_ev_quantity!(dtr)
    else
        @info "No Electric Vehicle aspects present in model."
        dtr.sets[:ElectricVehicles] = Array([])
        for s in [:AbsoluteEvDemand,:AbsoluteEvPower,:AbsoluteEvCapacity,:EvFuel,:EvType]
            dtr.parameters[s] = Dict()
        end
    end

    ## Heat
    if !(dtr.settings[:heat] |> ismissing)
        # e.g. fileDict["heat"] = joinpath(datapath,"heat","heat.csv")
        dfDict["heat"] = parse_file(fileDict["heat"]; dataname=dataname)
        parse_heat!(dtr, dfDict["heat"])

        # e.g. fileDict["heat_technologies"] = joinpath(datapath,"heat","heat_technologies.csv")
        dfDict["heat_technologies"] = parse_file(fileDict["heat_technologies"]; dataname=dataname)
        parse_heat_technologies!(dtr, dfDict["heat_technologies"])

        # e.g. fileDict["buildings"] = joinpath(datapath,"heat","buildings.csv")
        dfDict["buildings"] = parse_file(fileDict["buildings"]; dataname=dataname)
        parse_buildings!(dtr, dfDict["buildings"])

        # e.g. fileDict["temperature"] = joinpath(datapath,"heat","temperature.csv")
        dfDict["temperature"] = parse_file(fileDict["temperature"]; dataname=dataname)
        parse_temperature!(dtr, dfDict["temperature"])

        # e.g. fileDict["heat_demand"] = joinpath(datapath,"heat","heat_demand.csv")
        dfDict["heat_demand"] = parse_file(fileDict["heat_demand"]; dataname=dataname)
        parse_heat_demand!(dtr, dfDict["heat_demand"])

        # e.g. fileDict["dhw_demand"] = joinpath(datapath,"heat","dhw_demand.csv")
        dfDict["dhw_demand"] = parse_file(fileDict["dhw_demand"]; dataname=dataname)
        parse_dhw_demand!(dtr, dfDict["dhw_demand"])

        calc_hp_cop!(dtr)
        calc_heat_demand!(dtr)
    else
        @info "No Heat aspects present in model."
        dtr.sets[:BuildingType] = Array([])
        dtr.sets[:HeatingType]  = Array([])
        for s in [:HeatConsumption, :HeatShare, :MaxLevel, :HeatMaxPower, :StaticEfficiency, :COP]
            dtr.parameters[s] = Dict()
        end
    end

    ## Hydrogen (h2)
    if !(dtr.settings[:h2] |> ismissing)
        # e.g. fileDict["h2"] = joinpath(datapath,"h2","h2_technologies.csv")
        dfDict["h2"] = parse_file(fileDict["h2"]; dataname=dataname)
        parse_h2_technologies!(dtr, dfDict["h2"])

        calc_inv_gas!(dtr)
    else
        @info "No Power-to-Gas/Hydrogen aspects present in model."
        dtr.sets[:H2Technologies] = Array([])
        dtr.sets[:P2G] = Array([])
        dtr.sets[:G2P] = Array([])
        dtr.sets[:GasStorages] = Array([])
        # for s in []
        #     dtr.parameters[s] = Dict()
        # end
    end

    return dfDict
end

# "Build the data for the model"
# function parse_data_to_model(dtr::AbstractDieterModel, datapath::AbstractString)

#        parse_base_technologies!(dtr, df_base)

        # parse_storages!(dtr, joinpath(path,"base","storages.csv"))
        # parse_load!(dtr, joinpath(path,"base","load.csv"))
        # parse_availibility!(dtr,joinpath(path,"base","availability.csv"))
        # calc_base_parameters!(dtr)
        #
        #
        # parse_ev_technologies!(dtr, joinpath(path,"ev","ev.csv"))
        # parse_ev_demand!(dtr, joinpath(path,"ev","ev_demand.csv"))
        # parse_ev_power!(dtr, joinpath(path,"ev","ev_power.csv"))
        # calc_quantity!(dtr, ev)
        #
        #
        #
        # parse_heat!(dtr, joinpath(path,"heat","heat.csv"))
        # parse_heat_technologies!(dtr, joinpath(path,"heat","heat_technologies.csv"))
        # parse_buildings!(dtr, joinpath(path,"heat","buildings.csv"))
        # parse_temperature!(dtr, joinpath(path,"heat","temperature.csv"))
        # parse_heat_demand!(dtr, joinpath(path,"heat","heat_demand.csv"))
        # parse_dhw_demand!(dtr, joinpath(path,"heat","dhw_demand.csv"))
        # calc_hp_cop!(dtr)
        # calc_heat_demand!(dtr, heat)
        #
        #
        #
        # parse_h2_technologies!(dtr, joinpath(path,"h2","h2_technologies.csv"))
        # calc_inv_gas!(dtr)
# end


# mutable struct DieterModel
    # sets::Dict{Symbol, Array{String,1}}
    # parameters::Dict{Symbol, Any}
    # settings::Dict{Symbol, Any}
    # results::Dict{Symbol, DataFrame}
#=
    function DieterModel(path::AbstractString;
        settings::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        ev::T1=0,
        heat::T2=0,
        h2::T3=0,
        res_share::T4=0,
        cu_cost::T5=0,
        scen::T6=0,
        ) where {T1<:Number, T2<:Number, T3<:Number ,T4<:Number, T5<:Number, T6<:Union{String,Symbol, Int}}


        self = new()
        self.settings = Dict{Symbol, Any}(
            :datapath => path,
            :interest => 0.04,
            :co2 => 71,
            :ev => ev,
            :heat => heat,
            :h2 => h2,
            :min_res => res_share,
            :cu_cost => cu_cost,
            :scen => scen)

        merge!(self.settings, settings)

        self.sets = Dict{Symbol, Array{String,1}}()
        self.parameters = Dict{Symbol, Any}()
        self.results = Dict{Symbol, Any}()

        return self
    end

    DieterModel() = new()
end
=#
