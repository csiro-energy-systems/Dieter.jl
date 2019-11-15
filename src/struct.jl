
using JuMP

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

    results::Dict{String,<:Any}

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
            :co2 => 71,   # Carbon Price
            :ev => 0,
            :heat => 0,
            :H2 => 0,
            :min_res => 0,  # Minimum Renewable Share (%)
            :cu_cost => 0   # Curtailment Cost
        )
        default_settings[:datapath] = datapath
        # default_settings[:interest] = interest


        initial_settings = merge(default_settings, settings)

        # initialise results

        results = Dict{String,Any}()
        # results = Dict{Symbol,DataFrame}()

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
"Build the data for the model, returning intermediate DataFrames containing parsed data"
function parse_data_to_model!(dtr::AbstractDieterModel)  # datapath::AbstractString

    datapath = dtr.settings[:datapath]

    files = dtr.data["files"]

    dfDict = Dict{String,DataFrame}()

    ## Base data
    # e.g files["tech"] = joinpath(datapath,"base","technologies.csv")
    dfDict["tech"] = parse_file(files["tech"])
    parse_base_technologies!(dtr, dfDict["tech"])

    # e.g. files["storage"] = joinpath(datapath,"base","storages.csv")
    dfDict["storage"] = parse_file(files["storage"])
    parse_storages!(dtr, dfDict["storage"])

    # e.g. files["load"] = joinpath(datapath,"base","load.csv")
    dfDict["load"] = parse_file(files["load"])
    parse_load!(dtr, dfDict["load"])

    # e.g. files["avail"] = joinpath(datapath,"base","availability.csv")
    dfDict["avail"] = parse_file(files["avail"])
    parse_availibility!(dtr,dfDict["avail"])

    calc_base_parameters!(dtr)

    ## EV
    # e.g. files["ev"] = joinpath(datapath,"ev","ev.csv")
    dfDict["ev"] = parse_file(files["ev"])
    parse_ev_technologies!(dtr, dfDict["ev"])

    # e.g. files["ev_demand"] = joinpath(datapath,"ev","ev_demand.csv")
    dfDict["ev_demand"] = parse_file(files["ev_demand"])
    parse_ev_demand!(dtr, dfDict["ev_demand"])

    # files["ev_power"] = joinpath(datapath,"ev","ev_power.csv")
    dfDict["ev_power"] = parse_file(files["ev_power"])
    parse_ev_power!(dtr, dfDict["ev_power"])

    calc_ev_quantity!(dtr)

    return dfDict
end

#=
    parse_heat!(dtr, joinpath(path,"heat","heat.csv"))
    parse_heat_technologies!(dtr, joinpath(path,"heat","heat_technologies.csv"))
    parse_buildings!(dtr, joinpath(path,"heat","buildings.csv"))
    parse_temperature!(dtr, joinpath(path,"heat","temperature.csv"))
    parse_heat_demand!(dtr, joinpath(path,"heat","heat_demand.csv"))
    parse_dhw_demand!(dtr, joinpath(path,"heat","dhw_demand.csv"))
    calc_hp_cop!(dtr)
    calc_heat_demand!(dtr, heat)



    parse_h2_technologies!(dtr, joinpath(path,"h2","h2_technologies.csv"))
    calc_inv_gas!(dtr)
=#

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
