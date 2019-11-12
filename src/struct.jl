
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
A macro for adding the basic fields to an AbstractDieterModel type definition
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

mutable struct BaseDieterModel <: AbstractDieterModel @dm_fields end

#=
 Note: We avoid defining a custom inner constructor method for subtypes of AbstractDieterModel
 while defining the basic fields using a custom macro @dm_fields.
 If we want to enforce constraints or invariants on the fields data, we may need an inner constructor.
 This gives the option to define outer constructors as needed for AbstractDieterModel subtype.
=#

const current_path = @__DIR__

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

        default_settings = Dict{Symbol, Any}(
            :datapath => "",
            :interest => 0.04,
            :co2 => 71,
            :ev => 0,
            :heat => 0,
            :H2 => 0,
            :min_res => 0,
            :cu_cost => 0,
            :scen => "Default"
        )
        default_settings[:datapath] = datapath
        # default_settings[:interest] = interest


        initial_settings = merge(default_settings, settings)

        # initialise results

        results = Dict{String,Any}()
        # results = Dict{Symbol,DataFrame}()

        dtrm = ModelType(
            jump_model,
            data,
            sets,
            parameters,
            initial_settings,
            results
        )

    return dtrm
end

data_instance = Dict{String,Any}()
dm = InitialiseDieterModel(BaseDieterModel, data_instance)


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

#=
function parse_data_to_model(dm:AbstractDieterModel)
        parse_base_technologies!(dm, joinpath(path,"base","technologies.csv"))
        parse_storages!(dm, joinpath(path,"base","storages.csv"))
        parse_load!(dm, joinpath(path,"base","load.csv"))
        parse_availibility!(dm,joinpath(path,"base","availability.csv"))
        calc_base_parameters!(dm)


        parse_ev_technologies!(dm, joinpath(path,"ev","ev.csv"))
        parse_ev_demand!(dm, joinpath(path,"ev","ev_demand.csv"))
        parse_ev_power!(dm, joinpath(path,"ev","ev_power.csv"))
        calc_quantity!(dm, ev)



        parse_heat!(dm, joinpath(path,"heat","heat.csv"))
        parse_heat_technologies!(dm, joinpath(path,"heat","heat_technologies.csv"))
        parse_buildings!(dm, joinpath(path,"heat","buildings.csv"))
        parse_temperature!(dm, joinpath(path,"heat","temperature.csv"))
        parse_heat_demand!(dm, joinpath(path,"heat","heat_demand.csv"))
        parse_dhw_demand!(dm, joinpath(path,"heat","dhw_demand.csv"))
        calc_hp_cop!(dm)
        calc_heat_demand!(dm, heat)



        parse_h2_technologies!(dm, joinpath(path,"h2","h2_technologies.csv"))
        calc_inv_gas!(dm)
end
=#
