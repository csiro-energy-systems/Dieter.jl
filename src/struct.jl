mutable struct DieterModel
    sets::Dict{Symbol, Array{String,1}}
    parameters::Dict{Symbol, Any}
    settings::Dict{Symbol, Any}
    results::Dict{Symbol, DataFrame}

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

        parse_base_technologies!(self, joinpath(path,"base","technologies.csv"))
        parse_storages!(self, joinpath(path,"base","storages.csv"))
        parse_load!(self, joinpath(path,"base","load.csv"))
        parse_availibility!(self,joinpath(path,"base","availability.csv"))
        calc_base_parameters!(self)


        parse_ev_technologies!(self, joinpath(path,"ev","ev.csv"))
        parse_ev_demand!(self, joinpath(path,"ev","ev_demand.csv"))
        parse_ev_power!(self, joinpath(path,"ev","ev_power.csv"))
        calc_quantity!(self, ev)



        parse_heat!(self, joinpath(path,"heat","heat.csv"))
        parse_heat_technologies!(self, joinpath(path,"heat","heat_technologies.csv"))
        parse_buildings!(self, joinpath(path,"heat","buildings.csv"))
        parse_temperature!(self, joinpath(path,"heat","temperature.csv"))
        parse_heat_demand!(self, joinpath(path,"heat","heat_demand.csv"))
        parse_dhw_demand!(self, joinpath(path,"heat","dhw_demand.csv"))
        calc_hp_cop!(self)
        calc_heat_demand!(self, heat)



        parse_h2_technologies!(self, joinpath(path,"h2","h2_technologies.csv"))
        calc_inv_gas!(self)


        return self
    end

    DieterModel() = new()
end
