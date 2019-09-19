
function map_idcol(df::DataFrame; id_cols::Int=1, skip=Symbol[])
    push!(skip, names(df)[id_cols])

    if id_cols > 1
        ids = zip([col for col in eachcol(df[1:id_cols], false)]...) |>
            collect
    else
        ids = df[1]
    end

    non_id = df[id_cols+1:end]
    dict = Dict{Symbol, Dict}()
    for col in eachcol(non_id, true)
        if !(col in skip)
            vals = Dict(collect(zip(ids,col[2])))
            dict[col[1]] = vals
        end
    end
    return dict
end

map_dfheader_to_col(df::DataFrame) = Dict(String(col[1]) => disallowmissing(col[2])
    for col in eachcol(df, true))

function update_dict!(dict::Dict, key, val)
    if haskey(dict, key)
        merge!(dict[key], val)
    else
        push!(dict, key => val)
    end
end

get_variable_names(m::JuMP.Model) = [v.name for (k,v) in m.varData]



function convert_jump_container_to_df(var::JuMP.Containers.DenseAxisArray;
    dim_names::Vector{Symbol}=Vector{Symbol}(),
    value_col::Tuple{Symbol}=(:Value,))

    (length(dim_names) == 0) && (dim_names = [Symbol("x$i") for i in 1:length(var.axes)])

    if length(dim_names) != length(var.axes)
        throw(ArgumentError("Length of given name list does not fit the number of variable dimensions"))
    end

    tup_dim = (dim_names...,)

    ind = [collect(k[i] for i in 1:length(dim_names))
        for k in Base.Iterators.product(var.axes...)][:]

    df = [merge(NamedTuple{tup_dim}(ind[i]), NamedTuple{value_col}(var[(ind[i]...,)...]))
        for i in 1:length(ind)] |>
        DataFrame

    return df
end



function convert_jump_container_to_df(var::JuMP.Containers.SparseAxisArray;
    dim_names::Vector{Symbol}=Vector{Symbol}(),
    value_col::Tuple{Symbol}=(:Value,))

    len_dim = length(dim_names)
    tup_dim = (dim_names..., value_col...)

    df = [NamedTuple{tup_dim}(([k[i] for i in 1:len_dim]..., v))
        for (k,v) in var.data] |>
        DataFrame

    return df
end
