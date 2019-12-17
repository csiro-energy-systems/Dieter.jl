# This file is part of Dieter: Utility functions

"Update a data Dict, merging the data if present, adding data if not present"
function update_dict!(dict::Dict{Symbol,Any}, key, val::Dict)
    # if !isa(val,Dict)
        # @warn "This function expects the value $val to be of type Dict"
    # end
    if haskey(dict, key)
        merge!(dict[key], val) # Error if dict[key] & val are not Dict types
        # merge!(dict, Dict(key => val))
    else
        push!(dict, key => val)
    end
end

"""
Return a dictionary of `Dict`s from a `DataFrame` with column-names as keys,
and each `Dict` has a key-value pair of id => data-value.
Column symbols listed in `skip` are skipped.

# Examples

```julia
julia> df = DataFrame(x=1:4, y=11:14, z=2:5)
4×3 DataFrame
│ Row │ x     │ y     │ z     │
│     │ Int64 │ Int64 │ Int64 │
├─────┼───────┼───────┼───────┤
│ 1   │ 1     │ 11    │ 2     │
│ 2   │ 2     │ 12    │ 3     │
│ 3   │ 3     │ 13    │ 4     │
│ 4   │ 4     │ 14    │ 5     │

julia> map_idcol(df, skip=[:z], id_cols=1)
Dict{Symbol,Dict} with 1 entry:
  :y => Dict(4=>14,2=>12,3=>13,1=>11)
```
"""
function map_idcol(df::DataFrame; id_cols::Int=1, skip=Symbol[])
    push!(skip, names(df)[id_cols])

    if id_cols > 1
        ids = zip([col for col in eachcol(df[:,1:id_cols], false)]...) |>
            collect
    else
        ids = df[!, 1]
    end

    non_id = df[:, id_cols+1:end]
    dict = Dict{Symbol, Dict}()
    for col in eachcol(non_id, true)
        if !(col[1] in skip)
            vals = Dict(collect(zip(ids,col[2])))
            dict[col[1]] = vals
        end
    end
    return dict
end

"""
Create a `Dict` from a `DataFrame` with the `DataFrame` column-names as keys and columns as values.
"""
map_dfheader_to_col(df::DataFrame) = Dict(String(col[1]) => disallowmissing(col[2]) for col in eachcol(df, true))

get_variable_names(m::JuMP.Model) = [v.name for (k,v) in m.varData]  # varData requires instantiation for JuMP Model


"""
Returns a `DataFrame` with the values of the variables from the JuMP container `var`.
The column names of the `DataFrame` can be specified for the indexing columns in `dim_names`,
and the name of the data value column by a Symbol `value_col` e.g. :Value
"""
function convert_jump_container_to_df(var::JuMP.Containers.DenseAxisArray;
    dim_names::Vector{Symbol}=Vector{Symbol}(),
    value_col::Symbol=:Value)

    (length(dim_names) == 0) && (dim_names = [Symbol("dim$i") for i in 1:length(var.axes)])

    if length(dim_names) != length(var.axes)
        throw(ArgumentError("Length of given name list does not fit the number of variable dimensions"))
    end

    tup_dim = (dim_names...,)

    # With a product over all axis sets of size M, form an Mx1 Array of all indices to the JuMP container `var`
    ind = reshape([collect(k[i] for i in 1:length(dim_names)) for k in Base.Iterators.product(var.axes...)],:,1)

    var_val  = value.(var)

    df = [merge(NamedTuple{tup_dim}(ind[i]), NamedTuple{(value_col,)}(var_val[(ind[i]...,)...]))
        for i in 1:length(ind)] |>
        DataFrame

    return df
end



function convert_jump_container_to_df(var::JuMP.Containers.SparseAxisArray;
    dim_names::Vector{Symbol}=Vector{Symbol}(),
    value_col::Symbol=:Value)

    var_data = (value.(var)).data

    num_dim = length(collect(keys(var_data))[1])  # Roundabout way of finding the dimension of the var array

    (length(dim_names) == 0) && (dim_names = [Symbol("dim$i") for i in 1:num_dim])

    if length(dim_names) != num_dim
        throw(ArgumentError("Length of given name list does not fit the number of variable dimensions"))
    end

    tup_dim = (dim_names..., value_col)

    df = [NamedTuple{tup_dim}(([k[i] for i in 1:length(dim_names)]..., v)) for (k,v) in var_data] |>
        DataFrame

    return df
end

function convert_jump_container_to_df(var::Array{VariableRef};
    dim_names::Vector{Symbol}=Vector{Symbol}(),
    value_col::Symbol=:Value)

    num_dim = size(var,2)

    (length(dim_names) == 0) && (dim_names = [Symbol("dim$i") for i in 1:num_dim])

    if length(dim_names) != num_dim
        throw(ArgumentError("Length of given name list does not fit the number of variable dimensions"))
    end

    var_val = value.(var)

    df = DataFrame([var_val[:,i] for i in 1:num_dim], dim_names)

    return df
end
