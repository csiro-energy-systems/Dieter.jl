# This file is part of Dieter: Utility functions

"Given a Dict of keys => filenames, check at files exist"
function check_files_exist(fileDict::Dict{String,String})
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
        df = DBInterface.execute(db,queryString) |> DataFrame
    else
        error("Unrecognised filetype (as parsed from the file extension).")
    end

    return df
end

tuple2_filter(func,A,B) = filter(x->func(x[1],x[2]), [(a,b) for a in A for b in B])
# tuple3_filter(func,A,B,C) = filter(x->func(x[1],x[2],x[3]), [(a,b,c) for a in A for b in B for c in C])
# tupleN_filter(func,X...) = filter(x->func(x...), "CartesianProduct"(X...))

"""
Using a DataFrame, `create_relation` returns a function that can be used to determine
a 2-valued relation, that is, a set-to-set correspondence.
DataFrame column-names (Symbols) are given that define the set tuples,
and a third Indicator column-name is used to indicate set membership.
The return type is a function that can be used with two "sets" A and B in
`tuple2_filter(func,A,B) = filter(x->func(x[1],x[2]), [(a,b) for a in A for b in B])`
"""
function create_relation(df::DataFrame,First::Symbol,Second::Symbol,
                Indicator::Symbol) # ; criterion=testfunc_returns_bool

    if !(String(First) in names(df))
        @error "The symbol $First must be a valid name of the DataFrame."
    elseif !(String(Second) in names(df))
        @error "The symbol $Second must be a valid name of the DataFrame."
    elseif !(String(Indicator) in names(df))
        @error "The symbol $Indicator must be a valid name of the DataFrame."
    end

    # FirstType = eltype(df[!,:($First)])
    # SecondType   = eltype(df[!,:($Second)])

    func(x,y) =
        begin
            r = @where(df, :($First) .== x, :($Second) .== y)
            return sum(r[:,:($Indicator)]) > 0  ## This could be some other criteron for a relation existing:
            #  TODO allow user to pass the criterion function as an argument
        end

    return func
end

"Update a data Dict, merging the data if present, adding data if not present"
function update_dict!(dict::Dict{Symbol,Dict}, key::Symbol, val::Dict)
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
and each `Dict` has a key-value pair of id_cols => data-value.
Column symbols listed in `skip_cols` are skipped.

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

julia> map_idcol(df, 1, skip_cols=[:z])
Dict{Symbol,Dict} with 1 entry:
  :y => Dict(4=>14,2=>12,3=>13,1=>11)
```
"""
function map_idcol(df::DataFrame, id_cols::Int; skip_cols=Symbol[])

    push!(skip_cols, propertynames(df)[id_cols])

    if id_cols > 1
        ids = collect(zip([col for col in eachcol(df[:,1:id_cols])]...))
    else
        ids = df[!, 1]
    end

    non_id = df[:, id_cols+1:end]

    dict = Dict{Symbol, Dict}()

    for col in pairs(eachcol(non_id))
        if !(col[1] in skip_cols)
            vals = Dict(collect(zip(ids,col[2])))
            dict[col[1]] = vals
        end
    end
    return dict
end

map_idcol(df::DataFrame; skip_cols=Symbol[]) = map_idcol(df, 1, skip_cols=skip_cols)
map_idcol(df::DataFrame, id_col::Symbol; skip_cols=Symbol[]) = map_idcol(df, [id_col], skip_cols=skip_cols)

function map_idcol(df::DataFrame, id_cols::Array{Symbol,1}; skip_cols=Symbol[])

    for sym in skip_cols
        if sym in id_cols
            @warn "Skipped columm $sym is a member of the given id key columns.
                    It will be ignored."
        elseif !(sym in propertynames(df))
            @warn "Skipped columm $sym is a not column name in the DataFrame.
                    It will be ignored."
        end
    end

    push!(skip_cols, id_cols...)

    if length(id_cols) == 1 ## The Dict id will be a singleton
        Keys = df[!,:($(id_cols[1]))]   # Returns a Vector
    else # Make the Dict keys into Tuples
        Keys = [Tuple(row[id_cols]) for row in eachrow(df)]
    end

    non_id = setdiff(propertynames(df),skip_cols)  # Remove skipped columms

    dict = Dict{Symbol, Dict}()

    for col in pairs(eachcol(df))
        if col[1] in non_id
            dict[col[1]] = Dict(collect(zip(Keys,col[2])))
        end
    end

    return dict
end


"""
Create a `Dict` from a `DataFrame` with the `DataFrame` column-names as keys and columns as values.
"""
map_dfheader_to_col(df::DataFrame) = Dict(String(col[1]) => disallowmissing(col[2]) for col in pairs(eachcol(df)))

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

function convert_jump_container_to_df(var::VariableRef;
    dim_names::Vector{Symbol}=Vector{Symbol}(),
    value_col::Symbol=:Value)

    return convert_jump_container_to_df([var],dim_names=dim_names,value_col=value_col)
end

"""
Function to convert an Array{Any,2} with a header row of column name strings to a DataFrame
 - Thanks: https://stackoverflow.com/questions/25894634/dataframe-from-array-with-header
"""
table2df(mat) = DataFrame(mat[2:end,:], Symbol.(mat[1,:]))

"""
    tableToDict(tab::Array{Any,2}; keycols::Array{Int64,1}=[1])

Create a dictionary from an array in table form, where the table `tab` has:
 - key-index or indices given in the `keycol` array (integer column position), and
 - the column names for data appear in the first row.
"""
function tableToDict(tab::Array{<:Any,2}; keycols::Array{Int64,1}=[1])
    DataDict = Dict{Union{String,Tuple{String,Vararg{String}}},Any}() # initialise the dictionary that the function returns

    KeysName    = tab[1,keycols]      # the column label on the key column
    ColumnNames = tab[1,:]            # all the column names
    KeyArray    = tab[2:end,keycols]  # the keys for each row
    TableData   = tab[2:end,:]        # the data in the input array

    NumberRows    = size(TableData)[1]
    NumberColumns = size(TableData)[2]
    # Index columns with the key column removed:
    DataColumnIndices = setdiff(1:NumberColumns,keycols)

    if length(keycols) == 1 # # Make the DataDict key a singleton
        Keys = KeyArray
    else # Make the DataDict keys into Tuples
        Keys = [Tuple(KeyArray[i,:]) for i in 1:NumberRows]
    end

    for (k, key_value) in enumerate(Keys)
        DataDict[key_value] = Dict(ColumnNames[col] => TableData[k,col]
                                        for col in DataColumnIndices)
    end
    return DataDict
end

function SQLqueryToDict(sqlquery::SQLite.Query; keycols::Array{Int64,1}=[1])
    df = DataFrame(sqlquery)
    A = convert(Matrix, df)
    # or:
    # A = Tables.matrix(tb)s
    # Convert DataFrame to a table with column headers
    tableArray = [reshape(string.(names(df)),1,size(A)[2]); A]

    DataDict = tableToDict(tableArray; keycols=keycols)

    return DataDict
end



"""
`split_df_tuple` is a function that transforms a DataFrame with a Tuple{N,M...} type column.
It outputs an augmented DataFrame with the tuples split into two columns of type N, M... resp.

The function splits the column `InputCol` of type `Tuple` into individual columns named with `OutputCols` vector of Symbols.
Note that the `Symbol`s in `OutputCols` need to be of the name length and order to correspond to the splitting `Tuple`s.

 ### Arguments

 * `df` : a `DataFrame`
 * `InputCol` : a `Symbol` for a column of type `Tuple{Types...}`
 * `OutputCols`: a `Vector` of `Symbol`s for new column names,

 ### Returns

 * `::DataFrame`

 ### Examples

 ```
 julia> using DataFramesMeta, DataFrames

 julia> df = DataFrame(A_N = [("a",1),("b",2)], V = [8,9]);

 julia> split_df_tuple(df, :A_N, [:A, :N])
2×4 DataFrame
│ Row │ A_N      │ V     │ A      │ N     │
│     │ Tuple…   │ Int64 │ String │ Int64 │
├─────┼──────────┼───────┼────────┼───────┤
│ 1   │ ("a", 1) │ 8     │ a      │ 1     │
│ 2   │ ("b", 2) │ 9     │ b      │ 2     │
 ```

"""
function split_df_tuple(df::DataFrame, InputCol::Symbol,OutputCols::Vector{Symbol})
      return split_df_tuple!(copy(df),InputCol,OutputCols)
end

function split_df_tuple!(df::DataFrame, InputCol::Symbol,OutputCols::Vector{Symbol})
      for (i, newcol) in enumerate(OutputCols)
            df[!,newcol] = [tup[i] for tup in df[!,InputCol]]
      end
      return df
end

"Check if a `Dict`ionary has all missing values; if so, return `true`"
function isDictAllMissing(dict::Dict{T,S} where {T, S <: Union{Missing, Any}})
    return filter(x -> !ismissing(x.second), dict) |> isempty
end

# Filter a Dictionary for a matching string in the key.
dkeymatch(d::Dict,s::Regex) = filter(x -> (!ismissing(x.second) && occursin(s,x.first)),d)
# Filter a Dictionary for a matching string in the `pos`-th position in a Tuple key.
dkeymatch(d::Dict,s::Regex,pos::Int) = filter(x -> (!ismissing(x.second) && occursin(s,x.first[pos])),d)
# Filter a Dictionary for a matching string in the value.
dvalmatch(d::Dict,s::Regex) = filter(x -> (!ismissing(x.second) && occursin(s,x.second)),d)
