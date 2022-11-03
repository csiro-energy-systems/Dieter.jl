## File based in PowerSystems.jl documentation `make.jl` file.

using Documenter, Dieter
import DataStructures: OrderedDict
using Literate

pages = OrderedDict(
        "Welcome Page" => "index.md",
        "Quick Start Guide" => "quick_start_guide.md",
        "Tutorials" =>  
            Any[
                # "tutorials/intro_page.md",
                "tutorials/optmodel.md"
            ],
        # "Model Developer Guide" =>
        #     Any[
        #     ],
        # "Model Library" => Any[],
        "Public API Reference" => "api/public.md",
        "Internal API Reference" => "api/internal.md"
)

#= ## Commented out until we use a sub-folder structure.

# postprocess function to insert md
function insert_md(content)
    m = match(r"APPEND_MARKDOWN\(\"(.*)\"\)", content)
    if !isnothing(m)
        md_content = read(m.captures[1], String)
        content = replace(content, r"APPEND_MARKDOWN\(\"(.*)\"\)" => md_content)
    end
    return content
end

# This code performs the automated addition of Literate - Generated Markdowns. The desired
# section name should be the name of the file for instance network_matrices.jl -> Network Matrices
julia_file_filter = x -> occursin(".jl", x)
folders = Dict(
    "Model Library" => filter(julia_file_filter, readdir("docs/src/model_library")),
    "Modeler Guide" => filter(julia_file_filter, readdir("docs/src/modeler_guide")),
    "Model Developer Guide" => filter(julia_file_filter, readdir("docs/src/model_developer_guide")),
    "Code Base Developer Guide" => filter(julia_file_filter, readdir("docs/src/code_base_developer_guide")),
)

for (section, folder) in folders
    for file in folder
        section_folder_name = lowercase(replace(section, " " => "_"))
        outputdir = joinpath(pwd(), "docs", "src", "$section_folder_name")
        inputfile = joinpath("$section_folder_name", "$file")
        infile_path = joinpath(pwd(), "docs", "src", inputfile)
        outputfile = string("generated_", replace("$file", ".jl" => ""))
        execute = occursin("EXECUTE = TRUE", uppercase(readline(infile_path))) ? true : false
        execute && include(infile_path)
        Literate.markdown(infile_path,
                          outputdir;
                          name = outputfile,
                          credit = false,
                          flavor = Literate.DocumenterFlavor(),
                          documenter = true,
                          postprocess = insert_md,
                          execute = execute)
        subsection = titlecase(replace(split(file, ".")[1], "_" => " "))
        push!(pages[section], ("$subsection" =>  joinpath("$section_folder_name", "$(outputfile).md")))
    end
end
=#
makedocs(
    modules = [Dieter],
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),),
    sitename = "Dieter.jl",
    authors = "Alexander Zerrahn, Wolf-Peter Schill, Mario Kendziorski, James Foster",
    pages = [p for p in pages]
)

deploydocs(
    target = "build",
    # dirname = "",
    repo = "github.com/csiro-energy-systems/Dieter.jl.git",
    branch = "gh-pages",
    devbranch = "main",
    devurl = "dev",
    # versions = ["v#.#"],
    # versions = ["stable" => "v^", "v#.#"],
    # deploy_config = Documenter.GitHubActions(),
    push_preview=false,
)
