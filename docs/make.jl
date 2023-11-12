using Documenter, DocStringExtensions, Literate
using Rible

cd(@__DIR__)
function replace_includes(str)

    included = [
        "1Dfield_temporalprediction.jl",
        "2Dfield_crossprediction.jl", 
        "2Dfield_temporalprediction.jl"
    ]

    path = dirname(dirname(pathof(Rible)))*"/examples/"

    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    return str
end

# Literate it:
Literate.markdown(
    "src/stexamples.jl", 
    "src/";
    name = "stexamples", 
    preprocess = replace_includes,
    documenter = true,
)

Literate.markdown(
    "examples/tail/dynamics.jl", 
    "src/";
    name="tail"
)

#      
makedocs(
    root = @__DIR__,
    source = "src", #where the markdown source files are read from
    build = "build", #into which generated files and folders are written
    clean = true, #whether to do clean build
    doctest = true,
    modules = [Rible],
    # repo = Documenter.Remotes.GitHub("jacobleft/Rible.jl"),# concrete Remotes.Remote object
    repo = Documenter.Remotes.GitLab("robotgroup/Rible.jl"),# concrete Remotes.Remote object
    # remotes = [
    #     path => remote,
    # ],
    highlightsig = true,
    sitename = "Rible.jl",
    pages = [
        "index.md",
        "setup.md",
        "Modeling" => [
            "naturalcoordinates.md",
            "rigidbody.md",
            # "cable.md",
        ],
        # "tensegrity.md",
        # "control.md",
        # "Linearization" => [
        #     "linearization.md"
        # ],
        # "Statics" => [
        #     "inverse_statics.md"
        # ],
        "Dynamics" => [
            "solvers.md"
        ],
        # "Examples" => [
        #     # "tail.md"
        # ]
    ],
    pagesonly = true,
    draft = false,
    checkdocs = :all,
    linkcheck = true,
    warnonly = Documenter.except(
        # :autodocs_block, :cross_references, 
        # :docs_block, :doctest, 
        # :eval_block, :example_block, 
        # :footnote, 
        # :linkcheck_remotes, :linkcheck, 
        # :meta_block, :missing_docs, 
        # :parse_error, 
        # :setup_block
    ),
    workdir = joinpath(@__DIR__, "../yard"), # the working directory where @example and @repl code blocks are executed. 
    format = Documenter.HTML(),
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# deploydocs(
#     root = "<current-directory>",
#     target = "build",
#     dirname = "",
#     repo = "<required>",
#     branch = "gh-pages",
#     deps = nothing | <Function>,
#     make = nothing | <Function>,
#     devbranch = nothing,
#     devurl = "dev",
#     versions = ["stable" => "v^", "v#.#", devurl => devurl],
#     forcepush = false,
#     deploy_config = auto_detect_deploy_system(),
#     push_preview = false,
#     repo_previews = repo,
#     branch_previews = branch,
#     tag_prefix = "",
# )
