using Documenter, DocStringExtensions, Literate, DemoCards
using Rible

cd(@__DIR__)

# 1. generate demo files
# demopage, postprocess_cb, demo_assets = makedemos(""; # this is the relative path to docs/
#     root = @__DIR__,
#     src = "src",
#     build = "build",
#     branch = "gh-pages",
#     edit_branch = "main",
#     credit = true,
#     throw_error = false,
# )
# if there are generated css assets, pass it to Documenter.HTML
# assets = []
# isnothing(demo_assets) || (push!(assets, demo_assets))

# 2. normal Documenter usage
makedocs(
    root = @__DIR__,
    source = "src", #where the markdown source files are read from
    build = "build", #into which generated files and folders are written
    clean = true, #whether to do clean build
    doctest = true,
    modules = [Rible],
    repo = Documenter.Remotes.GitHub("Rible-Sim/Rible.jl"),# concrete Remotes.Remote object
    # remotes = [
    #     path => remote,
    # ],
    highlightsig = true,
    sitename = "Rible.jl",
    pages = [
        "index.md",
        "setup.md",
        ## "Modeling" => [
        ##     "naturalcoordinates.md",
        ##     "rigidbody.md",
        ##     # "cable.md",
        ## ],
        # "tensegrity.md",
        # "control.md",
        # "Linearization" => [
        #     "linearization.md"
        # ],
        # "Statics" => [
        #     "inverse_statics.md"
        # ],
        ## "Demos" => [
        ##     "pointmass.md",
        ## ],
        ## "Dynamics" => [
        ##     "solvers.md"
        ## ],
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
    ## workdir = joinpath(@__DIR__, ""), # the working directory where @example and @repl code blocks are executed. 
    format = Documenter.HTML(
        # mathengine = MathJax3(
        #     Dict(
        #         :loader => Dict("load" => ["[tex]/physics"]),
        #         :tex => Dict(
        #             "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
        #             "tags" => "ams",
        #             "packages" => ["base", "ams", "autoload", "physics"],
        #         ),
        #     )
        # )
    ),
)
# 3. postprocess after makedocs
# postprocess_cb()

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
