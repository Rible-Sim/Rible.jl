using Documenter
using DocumenterVitepress
using DocumenterCitations
using DocStringExtensions
using Literate
using DemoCards
using Rible

# Draft mode: set DOCUMENTER_DRAFT=true to skip @example execution
const draft = get(ENV, "DOCUMENTER_DRAFT", "false") == "true"

DocMeta.setdocmeta!(Rible, :DocTestSetup, :(using Rible); recursive=true)

# --- Bibliography ---
bib = CitationBibliography(joinpath(@__DIR__, "refs.bib"))

# --- Literate: generate markdown from .jl sources ---
Literate.markdown(
    joinpath(@__DIR__, "src", "en", "get_started.jl"),
    joinpath(@__DIR__, "src", "en");
    flavor = Literate.DocumenterFlavor(),
    credit = true,
)

pages = [
    "Home" => "index.md",
    "Get Started" => "get_started.md",
    "Tutorial" => [
        "Modeling Basics" => [
            "Bodies and Apparatuses" => "body_and_apparatus.md",
            "Structure" => "structure.md",
            "Control Hub" => "hub.md",
        ],
        "Robot Dynamics" => "dynamics.md",
    ],
    "Advanced" => [
        "Errors, Costs, and Sensitivity" => "adjoint/index.md",
        "Visualization" => "vis.md",
    ],
    "Extensions" => [
        "Coordinates" => [
            "Quaternion Coordinates" => "packages/RibleQCF.md",
        ],
        "Bodies and Structures" => [
            "Tensegrity Structures" => "packages/RibleTensegrity.md",
        ],
        "Integrators" => [
            "Extra Integrators" => "packages/RibleExtraIntegrators.md",
        ]
    ],
    "Development" => [
        "Natural Coordinates" => "dev/NCF.md"
    ],
    "API" => "api.md",
]

deploy_decision = Documenter.DeployDecision(;
    all_ok = true, #Should documentation be deployed?
    branch = "pages", #The branch to which documentation should be pushed
    is_preview = false, #Is this documentation build a pull request?
    repo = "github.com/Rible-Sim/Rible.jl", #The repo to which documentation should be pushed
    subfolder = "" # unused — deploy is by make_combined.sh
)

makedocs(;
    sitename = "Rible.jl",
    format = DocumenterVitepress.MarkdownVitepress(;
         repo = "github.com/Rible-Sim/Rible.jl", # The full URL of the repository
        devbranch = "dev", # The name of the development branch, like master or main.
        devurl = "dev", #The URL path to the development site, like dev or dev-branch.
        build_vitepress = get(ENV, "SKIP_VITEPRESS", "") != "true", # `false` to view locally, `true` for deploy
        deploy_url = "https://docs.rible.dev", #cannot subpath, If you are deploying from a custom URL,
        assets = "../assets", #A list of assets, the same as what is provided to Documenter's HTMLWriter.
        # keep, Sets the granularity of versions which should be kept. Options are :patch, :minor or :breaking (the default).
        # inventory_version, A version string, Defaults to the version defined in the Project.toml file in the parent folder of the documentation root
        deploy_decision,#: DeployDecision from Documenter.jl. This is used to determine whether to deploy the documentation or not. Options are: 
        # nothing: Default. Automatically determine whether to deploy the documentation.
        # Documenter.DeployDecision: Override the automatic decision and deploy based on the passed config.
    ),
    repo = Documenter.Remotes.GitHub("Rible-Sim", "Rible.jl"), # repo url
    # remotes = nothing, # to declare a list additional path::AbstractString => remote pairs that are used to determine the remote repository URLs 
    source = "src/en", #where the markdown source files are read from
    build = "build_en", #into which generated files and folders are written
    clean = true, #whether to do clean build
    doctest = true,
    modules = [Rible],
    authors = "Jiahui Luo and contributors",
    highlightsig = true,
    pagesonly = true,
    draft, # set DOCUMENTER_DRAFT=true for fast syntax-only builds
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
    # root    = "<current-directory>", the directory from which makedocs should run
    # workdir = @__DIR__, # the working directory where @example and @repl code blocks are executed. 
    # sitename = "", displayed in the title bar and/or the navigation menu when applicable.
    # expandfirst = [], allows some of the pages to be expanded (i.e. at-blocks evaluated etc.) before the others.
    # meta = Dict(:DocTestSetup => :(using MyPackages)),  can be used to provide default values for the @meta blocks executed on every page. 
    pages = pages,
    plugins = [bib],
)

# DocumenterVitepress moves rendered files from `build/final_site` into `build` on CI by default, but not when running locally

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
