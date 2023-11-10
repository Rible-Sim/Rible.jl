using Documenter
using DocStringExtensions
using Rible

makedocs(
    sitename = "Rible",
    format = Documenter.HTML(),
    modules = [Rible],
    workdir = @__DIR__,
    warnonly=:doctest,
    pages = [
        "index.md",
        "setup.md",
        "Modeling" => [
            "naturalcoords.md",
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
        "Examples" => [
            "tail.md"
        ]
    ]
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
