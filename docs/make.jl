using Documenter
using DocStringExtensions
using Rible

makedocs(
    sitename = "Rible",
    format = Documenter.HTML(),
    modules = [Rible],
    workdir = @__DIR__,
    pages = [
        "index.md",
        "setup.md",
        "Modeling" => [
            "naturalcoords.md",
            "rigidbody.md",
            "cable.md",
        ],
        "tensegrity.md",
        "control.md",
        "Linearization" => [
            "linearization.md"
        ],
        "Statics" => [
            "inverse_statics.md"
        ],
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
#=deploydocs(
    repo = "<repository url>"
)=#
