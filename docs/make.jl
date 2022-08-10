using Documenter
using DocStringExtensions
using TensegrityRobots

makedocs(
    sitename = "张拉整体机器人",
    format = Documenter.HTML(),
    modules = [TensegrityRobots],
    workdir = @__DIR__,
    pages = [
        "index.md",
        "setup.md",
        "建模基础" => [
            "naturalcoordinates.md",
            "rigidbody.md",
            "cable.md",
        ],
        "tensegrity.md",
        "control.md",
        "分析方法基础" => [
            "linearization.md"
        ],
        "静力学分析" => [
            "inverse_statics.md"
        ],
        "动力学分析" => [
            "solvers.md"
        ],
        "例子" => [
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
