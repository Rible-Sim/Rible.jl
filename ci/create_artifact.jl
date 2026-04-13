using Pkg.Artifacts, Pkg.PlatformEngines, SHA

git_tree_sha1 = create_artifact() do dir
    # get list of Tier 2 files
    files = split(read(`git diff --name-only --cached`, String), "\n", keepempty=false)
    for f in files
        rel_f = replace(f, "examples/assets/" => "")
        mkpath(dirname(joinpath(dir, rel_f)))
        cp(f, joinpath(dir, rel_f))
    end
end

# Update ASSETS_VERSION when creating a new artifact release
const ASSETS_VERSION = "assets-v3"
const TARBALL_NAME = "rible_extra_meshes.tar.gz"

PlatformEngines.package(artifact_path(git_tree_sha1), TARBALL_NAME)

sha256_hash = open(TARBALL_NAME, "r") do io
    bytes2hex(sha256(io))
end

bind_artifact!(
    "Artifacts.toml",
    "rible_extra_meshes",
    git_tree_sha1;
    download_info = [(
        "https://github.com/Rible-Sim/Rible.jl/releases/download/$(ASSETS_VERSION)/$(TARBALL_NAME)",
        sha256_hash
    )],
    lazy = true,
    force = true
)
