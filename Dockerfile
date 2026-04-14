# GitLab CI-optimized Dockerfile for Rible.jl
# - Uses BuildKit cache mounts so dependency downloads and precompilation persist between pipeline runs.
# - Copies the prepared depot into the final image so CI jobs can reuse precompiled packages.

ARG JULIA_VERSION=1.12
# Branch-specific cache key; override in CI with CACHE_SUFFIX=$CI_COMMIT_REF_SLUG
ARG CACHE_SUFFIX=ci

FROM docker.m.daocloud.io/library/julia:${JULIA_VERSION} AS base

WORKDIR /app
# Default depot path inside the running container
ENV JULIA_DEPOT_PATH=/app/depot

# -----------------------------------------------------------------------------
# deps: fetch & precompile with cached depot
# -----------------------------------------------------------------------------
FROM base AS deps

# Install CI tools into the image layer so they are always available
RUN julia --color=yes -e 'import Pkg; Pkg.add(["TestEnv", "LocalCoverage", "CoverageTools", "PkgBenchmark"]); Pkg.precompile()'

ARG CACHE_SUFFIX
# Build into a cached depot directory; we'll copy it into the image after precompiling
# Stack the cache depot on top of the image depot so we can use installed tools
ENV JULIA_DEPOT_PATH=/cache/julia-depot:/app/depot

# Copy only project metadata first to maximize cache hits for dependency downloads
COPY Project.toml ./
COPY RibleTensegrity/Project.toml RibleTensegrity/
COPY RibleQCF/Project.toml RibleQCF/
COPY RibleExtraIntegrators/Project.toml RibleExtraIntegrators/
COPY docs/Project.toml docs/

# Preinstall LocalCoverage in the depot for CI coverage jobs without mutating project manifests
# Preinstall LocalCoverage in the depot for CI coverage jobs without mutating project manifests
RUN --mount=type=cache,id=julia-depot-${JULIA_VERSION}-${CACHE_SUFFIX},target=/cache/julia-depot,sharing=locked \
    julia --color=yes -e 'using InteractiveUtils; versioninfo()'

# Create dummy source files so Pkg.instantiate() can resolve local paths without full source
# This allows caching of external dependencies before the actual source code is copied
RUN mkdir -p src \
    RibleTensegrity/src \
    RibleQCF/src \
    RibleExtraIntegrators/src \
    ext/RibleRungeKuttaExt \
    && echo 'module Rible end' > src/Rible.jl \
    && echo 'module RibleTensegrity end' > RibleTensegrity/src/RibleTensegrity.jl \
    && echo 'module RibleQCF end' > RibleQCF/src/RibleQCF.jl \
    && echo 'module RibleExtraIntegrators end' > RibleExtraIntegrators/src/RibleExtraIntegrators.jl \
    && echo 'module RibleGLMakieExt end' > ext/RibleGLMakieExt.jl \
    && echo 'module RibleApproxFunExt end' > ext/RibleApproxFunExt.jl \
    && echo 'module RibleBSplineKitExt end' > ext/RibleBSplineKitExt.jl \
    && echo 'module RibleRungeKuttaExt end' > ext/RibleRungeKuttaExt/RibleRungeKuttaExt.jl

# Warm the package resolver/download cache (no source yet, so this layer stays valid)
RUN --mount=type=cache,id=julia-depot-${JULIA_VERSION}-${CACHE_SUFFIX},target=/cache/julia-depot,sharing=locked \
    set -eux; \
    julia  --project=.                     -e 'using TestEnv; TestEnv.activate(); import Pkg; Pkg.instantiate()'; \
    julia  --project=RibleTensegrity       -e 'using TestEnv; TestEnv.activate(); import Pkg; Pkg.instantiate()'; \
    julia  --project=RibleQCF              -e 'using TestEnv; TestEnv.activate(); import Pkg; Pkg.instantiate()'; \
    julia  --project=RibleExtraIntegrators -e 'using TestEnv; TestEnv.activate(); import Pkg; Pkg.instantiate()'; \
    julia  --project=docs -e 'import Pkg; Pkg.instantiate()'


# Bring in source code needed for precompilation and tests
COPY src src
COPY ext ext
COPY test test
COPY RibleTensegrity/src RibleTensegrity/src
COPY RibleTensegrity/test RibleTensegrity/test
COPY RibleQCF/src RibleQCF/src
COPY RibleQCF/test RibleQCF/test
COPY RibleExtraIntegrators/src RibleExtraIntegrators/src
COPY RibleExtraIntegrators/test RibleExtraIntegrators/test
COPY docs/src docs/src
COPY docs/make*.jl docs/

# Precompile into the cached depot and copy it into place for the runtime image
RUN --mount=type=cache,id=julia-depot-${JULIA_VERSION}-${CACHE_SUFFIX},target=/cache/julia-depot,sharing=locked \
    set -eux; \
    julia --code-coverage=@                      --color=yes --check-bounds=yes --warn-overwrite=yes --depwarn=yes --inline=yes --startup-file=no --track-allocation=none --project=. -e 'using TestEnv; TestEnv.activate(); import Pkg; Pkg.precompile()'; \
    julia --code-coverage=@RibleTensegrity       --color=yes --check-bounds=yes --warn-overwrite=yes --depwarn=yes --inline=yes --startup-file=no --track-allocation=none --project=RibleTensegrity -e 'using TestEnv; TestEnv.activate(); import Pkg; Pkg.precompile()'; \
    julia --code-coverage=@RibleQCF              --color=yes --check-bounds=yes --warn-overwrite=yes --depwarn=yes --inline=yes --startup-file=no --track-allocation=none --project=RibleQCF -e 'using TestEnv; TestEnv.activate(); import Pkg; Pkg.precompile()'; \
    julia --code-coverage=@RibleExtraIntegrators --color=yes --check-bounds=yes --warn-overwrite=yes --depwarn=yes --inline=yes --startup-file=no --track-allocation=none --project=RibleExtraIntegrators -e 'using TestEnv; TestEnv.activate(); import Pkg; Pkg.precompile()'; \
    julia --color=yes --project=docs -e 'import Pkg; Pkg.precompile()'; \
    cp -a /cache/julia-depot/. /app/depot; \
    rm -rf /app/depot/logs /app/depot/downloads

# -----------------------------------------------------------------------------
# runtime: reuse deps stage filesystem to avoid extra copies
# -----------------------------------------------------------------------------
FROM deps AS runtime

# Set runtime depot path
ENV JULIA_DEPOT_PATH=/app/depot

# Default to a status command; CI jobs can override
CMD ["julia", "--project=.", "-e", "using Pkg; Pkg.status()"]
