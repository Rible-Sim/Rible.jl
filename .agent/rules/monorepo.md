---
trigger: always_on
---

# Agent Guidelines for Rible.jl

This document outlines the standards and workflows for AI agents working on the Rible.jl monorepo.

## Basic Rules
- Explain before editting
- Request for permission before creating files
- Ask for current setup if uncertain

## Project Structure

Rible.jl is a monorepo containing the main `Rible` package and several subpackages located in subdirectories (e.g., `RibleQCF`, `RibleTensegrity`). Subpackages often depend on the main `Rible` package via relative paths in their `Project.toml` files (under `[sources]` section).

### Folder structures
The main pkg and subpkgs has similar folder structures. Other than standard ones ( `src/`, `test/`, `docs/`, etc), there are `examples/` under which locate `bodies/`, `robots/`, `assets/` and `cases`, whose names explain themselves. 

## Development Workflows
AI agents should either engage in interactive developments OR reproduction  

### Interactive developments
Most of the time, the user has a running Julia REPL with `Revise` pkg loaded to automatically apply code changes. He is either working on a activated folder under `cases/` or one of the `tests/` files using `TestEnv` to setup to test env. Therefore, AI agents should examine and editting files therearound, rather creating files elsewhere. 

### Reproduction Scripts
When creating a reproduction script (repro) to isolate a bug or demonstrate a feature upon user request, it is critical to correctly set up the environment to use the local versions of the packages and other dependencies.

1.  **Self-Contained**: The script should handle its own environment setup.
2.  **Location and Environment**: The script should be located under `examples/repro` with its own folder, such as `examples/repro/issue_name`, which contains its own environment `Project.toml` to avoid polluting the global or project environment.
3.  **Local Dependencies**: Use `[sources]` in `Project.toml` to load the local packages using relative paths. This ensures the repro runs against the current codebase, not registered versions.
4.  **Other Dependencies**: Existing examples can be used for references as how to properly setup the repro script. These examples can be found in `examples/cases` or `{RibleSubPkg}/examples/cases`. You may copy the `[deps]` section from any of these examples to your repro folder to replicate their dependencies.
5.  **Running the Repro**: Instantiate (or update and resolve, if needed) the repro project environment
```bash
julia --project=examples/repro/issue_name -e "import Pkg; Pkg.instantiate()"
```
before running the script from the its own folder:
```bash
julia --project=examples/repro/issue_name repro_issue_name.jl
```
6. **Finishing**: Delete the repro folder and script after you are done with it.