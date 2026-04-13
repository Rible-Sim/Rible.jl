#!/usr/bin/env julia

import BenchmarkTools
import PkgBenchmark
using PkgBenchmark: BenchmarkConfig, export_markdown, judge, readresults, writeresults, BenchmarkResults
import Rible
using Dates

project_root = normpath(joinpath(@__DIR__, "..", ".."))

bench_script = joinpath(@__DIR__, "benchmarks.jl")
baseline_path = joinpath(@__DIR__, "baseline.json")
tune_path = joinpath(@__DIR__, "tune.json")
json_path = joinpath(@__DIR__, "benchmark_results.json")
md_path = joinpath(@__DIR__, "benchmark_results.md")
judge_md_path = joinpath(@__DIR__, "benchmark_judge.md")

results_json = PkgBenchmark._runbenchmark(
    bench_script,
    json_path,
    BenchmarkConfig(),
    tune_path;
    runoptions = (verbose = false,),
    custom_loadpath = "",
)

results = BenchmarkResults(
    "Rible",
    results_json["juliasha"],
    BenchmarkTools.load(
        IOBuffer(results_json["results"])
    )[1],
    now(),
    results_json["juliasha"],
    results_json["vinfo"],
    BenchmarkConfig()
)

export_markdown(md_path, results)

if isfile(baseline_path)
    baseline = readresults(baseline_path)
    judgement = judge(results, baseline)
    export_markdown(judge_md_path, judgement; export_invariants = true)
    if BenchmarkTools.isregression(judgement)
        @warn "Benchmark regression detected. See benchmark_judge.md for details."
        exit(1)
    end
else
    @warn "No baseline found. Writing current results to $(basename(baseline_path))."
    writeresults(baseline_path, results)
end

println("Benchmark results saved to:")
println("  $md_path")
println("  $json_path")
isfile(judge_md_path) && println("Benchmark comparison saved to: $judge_md_path")
