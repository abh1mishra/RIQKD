#!/bin/bash
# Run every job of a job list (one line of arguments of scripts/run_job.jl per job) with N parallel workers,
# one Julia process per job:
#     scripts/run_batch.sh jobs/keyrates.txt N [ROOT]
# Each job writes <out>/<jobid>.json, where <out> is the out= argument of the job line (default results). With the
# optional ROOT, the leading "results" of <out> is replaced by ROOT, e.g. ROOT=rerun writes rerun/ and rerun/curves/
# and leaves the shipped results untouched. Jobs whose result file already has "status": "done" are skipped, so an
# interrupted batch can simply be relaunched; an unfinished job restarts from the beginning. The log of each job is
# written to logs/. The exit status is non-zero if a job failed. The Julia executable can be set with the
# environment variable JULIA.
set -u
cd "$(dirname "$0")/.."
JOBS=$1
N=${2:-4}
ROOT=${3:-results}
[ -r "$JOBS" ] || { echo "cannot read the job list $JOBS" >&2; exit 2; }
JULIA=${JULIA:-julia}
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p logs

run_one() {
    read -r -a args <<< "$1"          # split the job line into arguments, without pathname expansion
    # output directory of the job and name of its result file, as in scripts/run_job.jl
    target=$(printf '%s\n' "$ROOT" "${args[@]}" | "$JULIA" --startup-file=no -e '
        root = readline()
        d = Dict("grid" => "g26", "dG" => "0.00", "out" => "results")
        for a in eachline(); k, v = split(a, "="; limit = 2); d[k] = v; end
        out = startswith(d["out"], "results") ? root * d["out"][8:end] : d["out"]
        id = replace("$(d["prot"])__$(d["rel"])__dG$(d["dG"])__$(d["sweep"])__$(d["grid"])", r"[^A-Za-z0-9_.@+-]" => "_")
        print(joinpath(out, id))')
    if grep -q '"status": "done"' "$target.json" 2>/dev/null; then
        echo "$(date '+%F %T') skip  $target (done)"
        return 0
    fi
    log="logs/$(echo "$target" | tr '/' '_').log"
    echo "$(date '+%F %T') start $target"
    "$JULIA" --project=. --startup-file=no scripts/run_job.jl "${args[@]}" "out=$(dirname "$target")" > "$log" 2>&1
    rc=$?
    echo "$(date '+%F %T') end   $target (exit $rc)"
    return $rc
}
export -f run_one
export JULIA ROOT

grep -v '^\s*#' "$JOBS" | grep -v '^\s*$' | xargs -d '\n' -P "$N" -I{} bash -c 'run_one "$1"' _ {}
rc=$?
echo "$(date '+%F %T') batch $JOBS finished (exit $rc)"
exit $rc
