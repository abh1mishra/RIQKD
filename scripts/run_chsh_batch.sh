#!/bin/bash
# Run the CHSH scans of a job list (one line of arguments of scripts/chsh_scan.jl per scan) with N parallel workers:
#     scripts/run_chsh_batch.sh jobs/chsh.txt N [ROOT]
# Each scan writes <out>/chsh__<rel>__G<G>__<grid>.json, where <out> is the out= argument of the line (default
# results/chsh). With the optional ROOT, the leading "results" of <out> is replaced by ROOT, e.g. ROOT=rerun writes
# rerun/chsh/ and leaves the shipped results untouched. Scans whose result file already has "status": "done" are
# skipped. The log of each scan is written to logs/. The exit status is non-zero if a scan failed. The Julia
# executable can be set with the environment variable JULIA.
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
    read -r -a args <<< "$1"          # split the line into arguments, without pathname expansion
    # output directory of the scan and name of its result file, as in scripts/chsh_scan.jl
    target=$(printf '%s\n' "$ROOT" "${args[@]}" | "$JULIA" --startup-file=no -e '
        root = readline()
        d = Dict("rel" => "lean:4", "grid" => "g26", "out" => "results/chsh", "Slist" => "")
        for a in eachline(); k, v = split(a, "="; limit = 2); d[k] = v; end
        out = startswith(d["out"], "results") ? root * d["out"][8:end] : d["out"]
        src = get(d, "path", "") == "source"
        id = replace("chsh__$(d["rel"])__" * (src ? "source" : "G$(d["G"])") * "__$(d["grid"])" * (isempty(d["Slist"]) ? "" : "__test"), r"[^A-Za-z0-9_.@+-]" => "_")
        print(joinpath(out, id))')
    if grep -q '"status": "done"' "$target.json" 2>/dev/null; then
        echo "$(date '+%F %T') skip  $target (done)"
        return 0
    fi
    log="logs/$(echo "$target" | tr '/' '_').log"
    echo "$(date '+%F %T') start $target"
    "$JULIA" --project=. --startup-file=no scripts/chsh_scan.jl "${args[@]}" "out=$(dirname "$target")" > "$log" 2>&1
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
