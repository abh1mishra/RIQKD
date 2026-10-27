# One job: a key-rate curve and its thresholds, for one protocol, relaxation, bound G and type of noise.
#
#   julia --project=. scripts/run_job.jl prot=43 rel=lean:3 dG=0.00 sweep=v [grid=g26] [step=0.01] [tol=5e-4]
#                                        [xmax=1.0] [xmin=0.5] [threads=1] [out=results]
#
#   prot    protocol (src/protocols.jl): 42, 43, 44, 6sym, 3st_ZBB, 63, BB84+1, 3st_ZB, 3st_ZXB, 3st_ZX, 42@φ, 43@φ
#   rel     relaxation (src/relaxation.jl), e.g. lean:3, lean:2+B*B*E, lean:2s+B*B*E
#   dG      G = G_honest + dG, where G_honest is the guessing probability of the target states;
#           dG=source sets G to the guessing probability of the depolarised source (sweep=v only)
#   sweep   v (visibility, at η = 1) or eta (detection efficiency, at v = 1)
#   grid    entropy partition (src/entropy.jl): g26, uniform10 or fine120
#
# The rate is evaluated at x = xmax, xmax − step, xmax − 2 step, … down to the first point where the von Neumann
# rate is ≤ 0. If the rate is already ≤ 0 at xmax < 1, the sweep steps upwards instead. The zero crossings of the
# von Neumann and of the min-entropy rates are then bisected to `tol`. A min-entropy threshold is found only if
# the sweep crosses it. If a program fails during a bisection, the bisection stops there, its bracket is kept and
# the point is recorded in `bisection_stopped_vN` or `bisection_stopped_min`. The result is written to
# <out>/<jobid>.json after every point.

include(joinpath(@__DIR__, "..", "src", "RIQKDNumerics.jl"))
using .RIQKDNumerics, JSON, Dates, Printf

function parse_args(args)
    d = Dict{String,String}("grid" => "g26", "step" => "0.01", "tol" => "5e-4", "xmin" => "0.5", "xmax" => "1.0",
        "threads" => "1", "out" => "results", "dG" => "0.00")
    for a in args
        k, v = split(a, '='; limit = 2)
        d[k] = v
    end
    return d
end

jobid(d) = replace("$(d["prot"])__$(d["rel"])__dG$(d["dG"])__$(d["sweep"])__$(d["grid"])",
    r"[^A-Za-z0-9_.@+-]" => "_")

function main(args)
    d = parse_args(args)
    p = protocol(d["prot"])
    rel = parse_relaxation(d["rel"])
    sweep = d["sweep"]
    sweep in ("v", "eta") || error("sweep must be v or eta")
    α, β = ks_coeffs(GRIDS[d["grid"]])
    step, tol, xmin = parse(Float64, d["step"]), parse(Float64, d["tol"]), parse(Float64, d["xmin"])
    xmax = parse(Float64, d["xmax"])
    Ghon = honest_G(p)
    Gof(x) = d["dG"] == "source" ? (sweep == "v" ? honest_G(p; v = x) : Ghon) : Ghon + parse(Float64, d["dG"])
    mkpath(d["out"])
    id = jobid(d)
    outfile = joinpath(d["out"], id * ".json")

    S = build_sdp(p, rel; threads = parse(Int, d["threads"]))
    @printf("[%s] %s built: %d words, blocks %s, %d moments, %.0f s\n", now(), id, S.nwords,
        string(unique(S.blocksizes)), S.nmoments, S.build_seconds)
    flush(stdout)

    res = Dict{String,Any}(
        "job" => d, "jobid" => id, "protocol" => p.name, "description" => p.description,
        "theta_x" => p.θx, "theta_y" => p.θy, "keypairs" => p.keypairs, "pkey" => pkey(p),
        "G_honest" => Ghon, "relaxation" => string(rel), "nwords" => S.nwords, "blocksizes" => S.blocksizes,
        "nmoments" => S.nmoments, "build_seconds" => S.build_seconds, "grid" => GRIDS[d["grid"]],
        "merged_k0" => true, "host" => gethostname(), "started" => string(now()), "versions" => versions(),
        "curve" => Any[], "bisection" => Any[], "bisection_minentropy" => Any[], "status" => "running")
    # JSON has no NaN or Inf: failed programs are written as null
    clean(x::AbstractFloat) = isfinite(x) ? x : nothing
    clean(x::AbstractDict) = Dict(k => clean(v) for (k, v) in x)
    clean(x::AbstractVector) = [clean(v) for v in x]
    clean(x::Tuple) = [clean(v) for v in x]
    clean(x) = x
    save() = begin
        open(outfile * ".tmp", "w") do io
            JSON.print(io, clean(res), 1)
        end
        mv(outfile * ".tmp", outfile; force = true)
    end

    at(x) = sweep == "v" ? (v = x, η = 1.0) : (v = 1.0, η = x)
    full(x) = keyrate_point(S, α, β; at(x)..., G = Gof(x))
    minonly(x) = minentropy_point(S; at(x)..., G = Gof(x))
    record!(list, x, r) = push!(list, merge(Dict{String,Any}("x" => x), Dict(String(k) => v for (k, v) in pairs(r))))
    report(tag, x, r) = begin
        @printf("[%s] %s %s x=%.4f H=%.5f Hmin=%.5f HAB=%.5f K=%+.5f Kmin=%+.5f (%.0f s)%s\n", now(), id, tag, x,
            r.H, r.Hmin, r.HAB, r.K, r.Kmin, r.seconds,
            isempty(r.nonoptimal) ? "" : " nonoptimal=$(length(r.nonoptimal)):$(r.nonoptimal[1])")
        flush(stdout)
    end

    # sweep from x = xmax down to the first point with K ≤ 0
    last_pos = first_neg = last_pos_min = first_neg_min = nothing
    xs = collect(xmax:-step:xmin)
    for x in xs
        r = full(x)
        record!(res["curve"], x, r)
        report("curve", x, r)
        save()
        if isfinite(r.Kmin) && first_neg_min === nothing
            r.Kmin > 0 ? (last_pos_min = x) : (last_pos_min !== nothing && (first_neg_min = x))
        end
        if isfinite(r.K)
            if r.K > 0
                last_pos = x
            else
                first_neg = x
                break
            end
        end
    end

    # threshold above the starting point: step upwards until the rate is positive
    if last_pos === nothing && first_neg == xs[1] && xs[1] < 1.0
        x = xs[1]
        while x + step <= 1.0 + 1e-9
            x = min(x + step, 1.0)
            r = full(x)
            record!(res["curve"], x, r)
            report("up", x, r)
            save()
            isfinite(r.K) || continue          # a failed point carries no information on the sign of the rate
            if r.K > 0
                last_pos = x
                break
            end
            first_neg = x
        end
    end

    # Bisection, assuming that the rate decreases with the noise: rate(hi) > 0 ≥ rate(lo). A point where a
    # program failed carries no information on the sign of the rate; the bisection then stops, and the bracket
    # [lo, hi] of the last two successful evaluations is kept.
    function bisect(lo, hi, f, which, list)
        while hi - lo > tol
            mid = (lo + hi) / 2
            r = f(mid)
            record!(res[list], mid, r)
            report(string(which), mid, r)
            save()
            val = which == :vN ? r.K : r.Kmin
            if !isfinite(val)
                res["bisection_stopped_$(which)"] = mid
                break
            end
            val > 0 ? (hi = mid) : (lo = mid)
        end
        return lo, hi
    end
    res["threshold"] = res["threshold_minentropy"] = nothing
    if last_pos !== nothing && first_neg !== nothing
        lo, hi = bisect(first_neg, last_pos, full, :vN, "bisection")
        res["threshold"], res["threshold_bracket"] = (lo + hi) / 2, [lo, hi]
    end
    if last_pos_min !== nothing && first_neg_min !== nothing
        lo, hi = bisect(first_neg_min, last_pos_min, minonly, :min, "bisection_minentropy")
        res["threshold_minentropy"], res["threshold_minentropy_bracket"] = (lo + hi) / 2, [lo, hi]
    end
    # true if no evaluated point certified a positive von Neumann rate
    res["no_key_at_ideal_point"] = last_pos === nothing
    res["status"] = "done"
    res["finished"] = string(now())
    save()
    @printf("[%s] %s DONE threshold=%s threshold_minentropy=%s\n", now(), id, string(res["threshold"]),
        string(res["threshold_minentropy"]))
end

main(ARGS)
