# Key rate at fixed noise as a function of the relaxation δG of the bound, G = G_honest + δG (paper, Sec. VI.B):
#
#   julia --project=. scripts/dG_scan.jl prot=43 rel=lean:3 [v=1.0] [eta=1.0] [step=0.01] [dGmax=0.3] [tol=5e-4]
#                                         [Ktol=1e-4] [grid=g26] [out=results/dG]
#
# δG is increased in steps of `step` from 0 until the von Neumann rate K = H(A|E) − H(A|B) is at most Ktol, and the
# largest δG with K > Ktol is then bisected to `tol`; the same is done for the min-entropy rate. Rates below Ktol are
# treated as zero: where the exact value of the relaxation is zero, MOSEK returns about 10⁻⁷ bits. At v = η = 1 this
# gives the largest relaxation of the bound for which the ideal behaviour still certifies key. Output:
# <out>/dG__<prot>__<rel>__v<v>__eta<eta>__<grid>.json.

include(joinpath(@__DIR__, "..", "src", "RIQKDNumerics.jl"))
using .RIQKDNumerics, JSON, Dates, Printf

function main(args)
    d = Dict{String,String}("v" => "1.0", "eta" => "1.0", "step" => "0.01", "dGmax" => "0.3", "tol" => "5e-4",
        "Ktol" => "1e-4", "grid" => "g26", "out" => "results/dG")
    for a in args
        k, x = split(a, '='; limit = 2)
        d[k] = x
    end
    p = protocol(d["prot"])
    v, η = parse(Float64, d["v"]), parse(Float64, d["eta"])
    step, dGmax, tol = parse(Float64, d["step"]), parse(Float64, d["dGmax"]), parse(Float64, d["tol"])
    Ktol = parse(Float64, d["Ktol"])
    S = build_sdp(p, parse_relaxation(d["rel"]))
    α, β = ks_coeffs(GRIDS[d["grid"]])
    Ghon = honest_G(p)
    mkpath(d["out"])
    id = replace("dG__$(d["prot"])__$(d["rel"])__v$(d["v"])__eta$(d["eta"])__$(d["grid"])", r"[^A-Za-z0-9_.@+-]" => "_")
    outfile = joinpath(d["out"], id * ".json")
    res = Dict{String,Any}("job" => d, "protocol" => p.name, "relaxation" => d["rel"], "G_honest" => Ghon,
        "nwords" => S.nwords, "nmoments" => S.nmoments, "grid" => GRIDS[d["grid"]], "host" => gethostname(),
        "started" => string(now()), "versions" => versions(), "points" => Any[], "status" => "running")
    clean(x::AbstractFloat) = isfinite(x) ? x : nothing
    clean(x::AbstractDict) = Dict(k => clean(y) for (k, y) in x)
    clean(x::AbstractVector) = [clean(y) for y in x]
    clean(x) = x
    save() = (open(outfile * ".tmp", "w") do io
        JSON.print(io, clean(res), 1)
    end; mv(outfile * ".tmp", outfile; force = true))
    @printf("[%s] %s built: %d words, %d moments, %.0f s\n", now(), id, S.nwords, S.nmoments, S.build_seconds)
    flush(stdout)

    function evaluate(dG; minonly = false)
        r = minonly ? minentropy_point(S; v = v, η = η, G = Ghon + dG) : keyrate_point(S, α, β; v = v, η = η, G = Ghon + dG)
        push!(res["points"], merge(Dict{String,Any}("dG" => dG, "minonly" => minonly), Dict(String(k) => x for (k, x) in pairs(r))))
        save()
        @printf("[%s] %s dG=%.5f H=%.5f Hmin=%.5f K=%+.5f Kmin=%+.5f\n", now(), id, dG, r.H, r.Hmin, r.K, r.Kmin)
        flush(stdout)
        return r
    end
    # Largest dG with a positive rate: scan upwards, then bisect. A failed program stops the search.
    function largest(rate; minonly)
        lo, hi = nothing, nothing
        dG = 0.0
        while dG <= dGmax + 1e-12
            k = rate(evaluate(dG; minonly = minonly))
            if !isfinite(k)                   # the programs can fail at the ideal point itself: skip it
                (lo === nothing && dG == 0.0) || return nothing
                dG = round(dG + step; digits = 8)
                continue
            end
            if k > Ktol
                lo = dG
            else
                hi = dG
                break
            end
            dG = round(dG + step; digits = 8)
        end
        (lo === nothing || hi === nothing) && return nothing
        while hi - lo > tol
            mid = (lo + hi) / 2
            k = rate(evaluate(mid; minonly = minonly))
            isfinite(k) || break
            k > Ktol ? (lo = mid) : (hi = mid)
        end
        return (lo + hi) / 2
    end
    res["dG_max"] = largest(r -> r.K; minonly = false)
    res["dG_max_minentropy"] = largest(r -> r.Kmin; minonly = true)
    res["status"] = "done"
    res["finished"] = string(now())
    save()
    @printf("[%s] %s DONE dG_max=%s dG_max_minentropy=%s\n", now(), id, string(res["dG_max"]), string(res["dG_max_minentropy"]))
end

main(ARGS)
