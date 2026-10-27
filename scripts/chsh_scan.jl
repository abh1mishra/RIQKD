# Entropy as a function of the CHSH value (paper, Sec. VI.D and App. C): H(A|E) certified by S_CHSH ≥ S and
# P_g ≤ G, for the four BB84 preparations and Bob's two test measurements (Z±X)/√2 of the (4,3) protocol.
# The key is Alice's Z-basis bit.
#
#   julia --project=. scripts/chsh_scan.jl G=0.5 [rel=lean:4] [Smin=0] [Smax=4] [step=0.02] [Slist=2.5,2.7]
#                                          [grid=g26] [out=results/chsh]
#   julia --project=. scripts/chsh_scan.jl path=source [vmin=0.6] [vmax=0.99] [vstep=0.01] ...
#       (along the depolarised source: S = 2√2 v and G = (1+v)/4)
#
# For each S the script computes lower bounds on
#   F(S, G)  = min H(A|E)  subject to  (1/2) Σ_{x,y} C[x,y] E(x,y) ≥ S  and  P_g ≤ G,
# and upper bounds on Pg(S, G), Eve's probability of guessing the key bit under the same constraints. C are the
# signs of the ideal correlators (`chsh_coeffs`), so that S is the CHSH value of the entanglement-based twin.
# Only S is constrained: the marginals and the individual correlators are free. Next to each point the script
# writes the device-independent values H_DI(S) = 1 − h(1/2 + (1/2)√((S/2)² − 1)) and
# Pg_DI(S) = 1/2 + (1/2)√(2 − S²/4) (0 and 1 for S ≤ 2). They are meaningful only up to the quantum maximum
# S = 2√2; beyond it the formulas are clamped to 1 and 1/2. The scan stops after three consecutive infeasible values
# of S, beyond the largest S compatible with G.

include(joinpath(@__DIR__, "..", "src", "RIQKDNumerics.jl"))
using .RIQKDNumerics, JSON, Dates, Printf

function main(args)
    d = Dict{String,String}("rel" => "lean:4", "Smin" => "0.0", "Smax" => "4.0", "step" => "0.02", "grid" => "g26",
        "out" => "results/chsh", "Slist" => "")
    for a in args
        k, v = split(a, '='; limit = 2)
        d[k] = v
    end
    source = get(d, "path", "") == "source"
    G = source ? NaN : parse(Float64, d["G"])
    p = protocol("4chsh")
    C = chsh_coeffs(p)
    S = build_sdp(p, parse_relaxation(d["rel"]); functional = C)
    α, β = ks_coeffs(GRIDS[d["grid"]])
    Svals = isempty(d["Slist"]) ? collect(parse(Float64, d["Smin"]):parse(Float64, d["step"]):parse(Float64, d["Smax"])) :
            parse.(Float64, split(d["Slist"], ','))
    vvals = collect(parse(Float64, get(d, "vmax", "0.99")):-parse(Float64, get(d, "vstep", "0.01")):parse(Float64, get(d, "vmin", "0.6")))
    pairs = source ? [(2 * sqrt(2) * v, (1 + v) / 4, v) for v in vvals] : [(s, G, NaN) for s in Svals]
    hDI(s) = s <= 2 ? 0.0 : 1 - hbin(0.5 + 0.5 * sqrt((s / 2)^2 - 1))
    pgDI(s) = s <= 2 ? 1.0 : 0.5 + 0.5 * sqrt(max(2 - s^2 / 4, 0.0))

    mkpath(d["out"])
    id = replace("chsh__$(d["rel"])__" * (source ? "source" : "G$(d["G"])") * "__$(d["grid"])" * (isempty(d["Slist"]) ? "" : "__test"),
        r"[^A-Za-z0-9_.@+-]" => "_")
    outfile = joinpath(d["out"], id * ".json")
    res = Dict{String,Any}("job" => d, "G" => G, "relaxation" => d["rel"], "protocol" => p.name,
        "description" => p.description, "C" => [C[x, :] for x in 1:size(C, 1)], "theta_x" => p.θx, "theta_y" => p.θy,
        "nwords" => S.nwords, "nmoments" => S.nmoments, "build_seconds" => S.build_seconds, "grid" => GRIDS[d["grid"]],
        "host" => gethostname(), "started" => string(now()), "versions" => versions(), "points" => Any[],
        "status" => "running")
    clean(x::AbstractFloat) = isfinite(x) ? x : nothing
    clean(x::AbstractDict) = Dict(k => clean(v) for (k, v) in x)
    clean(x::AbstractVector) = [clean(v) for v in x]
    clean(x) = x
    save() = (open(outfile * ".tmp", "w") do io
        JSON.print(io, clean(res), 1)
    end; mv(outfile * ".tmp", outfile; force = true))
    @printf("[%s] %s built: %d words, %d moments, %.0f s\n", now(), id, S.nwords, S.nmoments, S.build_seconds)
    flush(stdout)

    infeasible_run = 0
    seen_feasible = false
    for (s, Gs, v) in pairs
        r = functional_point(S, α, β; value = s, G = Gs)
        push!(res["points"], Dict{String,Any}("S" => s, "G" => Gs, "v" => v, "H" => r.H, "Hmin" => r.Hmin, "Pg" => r.Pg,
            "H_DI" => hDI(s), "Pg_DI" => pgDI(s), "nonoptimal" => r.nonoptimal, "pg_status" => r.pg_status,
            "node_values" => r.node_values, "seconds" => r.seconds))
        save()
        @printf("[%s] %s S=%.3f G=%.4f H=%.5f (DI %.5f)  Pg=%.5f (DI %.5f)  %.0f s%s\n", now(), id, s, Gs, r.H, hDI(s), r.Pg, pgDI(s),
            r.seconds, isempty(r.nonoptimal) ? "" : " nonoptimal=$(length(r.nonoptimal))")
        flush(stdout)
        if isfinite(r.H)
            seen_feasible = true
            infeasible_run = 0
        else
            infeasible_run += 1
            (seen_feasible && infeasible_run >= 3) && break     # beyond the largest S compatible with G
        end
    end
    res["status"] = "done"
    res["finished"] = string(now())
    save()
end

main(ARGS)
