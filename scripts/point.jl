# Entropy bounds and key rate at a single point, e.g. to compare relaxations or partitions (App. C):
#
#   julia --project=. scripts/point.jl prot=43 rel=lean:3 v=0.86 [eta=1.0] [dG=0.00] [grid=g26] [minentropy=true]
#
#   dG=source sets G to the guessing probability of the depolarised source.
# Prints the size of the program, the time per program, H(A|E), H_min(A|E), H(A|B) and the rates per key round.

include(joinpath(@__DIR__, "..", "src", "RIQKDNumerics.jl"))
using .RIQKDNumerics, Printf

d = Dict{String,String}("v" => "1.0", "eta" => "1.0", "dG" => "0.00", "grid" => "g26", "minentropy" => "true")
for a in ARGS
    k, val = split(a, '='; limit = 2)
    d[k] = val
end
p = protocol(d["prot"])
v, η = parse(Float64, d["v"]), parse(Float64, d["eta"])
G = d["dG"] == "source" ? honest_G(p; v = v) : honest_G(p) + parse(Float64, d["dG"])
S = build_sdp(p, parse_relaxation(d["rel"]))
α, β = ks_coeffs(GRIDS[d["grid"]])
nprog = length(α) + (d["minentropy"] == "true")
t = @elapsed r = keyrate_point(S, α, β; v = v, η = η, G = G, minentropy = d["minentropy"] == "true")
@printf("%s %s v=%.4f eta=%.4f G=%.10f grid=%s(%d nodes): %d words, %d moments, build %.0f s, %.2f s per program\n",
    p.name, d["rel"], v, η, G, d["grid"], length(GRIDS[d["grid"]]), S.nwords, S.nmoments, S.build_seconds, t / nprog)
@printf("H(A|E) = %.7f   Hmin(A|E) = %.7f   H(A|B) = %.7f   K = %+.7f   Kmin = %+.7f   non-optimal programs: %d\n",
    r.H, r.Hmin, r.HAB, r.K, r.Kmin, length(r.nonoptimal))
