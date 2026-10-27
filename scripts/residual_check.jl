# Residual correction of the dual solutions (paper, Sec. V.C):
#
#   julia --project=. scripts/residual_check.jl                                  # the four points of Sec. V.C
#   julia --project=. scripts/residual_check.jl prot=43 rel=lean:3 v=0.86 [eta=1.0] [grid=g26]
#
# For each program, the dual matrices returned by MOSEK are projected onto the positive semidefinite cone, and the
# residual r of the dual equality constraints is recomputed. Every moment of the state-weighted relaxation has one
# trace-class factor and is bounded, |Tr(ρ_x w)| ≤ 1 and |Tr(σ w)| ≤ G, so subtracting max(1, G) ‖r‖₁ from the dual
# objective gives a bound that no longer depends on the solver tolerances. This is a floating-point check, not an
# interval-arithmetic certificate. The script prints the entropy bound before and after the correction.

include(joinpath(@__DIR__, "..", "src", "RIQKDNumerics.jl"))
using .RIQKDNumerics, JuMP, LinearAlgebra, Printf

"Dual objective after projecting the PSD blocks and subtracting the weighted equality residual."
function corrected_objective(S::SDPModel, G)
    m = S.model
    val = Dict(v => value(v) for v in all_variables(m))
    for (F, T) in list_of_constraint_types(m)
        T <: MOI.PositiveSemidefiniteConeTriangle || continue
        for c in all_constraints(m, F, T)
            obj = constraint_object(c)
            d = obj.set.side_dimension
            X = Matrix{VariableRef}(undef, d, d)
            k = 0
            for j in 1:d, i in 1:j
                k += 1
                X[i, j] = X[j, i] = obj.func[k]
            end
            e = eigen(Symmetric(value.(X)))
            A = e.vectors * Diagonal(max.(e.values, 0.0) .+ 1e-12) * e.vectors'
            for i in 1:d, j in 1:d
                val[X[i, j]] = A[i, j]
            end
        end
    end
    val[S.μ] = max(val[S.μ], 0.0)
    ev(v) = val[v]
    res = sum(abs(value(ev, constraint_object(c).func) - normalized_rhs(c)) for c in values(S.con))
    return value(ev, objective_function(m)) - max(1.0, G) * res
end

function check_point(prot, relstr; v = 1.0, η = 1.0, grid = "g26")
    p = protocol(prot)
    S = build_sdp(p, parse_relaxation(relstr))
    G = honest_G(p)
    set_behaviour!(S, behaviour(p; v = v, η = η), G)
    α, β = ks_coeffs(GRIDS[grid])
    raw = corr = 0.0
    for k in eachindex(α)
        z, ts, ps, _ = solve_node!(S, RIQKDNumerics.node_objective(p, α[k], β[k]))
        isfinite(z) || error("node $k failed: $ts / $ps")
        raw += z
        corr += corrected_objective(S, G)
    end
    H, Hc = (1 + raw) / log(2), (1 + corr) / log(2)
    @printf("%-8s %-12s v=%.3f eta=%.3f  H(A|E) = %.9f   after correction %.9f   (difference %.1e)\n",
        prot, relstr, v, η, H, Hc, H - Hc)
    flush(stdout)
    return H, Hc
end

if isempty(ARGS)
    # the known-zero control (three-state protocol of Tamaki et al., RI-classical) and three noisy points
    check_point("3st_ZX", "lean:3"; v = 1.0)
    check_point("42", "lean:2"; v = 0.98)
    check_point("42", "lean:3"; v = 0.94)
    check_point("43", "lean:3"; v = 0.86)
else
    d = Dict{String,String}("v" => "1.0", "eta" => "1.0", "grid" => "g26")
    for a in ARGS
        k, x = split(a, '='; limit = 2)
        d[k] = x
    end
    check_point(d["prot"], d["rel"]; v = parse(Float64, d["v"]), η = parse(Float64, d["eta"]), grid = d["grid"])
end
