# Entropy bounds and key rates.
#
# Koßmann–Schwonnek bound (paper, Sec. V.A and App. A):
#     H(A|E) ≥ (1/ln 2) [ 1 + Σ_k min Σ_i Tr Θ^i_k ],
#     Θ^i_k = −(1/2ν) [ α_k (P^i_0 ρ_{x(i,0)} + P^i_1 ρ_{x(i,1)}) + β_k (P^i_0 + P^i_1)(ρ_{x(i,0)} + ρ_{x(i,1)}) ],
# averaged over the ν key bases i. Each node k is minimised separately with the same relaxation (decoupling over
# the nodes, Sec. V.C), which is a further relaxation and gives a valid lower bound.

"""
Partitions 0 < t₁ < … < t_r = 1 of the integration interval. The k = 0 term is always added by `ks_coeffs`.

* `g26`: the 26-node partition of the paper (Sec. V.C), used for all the results except the comparisons of partitions;
* `uniform10`: the uniform partition t_k = k/10 (App. C, Table III);
* `fine120`: 120 nodes, geometric below 0.05 and uniform above, for convergence checks.
"""
const GRIDS = Dict(
    "g26" => [0.005, 0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0],
    "uniform10" => collect(0.1:0.1:1.0),
    "fine120" => sort(unique(round.([exp.(range(log(1e-4), log(0.05), length = 25)); collect(range(0.05, 1.0, length = 96))];
        digits = 10))),
)

"""
    ks_coeffs(t; merge_k0 = true) -> (α, β)

Coefficients α_k, β_k of the entropy bound for the nodes `t` (paper, App. B): α_k = −c_k and β_k = c_k t_k for
k = 1, …, r, and α₀ = −1, β₀ = t₁ for the term that bounds the integral over [0, t₁]. That term has the same
operator as the term k = 1, so by default the two are merged, which saves one program per point.
"""
function ks_coeffs(t::AbstractVector; merge_k0::Bool = true)
    c = chord_coeffs(t)
    α, β = -c, t .* c
    merge_k0 || return [-1.0; α], [t[1]; β]
    return [-1.0 + α[1]; α[2:end]], [t[1] + β[1]; β[2:end]]
end

"Objective coefficients f[(j, x)] of Tr(E_j ρ_x) for one node (α, β)."
function node_objective(p::Protocol, α::Real, β::Real)
    ν = nbases(p)
    f = Dict{NTuple{2,Int},Float64}()
    add!(k, c) = (f[k] = get(f, k, 0.0) + c)
    for (i, (x0, x1, _)) in enumerate(p.keypairs)
        P0, P1 = 2i - 1, 2i
        add!((P0, x0), -0.5 * (α + β) / ν)
        add!((P1, x1), -0.5 * (α + β) / ν)
        add!((P0, x1), -0.5 * β / ν)
        add!((P1, x0), -0.5 * β / ν)
    end
    return f
end

"""
    entropy_bound(S, α, β) -> (H, node_values, statuses, seconds)

Lower bound on H(A|E) in bits for the behaviour and G currently set in `S`. It is `NaN` if any node failed.
"""
function entropy_bound(S::SDPModel, α::AbstractVector, β::AbstractVector)
    vals = Float64[]
    stats = String[]
    secs = 0.0
    for k in eachindex(α)
        v, ts, ps, t = solve_node!(S, node_objective(S.prot, α[k], β[k]))
        push!(vals, v)
        push!(stats, ts == "OPTIMAL" ? "OPTIMAL" : "$ts/$ps")
        secs += t
    end
    return (1 + sum(vals)) / log(2), vals, stats, secs
end

"""
    pguess(S) -> (Pg, status, seconds)

Upper bound on Eve's probability of guessing the key bit, averaged over the key bases, from the same relaxation.
Eve measures {E_{2i-1}, 𝟙 − E_{2i-1}} in basis i. The min-entropy bound is H_min(A|E) = −log₂ Pg.
"""
function pguess(S::SDPModel)
    ν = nbases(S.prot)
    f = Dict{NTuple{2,Int},Float64}()
    for (i, (x0, x1, _)) in enumerate(S.prot.keypairs)
        f[(2i - 1, x0)] = -1 / (2ν)
        f[(2i - 1, x1)] = 1 / (2ν)
    end
    v, ts, ps, t = solve_node!(S, f)
    return 0.5 - v, (ts == "OPTIMAL" ? ts : "$ts/$ps"), t
end

"Behaviour matrix P0[x, y] = p(0|x,y) for visibility v and efficiency η (no-click binned into b = 0)."
behaviour(p::Protocol; v = 1.0, η = 1.0) = [prob0(p, x, y; v = v, η = η) for x in 1:nstates(p), y in 1:nmeas(p)]

"""
    keyrate_point(S, α, β; v, η, G, minentropy = true) -> NamedTuple

Entropy bounds and key rates at one point. `K = H − H(A|B)` is the rate per key round and `R = p_key K` the
rate per signal for uniform inputs; `Kmin` and `Rmin` are the same with the min-entropy.
"""
function keyrate_point(S::SDPModel, α, β; v = 1.0, η = 1.0, G, minentropy::Bool = true)
    p = S.prot
    set_behaviour!(S, behaviour(p; v = v, η = η), G)
    H, vals, stats, secs = entropy_bound(S, α, β)
    Pg, pstat, psecs = minentropy ? pguess(S) : (NaN, "skipped", 0.0)
    hab = HAB(p; v = v, η = η)
    K = H - hab
    Hmin = -log2(Pg)
    return (v = v, η = η, G = G, H = H, Hmin = Hmin, Pg = Pg, HAB = hab, K = K, Kmin = Hmin - hab,
        R = pkey(p) * K, Rmin = pkey(p) * (Hmin - hab),
        nonoptimal = [s for s in stats if s != "OPTIMAL"], pg_status = pstat,
        node_values = vals, seconds = secs + psecs)
end

"Min-entropy only (one program per point), used to bisect the min-entropy threshold."
function minentropy_point(S::SDPModel; v = 1.0, η = 1.0, G)
    p = S.prot
    set_behaviour!(S, behaviour(p; v = v, η = η), G)
    Pg, pstat, secs = pguess(S)
    hab = HAB(p; v = v, η = η)
    Hmin = -log2(Pg)
    return (v = v, η = η, G = G, H = NaN, Hmin = Hmin, Pg = Pg, HAB = hab, K = NaN, Kmin = Hmin - hab,
        R = NaN, Rmin = pkey(p) * (Hmin - hab), nonoptimal = String[], pg_status = pstat,
        node_values = Float64[], seconds = secs)
end

"""
    functional_point(S, α, β; value, G) -> NamedTuple

Entropy bounds when the behaviour is replaced by the single constraint (1/2) Σ C E ≥ `value` (models built with
the keyword `functional`, see `build_sdp`), with the RI bound G.
"""
function functional_point(S::SDPModel, α, β; value::Real, G::Real)
    set_functional!(S, value, G)
    H, vals, stats, secs = entropy_bound(S, α, β)
    Pg, pstat, psecs = pguess(S)
    return (S = value, G = G, H = H, Hmin = -log2(Pg), Pg = Pg,
        nonoptimal = [s for s in stats if s != "OPTIMAL"], pg_status = pstat, node_values = vals, seconds = secs + psecs)
end
