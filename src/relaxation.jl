# Semidefinite relaxations of the key-rate problem (paper, Sec. V.B), built with PCPOP.
#
# Letters: B[y] = Bob's projector M_{0|y} (M_{1|y} = 𝟙 − B[y] is eliminated); E[2i-1], E[2i] = Eve's projectors
# P^i_0, P^i_1 for key basis i, which commute with every B[y]; BE[x] = ρ_x (x = 1..n) and BE[n+1] = σ, where
# σ ≥ ρ_x/n and Tr σ ≤ G encode P_g ≤ G (paper, Sec. II.B). Moments are L(w) = Tr(w). They are taken real, since
# the problem is invariant under complex conjugation, so that a word and its adjoint are identified.

"""
    Relaxation(kind, level, drop_free)

* `kind = :lean`: the state-weighted relaxation used for all the results (Sec. V.B). For each trace-class positive
  operator F in {ρ_x, σ − ρ_x/n}, the matrix [Tr(u† F v)]_{u,v ∈ W} is required to be positive semidefinite, where
  W is the PCPOP monomial set `level` in Bob's and Eve's projectors. There is no moment matrix, so every moment
  contains exactly one trace-class factor. (ρ_x ≥ 0 and Tr ρ_x = 1 already imply ρ_x ≤ 𝟙.) The level `"2s"` means
  all words of length at most two except the products of Eve's projectors belonging to different key bases;
  further PCPOP terms can be appended, e.g. `"2s+B*B*E"`.
* `kind = :tracial`: the tracial relaxation (Sec. V.B), used only by `scripts/validate.jl` to check the lemma of
  Sec. V.B. The moment matrix is indexed by {𝟙, ρ_x, σ}·S and the localising matrices of σ − ρ_x/n and ρ_x − ρ_x²
  by S, where S is the PCPOP monomial set `level` in all the letters, e.g. `"1+B*E"`. With `drop_free = true`, the
  rows of the moment matrix whose monomial contains neither a state nor σ are removed. Their moments
  (Tr 𝟙, Tr M, Tr MP, …) are unconstrained, so removing them does not change the value, but it makes the program
  well posed.
"""
struct Relaxation
    kind::Symbol
    level::String
    drop_free::Bool
end

"Parse `lean:2`, `lean:2s`, `lean:3`, `lean:2+B*B*E`, `tracial:1+B*E`, `tracial:1+B*E:keepfree`, ..."
function parse_relaxation(s::AbstractString)
    parts = split(s, ':')
    kind = Symbol(parts[1])
    kind in (:lean, :tracial) || error("unknown relaxation kind $(parts[1])")
    drop_free = !(length(parts) >= 3 && parts[3] == "keepfree")
    return Relaxation(kind, String(parts[2]), drop_free)
end

Base.string(r::Relaxation) = "$(r.kind):$(r.level)" * (r.kind == :tracial && !r.drop_free ? ":keepfree" : "")

"PCPOP level string of the state-weighted word set."
function lean_level_string(level::AbstractString, ν::Int)
    if startswith(level, "2s")        # "2s" or e.g. "2s+B*B*E"
        same = join(["E[$(2i-1):$(2i)]*E[$(2i-1):$(2i)]" for i in 1:ν], "+")
        return "1+B*B+B*E+" * same * level[3:end]
    end
    return String(level)      # e.g. "2", "3" (all words up to that length), or any PCPOP level string
end

"""
A built relaxation. The JuMP model is the dual of the moment relaxation, so its objective value at a feasible
point is a valid lower bound. The behaviour, G and the objective enter only through objective coefficients and
right-hand sides, so they can be changed without rebuilding the model.
"""
mutable struct SDPModel
    model::Model
    prot::Protocol
    rel::Relaxation
    λ::Matrix{VariableRef}         # multipliers of Tr(ρ_x M_{0|y}) = p(0|x,y)
    λn::Vector{VariableRef}        # multipliers of Tr ρ_x = 1
    μ::VariableRef                 # multiplier of Tr(−σ) ≥ −G
    con::Dict{Any,ConstraintRef}   # moment → dual equality constraint
    objkeys::Dict{NTuple{2,Int},Any}   # (Eve letter j, state x) → moment Tr(E_j ρ_x)
    nwords::Int
    blocksizes::Vector{Int}
    nmoments::Int
    build_seconds::Float64
    μS::Union{Nothing,VariableRef}   # multiplier of the functional constraint (keyword `functional` of build_sdp)
    Csum::Float64
end

"Add the affine expression `e` to the entry of moment `k`, identifying `k` with its adjoint."
function _addto!(total::Dict{Any,AffExpr}, k, e)
    if haskey(total, k)
        add_to_expression!(total[k], e)
        return
    end
    kc = conj(k)
    if haskey(total, kc)
        add_to_expression!(total[kc], e)
    else
        total[k] = copy(e)
    end
    return
end

"""
    build_sdp(p, rel; threads = 1, functional = nothing)

Build the relaxation `rel` for the protocol `p`, without behaviour, G or objective. These are set with
[`set_behaviour!`](@ref) and [`solve_node!`](@ref).

With `functional = C` (an n × m matrix), the behaviour is replaced by the single constraint
(1/2) Σ_{x,y} C[x,y] E(x,y) ≥ S on the correlators E(x,y) = 2 p(0|x,y) − 1, e.g. the CHSH value of Sec. VI.D;
S and G are then set with [`set_functional!`](@ref).
"""
function build_sdp(p::Protocol, rel::Relaxation; threads::Int = 1, functional::Union{Nothing,AbstractMatrix} = nothing)
    t0 = time()
    n, m, ν = nstates(p), nmeas(p), nbases(p)
    @pcmonoid Mon B[m, 0] E[2ν, 0] BE[n + 1, 0]
    @comms B E
    Projector.([B; E])
    build(Mon)
    ρ = BE[1:n]
    σ = BE[n+1]
    q = 1 / n

    blocks = Tuple{Any,Vector}[]
    nwords = 0
    if rel.kind == :lean
        W = mons_at_level(Vector{PCPOP.Variable}([B; E]), lean_level_string(rel.level, ν))
        nwords = length(W)
        for x in 1:n
            push!(blocks, (1 * ρ[x], W))
        end
        for x in 1:n
            push!(blocks, (σ - q * ρ[x], W))
        end
    else
        op_ge = vcat([σ - q * ρ[x] for x in 1:n], [ρ[x] - ρ[x] * ρ[x] for x in 1:n])
        S, P = get_monomials(0, rel.level; op_ge = op_ge, list_vars = Mon.vertices)
        if rel.drop_free
            states = BE
            P = filter(u -> any(v -> any(s -> s === v, states), variables(u)), P)
        end
        nwords = length(S)
        push!(blocks, (1, P))
        for g in op_ge
            push!(blocks, (g, S))
        end
    end

    model = Model(Mosek.Optimizer)
    set_silent(model)
    set_attribute(model, "MSK_IPAR_NUM_THREADS", threads)
    Z = [@variable(model, [1:length(r), 1:length(r)], PSD) for (_, r) in blocks]
    total = Dict{Any,AffExpr}()
    for (i, (g, r)) in enumerate(blocks)
        lmi = PCPOP.cyclic_npa_moment(r, Z[i]; cPoly = g)
        for (k, e) in lmi
            _addto!(total, k, e)
        end
    end
    mons = Set(keys(total))

    λn = @variable(model, [1:n])
    μ = @variable(model, lower_bound = 0)
    lin = Tuple{VariableRef,Any}[]
    μS = nothing
    Csum = 0.0
    if functional === nothing
        λ = @variable(model, [1:n, 1:m])
        for x in 1:n, y in 1:m
            push!(lin, (λ[x, y], 1 * (ρ[x] * B[y])))
        end
    else
        # (1/2) Σ C E ≥ S with E = 2 Tr(ρ_x M_{0|y}) − 1  ⟺  Σ C[x,y] Tr(ρ_x M_{0|y}) ≥ S + Σ C / 2   (multiplier μS ≥ 0)
        size(functional) == (n, m) || error("functional must be $(n)×$(m)")
        λ = Matrix{VariableRef}(undef, 0, 0)
        μS = @variable(model, lower_bound = 0)
        Csum = Float64(sum(functional))
        for x in 1:n, y in 1:m
            functional[x, y] == 0 && continue
            push!(lin, (μS, functional[x, y] * (ρ[x] * B[y])))
        end
    end
    for x in 1:n
        push!(lin, (λn[x], 1 * ρ[x]))
    end
    push!(lin, (μ, -1 * σ))
    for (var, poly) in lin
        red = PCPOP.trace_mons_reduce(mons, [poly])[1]
        for (k, c) in red
            add_to_expression!(total[k], c, var)
        end
    end

    objkeys = Dict{NTuple{2,Int},Any}()
    for (i, (x0, x1, _)) in enumerate(p.keypairs), j in (2i - 1, 2i), x in (x0, x1)
        red = PCPOP.trace_mons_reduce(mons, [1 * (E[j] * ρ[x])])[1]
        objkeys[(j, x)] = first(red.monomials)
    end

    con = Dict{Any,ConstraintRef}()
    for (k, e) in total
        con[k] = @constraint(model, e == 0.0)
    end
    @objective(model, Max, sum(λn) - μ)
    return SDPModel(model, p, rel, λ, λn, μ, con, objkeys, nwords,
        [length(r) for (_, r) in blocks], length(total), time() - t0, μS, Csum)
end

"Set the value S of the functional constraint (1/2) Σ C E ≥ S and the bound G (models built with `functional`)."
function set_functional!(S::SDPModel, value::Real, G::Real)
    S.μS === nothing && error("the model was not built with a functional constraint")
    set_objective_coefficient(S.model, S.μS, value + S.Csum / 2)
    set_objective_coefficient(S.model, S.μ, -G)
    return S
end

"Set the observed behaviour P0[x, y] = p(0|x,y) and the bound G."
function set_behaviour!(S::SDPModel, P0::AbstractMatrix, G::Real)
    n, m = size(S.λ)
    for x in 1:n, y in 1:m
        set_objective_coefficient(S.model, S.λ[x, y], P0[x, y])
    end
    set_objective_coefficient(S.model, S.μ, -G)
    return S
end

"""
    solve_node!(S, f) -> (value, termination_status, primal_status, seconds)

Minimise Σ f[(j, x)] Tr(E_j ρ_x) over the relaxation, with the current behaviour and G. The value returned is the
objective of the dual program at MOSEK's solution. It is `NaN` unless MOSEK reports a feasible or nearly feasible
point.
"""
function solve_node!(S::SDPModel, f::Dict{NTuple{2,Int},Float64})
    for (jx, k) in S.objkeys
        set_normalized_rhs(S.con[k], get(f, jx, 0.0))
    end
    t = @elapsed optimize!(S.model)
    ts, ps = termination_status(S.model), primal_status(S.model)
    ok = ps == MOI.FEASIBLE_POINT || ps == MOI.NEARLY_FEASIBLE_POINT
    return (ok ? objective_value(S.model) : NaN), string(ts), string(ps), t
end
