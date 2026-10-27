# Protocols, honest statistics, the error-correction term H(A|B) and the honest value of G.
#
# Conventions:
#   * an angle θ denotes the qubit state |θ⟩ = cos θ |0⟩ + sin θ |1⟩, i.e. Bloch angle 2θ in the XZ plane;
#   * the honest state x is  v |θ_x⟩⟨θ_x| + (1 − v) 𝟙/2, with visibility v (paper, Sec. VI.A);
#   * Bob's measurement y is {|θ_y⟩⟨θ_y|, 𝟙 − |θ_y⟩⟨θ_y|}, with detection efficiency η. In the semidefinite program
#     the no-click event is binned into the outcome b = 0:
#         M_{0|y} = η |θ_y⟩⟨θ_y| + (1 − η) 𝟙,   M_{1|y} = η (𝟙 − |θ_y⟩⟨θ_y|);
#     in H(A|B) the no-click outcome is kept separate.
# For the four-state protocols the states are numbered 1 → |0⟩, 2 → |1⟩, 3 → |+⟩, 4 → |−⟩. This differs from the
# labelling of the paper (Sec. IV.A), which is 1 → |0⟩, 2 → |+⟩, 3 → |1⟩, 4 → |−⟩. For `6sym` the numbering is that
# of the paper (Sec. IV.C).

"""
    Protocol(name, θx, θy, keypairs, description)

`θx[x]` and `θy[y]` are the preparation and measurement angles. `keypairs[i] = (x0, x1, y)` means that key basis
`i` encodes the key bit a = 0 in state `x0` and a = 1 in state `x1`, and that Bob measures `y` in the key rounds.
"""
struct Protocol
    name::String
    θx::Vector{Float64}
    θy::Vector{Float64}
    keypairs::Vector{NTuple{3,Int}}
    description::String
end

nstates(p::Protocol) = length(p.θx)
nmeas(p::Protocol) = length(p.θy)
nbases(p::Protocol) = length(p.keypairs)
"Fraction of rounds that generate key, for uniform inputs x and y."
pkey(p::Protocol) = 2 * nbases(p) / (nstates(p) * nmeas(p))

const FOUR = [0, π / 2, π / 4, 3π / 4]          # |0⟩, |1⟩, |+⟩, |−⟩
const THREE = [0, π / 2, π / 4]                 # |0⟩, |1⟩, |+⟩

const PROTOCOLS = Dict{String,Protocol}(
    # protocols of the main text (Secs. IV and VI)
    "42" => Protocol("42", FOUR, [0, π / 8], [(1, 2, 1)],
        "(4,2): BB84 states; Z and (Z+X)/√2; key from Z"),
    "43" => Protocol("43", FOUR, [0, π / 8, -π / 8], [(1, 2, 1)],
        "(4,3): BB84 states; Z, (Z±X)/√2; key from Z"),
    "44" => Protocol("44", FOUR, [0, π / 8, -π / 8, π / 4], [(1, 2, 1), (3, 4, 4)],
        "(4,4): BB84 states; Z, (Z±X)/√2, X; key from Z and X"),
    "6sym" => Protocol("6sym", [0, π / 6, π / 3, π / 2, 2π / 3, 5π / 6], [0, π / 6, π / 3],
        [(1, 4, 1), (2, 5, 2), (3, 6, 3)],
        "(6,3): three coplanar bases at Bloch angles 0°, 60°, 120°, each measured and used for key"),
    "3st_ZBB" => Protocol("3st_ZBB", THREE, [0, π / 8, -π / 8], [(1, 2, 1)],
        "(3,3): states |0⟩, |1⟩, |+⟩; Z, (Z±X)/√2; key from Z"),
    # variants of App. C (Table IV)
    "63" => Protocol("63", [0, π / 4, π / 8, π / 2, 3π / 4, 5π / 8], [0, π / 4, π / 8],
        [(1, 4, 1), (2, 5, 2), (3, 6, 3)],
        "(6,3) with the bases of Z, X, (Z+X)/√2, each measured and used for key"),
    "BB84+1" => Protocol("BB84+1", FOUR, [0, π / 4, π / 8], [(1, 2, 1), (3, 4, 2)],
        "BB84 states; Z, X, (Z+X)/√2; key from Z and X"),
    "3st_ZB" => Protocol("3st_ZB", THREE, [0, π / 8], [(1, 2, 1)],
        "states |0⟩, |1⟩, |+⟩; Z, (Z+X)/√2; key from Z"),
    "3st_ZXB" => Protocol("3st_ZXB", THREE, [0, π / 4, π / 8], [(1, 2, 1)],
        "states |0⟩, |1⟩, |+⟩; Z, X, (Z+X)/√2; key from Z"),
    "3st_ZX" => Protocol("3st_ZX", THREE, [0, π / 4], [(1, 2, 1)],
        "three-state protocol: |0⟩, |1⟩, |+⟩; Z, X; key from Z (RI-classical by the proposition of Sec. III.C: no key)"),
    # the CHSH scenario of Sec. VI.D: the four BB84 states and the two test measurements (Z±X)/√2 of (4,3)
    "4chsh" => Protocol("4chsh", FOUR, [π / 8, -π / 8], [(1, 2, 1)],
        "BB84 states; (Z±X)/√2; key = Z-basis bit (Bob's key measurement is not part of the constraints)"),
)

"""
    protocol(name)

Look up a protocol. `"42@φ"` is the variant of (4,2) whose test measurement is cos φ Z + sin φ X instead of
(Z+X)/√2, and `"43@φ"` the variant of (4,3) whose test measurements are cos φ Z ± sin φ X, with φ in degrees
(App. C).
"""
function protocol(name::AbstractString)
    if occursin('@', name)
        base, φ = split(name, '@')
        return tilted(PROTOCOLS[base], parse(Float64, φ))
    end
    return PROTOCOLS[name]
end

function tilted(p::Protocol, φdeg::Real)
    θ = deg2rad(φdeg) / 2
    θy = p.name == "42" ? [0, θ] : p.name == "43" ? [0, θ, -θ] : error("tilts are defined for 42 and 43 only")
    return Protocol("$(p.name)@$(φdeg)", p.θx, θy, p.keypairs, p.description * "; tilt $(φdeg)°")
end

ketθ(θ) = [cos(θ), sin(θ)]
projθ(θ) = ketθ(θ) * ketθ(θ)'
honest_state(p::Protocol, x; v = 1.0) = v * projθ(p.θx[x]) + (1 - v) * Matrix(I, 2, 2) / 2

"Probability of Bob's outcome 0, with the no-click event binned into it: the data of the semidefinite program."
function prob0(p::Protocol, x, y; v = 1.0, η = 1.0)
    M0 = η * projθ(p.θy[y]) + (1 - η) * Matrix(I, 2, 2)
    return tr(honest_state(p, x; v = v) * M0)
end

"H(A|B) of the key rounds in bits, with the no-click event kept as a third outcome, averaged over the key bases."
function HAB(p::Protocol; v = 1.0, η = 1.0)
    tot = 0.0
    for (x0, x1, y) in p.keypairs
        J = zeros(2, 3)
        for (a, x) in enumerate((x0, x1))
            ρ = honest_state(p, x; v = v)
            J[a, 1] = 0.5 * η * tr(ρ * projθ(p.θy[y]))
            J[a, 2] = 0.5 * η * tr(ρ * (I - projθ(p.θy[y])))
            J[a, 3] = 0.5 * (1 - η)
        end
        tot += conditional_entropy(abs.(J))
    end
    return tot / nbases(p)
end

"""
    honest_G(p; v = 1)

Guessing probability of x from the honest qubit states with the uniform prior,
max over POVMs {N_x} of (1/n) Σ_x Tr(ρ_x N_x), computed by semidefinite programming. This is the bound G implied by
the target states (v = 1), or by the depolarised states when the noise is attributed to the source. For all
the protocols above it equals (1 + v)/n up to the solver accuracy, about 10⁻⁸ (Sec. V.C).
"""
function honest_G(p::Protocol; v = 1.0)
    n = nstates(p)
    m = Model(Mosek.Optimizer)
    set_silent(m)
    N = [@variable(m, [1:2, 1:2], PSD) for _ in 1:n]
    @constraint(m, sum(N) .== Matrix(I, 2, 2))
    @objective(m, Max, sum(tr(honest_state(p, x; v = v) * N[x]) for x in 1:n) / n)
    optimize!(m)
    return objective_value(m)
end

"""
    chsh_coeffs(p) -> C

Signs of the honest correlators, C[x,y] = sign E(x,y), so that S = (1/2) Σ C[x,y] E(x,y) is the CHSH-type
expression maximised by the target states and measurements. For `4chsh`, S is the CHSH value of the
entanglement-based twin (paper, Sec. VI.D), and S = 2√2 for the ideal states.
"""
function chsh_coeffs(p::Protocol)
    E = [cos(2 * (p.θx[x] - p.θy[y])) for x in 1:nstates(p), y in 1:nmeas(p)]
    return [abs(e) < 1e-12 ? 0.0 : sign(e) for e in E]
end
