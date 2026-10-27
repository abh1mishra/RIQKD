# Shannon entropies and the coefficients of the Koßmann–Schwonnek bound.
# The expressions are those of utils.jl in A. Mishra's SDIQKD repository (https://github.com/abh1mishra/SDIQKD),
# written in the same order of floating-point operations, so that the values are bit-for-bit identical.

xlog2x(x) = (x < 0 || x > 1 || x == 0 || x == 1) ? 0.0 : -x * log2(x)

"Binary entropy h(x) in bits."
hbin(x) = xlog2x(x) + xlog2x(1 - x)

"""
    conditional_entropy(J)

Conditional Shannon entropy H(A|B) in bits of the joint distribution `J[a, b]`.
"""
function conditional_entropy(J::AbstractMatrix{Float64})
    nA, nB = size(J)
    H = 0.0
    PB = sum(J, dims = 1)
    for b in 1:nB
        PB[b] == 0.0 && continue
        for a in 1:nA
            J[a, b] == 0.0 && continue
            H -= J[a, b] * log2(J[a, b] / PB[b])
        end
    end
    return H
end

"""
    chord_coeffs(t) -> c

Weights c_k of the chord bound ∫_{t₁}^1 f(s) ds/s ≤ Σ_k c_k f(t_k) for a convex f on the nodes
0 < t₁ < … < t_r = 1 (paper, App. B).
"""
function chord_coeffs(t::AbstractVector)
    r = length(t)
    r >= 2 || error("the partition needs at least two nodes")
    (t[1] > 0 && t[end] == 1 && all(diff(t) .> 0)) || error("the nodes must satisfy 0 < t₁ < … < t_r = 1")
    c = zeros(r)
    c[1] = (1 + t[1] / (t[2] - t[1])) * log(t[2] / t[1]) - 1
    c[r] = 1 - (t[r-1] / (t[r] - t[r-1])) * log(t[r] / t[r-1])
    for k in 2:r-1
        c[k] = (1 + t[k] / (t[k+1] - t[k])) * log(t[k+1] / t[k]) - log(t[k] / t[k-1]) * t[k-1] / (t[k] - t[k-1])
    end
    return c
end
