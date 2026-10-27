"""
    RIQKDNumerics

Lower bounds on the asymptotic key rate of prepare-and-measure QKD under the restricted-information (RI)
assumption: the probability of guessing Alice's input from everything that leaves her laboratory is at most G.

This is the code of the numerical results of *Prepare-and-measure quantum key distribution from a bound on
information leakage* (A. Mishra, M. Masini, M. Balanzó-Juandó, A. Tavakoli, S. Pironio). The conditional entropy
H(A|E) is bounded with the Koßmann–Schwonnek relative-entropy bound on a fixed partition of [0, 1] (Sec. V.A and
App. B of the paper). Each node of that bound is a non-commutative polynomial program, relaxed by the
state-weighted moment relaxation of Sec. V.B. The relaxations are generated with PCPOP.jl
(https://github.com/abh1mishra/PCPOP.jl) and solved in dual form with MOSEK through JuMP.

Files:
* `utils.jl`: Shannon entropies and the coefficients of the entropy bound;
* `protocols.jl`: the protocols, their honest statistics, H(A|B) and the honest value of G;
* `relaxation.jl`: construction of the semidefinite programs;
* `entropy.jl`: entropy partitions, the entropy bound, the min-entropy and the key rate at one point.
"""
module RIQKDNumerics

using LinearAlgebra, Printf, Dates
using JuMP, Mosek, MosekTools
using PCPOP
import Pkg

include("utils.jl")
include("protocols.jl")
include("relaxation.jl")
include("entropy.jl")

export Protocol, PROTOCOLS, protocol, tilted, nstates, nmeas, nbases, pkey, honest_G, prob0, HAB, chsh_coeffs
export Relaxation, parse_relaxation, SDPModel, build_sdp, set_behaviour!, set_functional!, solve_node!
export GRIDS, ks_coeffs, entropy_bound, pguess, keyrate_point, minentropy_point, functional_point, behaviour
export hbin, versions

"Versions of Julia, MOSEK and the main packages, recorded in every result file."
function versions()
    d = Dict{String,Any}("julia" => string(VERSION))
    d["MOSEK"] = try
        join(Mosek.getversion(), ".")
    catch
        "unknown"
    end
    for (_, info) in Pkg.dependencies()
        info.name in ("PCPOP", "JuMP", "MathOptInterface", "Mosek", "MosekTools") || continue
        d[info.name] = info.is_tracking_repo ? "$(info.git_source)#$(info.git_revision)" : string(info.version)
    end
    return d
end

end # module
