# Checks of the relaxation (paper, Sec. V.B and App. C), a few minutes in total:
#
#   julia --project=. scripts/validate.jl
#
#   C1  the honest guessing probabilities are 2/n (Sec. II.B);
#   C2  the lemma of Sec. V.B: the tracial relaxation, without the rows that lack a trace-class factor, has the
#       same value as the state-weighted relaxation on the same words in Bob's and Eve's projectors. Checked for
#       (4,2) at v = 0.98 with the partition t_k = k/10, for the tracial level 1+B*E (words 1, M, P, MP) and for
#       1+B*E+B*B+E*E (all words of length ≤ 2), as in App. C;
#   C3  with the rows that lack a trace-class factor, the dual program has no strictly feasible point and MOSEK
#       returns a slightly higher value (App. C);
#   C4  with the words of length ≤ 1 the relaxation is vacuous: Eve guesses the key bit with probability 1 even for
#       the ideal (4,2) behaviour (App. C);
#   C5  values of the state-weighted relaxation at three points, as reproduced by an independent implementation
#       (Python, CVXPY with MOSEK) to within 10⁻⁵ (App. C).

include(joinpath(@__DIR__, "..", "src", "RIQKDNumerics.jl"))
using .RIQKDNumerics, Printf

fails = 0
function check(name, got, want, tol)
    ok = abs(got - want) <= tol
    global fails += !ok
    @printf("%-66s %.7f  (expected %.7f ± %.0e)  %s\n", name, got, want, tol, ok ? "OK" : "** FAILED **")
    flush(stdout)
    return got
end

function H_at(prot, relstr, grid, v)
    p = protocol(prot)
    S = build_sdp(p, parse_relaxation(relstr))
    α, β = ks_coeffs(GRIDS[grid])
    r = keyrate_point(S, α, β; v = v, η = 1.0, G = honest_G(p), minentropy = false)
    @printf("   [%s %s, %s, v = %.2f: %d words, %d moments, build %.0f s, %d non-optimal programs]\n",
        prot, relstr, grid, v, S.nwords, S.nmoments, S.build_seconds, length(r.nonoptimal))
    return r.H
end

println("C1 honest guessing probabilities")
check("   (4,2)", honest_G(protocol("42")), 1 / 2, 1e-7)
check("   (6,3)", honest_G(protocol("6sym")), 1 / 3, 1e-7)
check("   (3,3)", honest_G(protocol("3st_ZBB")), 2 / 3, 1e-7)

println("C2 tracial relaxation = state-weighted relaxation on the same words ((4,2), v = 0.98, t_k = k/10)")
a = check("   tracial:1+B*E, rows without a trace-class factor removed", H_at("42", "tracial:1+B*E", "uniform10", 0.98), 0.379291, 2e-6)
b = check("   lean:1+B*E (words 1, M, P, MP)", H_at("42", "lean:1+B*E", "uniform10", 0.98), 0.379291, 2e-6)
check("   difference", b - a, 0.0, 2e-6)
a = check("   tracial:1+B*E+B*B+E*E, rows without a trace-class factor removed", H_at("42", "tracial:1+B*E+B*B+E*E", "uniform10", 0.98), 0.428383, 2e-6)
b = check("   lean:2 (all words of length ≤ 2)", H_at("42", "lean:2", "uniform10", 0.98), 0.428383, 2e-6)
check("   difference", b - a, 0.0, 2e-6)

println("C3 tracial relaxation with the rows that lack a trace-class factor (ill-posed dual)")
c = H_at("42", "tracial:1+B*E:keepfree", "uniform10", 0.98)
@printf("   tracial:1+B*E:keepfree gives %.7f, i.e. %.1e above the well-posed value\n", c, c - 0.379291)

println("C4 words of length ≤ 1: vacuous relaxation")
let p = protocol("42")
    S = build_sdp(p, parse_relaxation("lean:1"))
    set_behaviour!(S, behaviour(p), honest_G(p))
    Pg, st, _ = pguess(S)
    check("   Eve's guessing probability of the key bit, ideal (4,2)", Pg, 1.0, 1e-6)
end

println("C5 reference values (26-node partition)")
check("   (4,2) lean:2, v = 0.98", H_at("42", "lean:2", "g26", 0.98), 0.43168, 1e-5)
check("   (4,2) lean:3, v = 0.94", H_at("42", "lean:3", "g26", 0.94), 0.24956, 1e-5)
check("   (4,3) lean:3, v = 0.86", H_at("43", "lean:3", "g26", 0.86), 0.40546, 1e-5)

println(fails == 0 ? "ALL CHECKS PASSED" : "$fails CHECK(S) FAILED")
exit(fails == 0 ? 0 : 1)
