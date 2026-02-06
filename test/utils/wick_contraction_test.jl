using Test
using SymPy

@testset "Wick Contraction Utils" begin
    @test creator(orbital_class("p")) ==
        Operator(:create, Sym("adag"), Sym("p"), Sym("a^{†}_{p}"), nothing)
    @test annihilator(orbital_class("q")) ==
        Operator(:annihilate, Sym("a"), Sym("q"), Sym("a_{q}"), nothing)
    @test creator(orbital_class("p_{α}")) ==
        Operator(:create, Sym("adag"), Sym("p_{α}"), Sym("a^{†}_{p_{α}}"), nothing)
    @test annihilator(orbital_class("q_{β}")) ==
        Operator(:annihilate, Sym("a"), Sym("q_{β}"), Sym("a_{q_{β}}"), nothing)
    @test one_body_operator(orbital_class("p_{α}"), orbital_class("q_{β}")) ==
        [Operator(:create, Sym("adag"), Sym("p_{α}"), Sym("a^{†}_{p_{α}}"), nothing),
         Operator(:annihilate, Sym("a"), Sym("q_{β}"), Sym("a_{q_{β}}"), nothing)]
    @test two_body_operator(orbital_class("p_{α}"), orbital_class("q_{β}"),
        orbital_class("r_{α}"), orbital_class("s_{β}")) ==
        [Operator(:create, Sym("adag"), Sym("p_{α}"), Sym("a^{†}_{p_{α}}"), nothing),
         Operator(:create, Sym("adag"), Sym("r_{α}"), Sym("a^{†}_{r_{α}}"), nothing),
         Operator(:annihilate, Sym("a"), Sym("s_{β}"), Sym("a_{s_{β}}"), nothing),
         Operator(:annihilate, Sym("a"), Sym("q_{β}"), Sym("a_{q_{β}}"), nothing)]

    @test two_body_element(orbital_class("p"), orbital_class("q"),
        orbital_class("r"), orbital_class("s")).name ==
        Sym("j_{pqrs}") - Sym("k_{pqrs}")
    @test two_body_element(orbital_class("p"), orbital_class("q"),
        orbital_class("r"), orbital_class("s"); swap=true).name ==
        Sym("j_{rspq}") - Sym("k_{rspq}")
    @test one_body_element(orbital_class("p"), orbital_class("q"); symbol="f").name ==
        Sym("f_{pq}")
    @test two_body_element(orbital_class("p"), orbital_class("q"),
        orbital_class("r"), orbital_class("s"); symbol_j="v", symbol_k="w").name ==
        Sym("v_{pqrs}") - Sym("w_{pqrs}")
    @test one_body_term(orbital_class("p"), orbital_class("q"); symbol="f") ==
        Sym("f_{pq}") * Sym("a^{†}_{p}") * Sym("a_{q}")
    @test two_body_term(orbital_class("p"), orbital_class("q"),
        orbital_class("r"), orbital_class("s"); symbol_j="v", symbol_k="w") ==
        Sym(1//4) * (Sym("v_{pqrs}") - Sym("w_{pqrs}")) *
        Sym("a^{†}_{p}") * Sym("a^{†}_{r}") * Sym("a_{s}") * Sym("a_{q}")

    ham = hamiltonian(orbital_class("p"), orbital_class("q"),
        orbital_class("r"), orbital_class("s"); c_ex=0.5)
    @test ham.h_pq.name == Sym("h_{pq}")
    @test ham.g_pqrs.name ==
        Sym("j_{pqrs}") - Sym(0.5) * Sym("k_{pqrs}")
    @test ham.c_ex == Sym(0.5)
    @test ham.one_body ==
        [Operator(:create, Sym("adag"), Sym("p"), Sym("a^{†}_{p}"), nothing),
         Operator(:annihilate, Sym("a"), Sym("q"), Sym("a_{q}"), nothing)]
    @test ham.two_body ==
        [Operator(:create, Sym("adag"), Sym("p"), Sym("a^{†}_{p}"), nothing),
         Operator(:create, Sym("adag"), Sym("r"), Sym("a^{†}_{r}"), nothing),
         Operator(:annihilate, Sym("a"), Sym("s"), Sym("a_{s}"), nothing),
         Operator(:annihilate, Sym("a"), Sym("q"), Sym("a_{q}"), nothing)]

    state = single_excitation_state(orbital_class("i_{α}"), orbital_class("a_{β}"))
    @test state.x_ai == Sym("x^a_{β}_i_{α}")
    @test state.y_ai == Sym("y^a_{β}_i_{α}")
    @test state.t_ai == Sym("x^a_{β}_i_{α}")
    @test state.op1 ==
        [Operator(:create, Sym("adag"), Sym("a_{β}"), Sym("a^{†}_{a_{β}}"), 0),
         Operator(:annihilate, Sym("a"), Sym("i_{α}"), Sym("a_{i_{α}}"), 1)]

    dstate = double_excitation_state(orbital_class("i"), orbital_class("j"),
        orbital_class("a"), orbital_class("b"))
    @test dstate.t_ai == Sym("t^a_i")
    @test dstate.t_aibj == Sym("t^ab_ij")
    @test dstate.op2 ==
        [Operator(:create, Sym("adag"), Sym("a"), Sym("a^{†}_{a}"), 0),
         Operator(:create, Sym("adag"), Sym("b"), Sym("a^{†}_{b}"), 0),
         Operator(:annihilate, Sym("a"), Sym("j"), Sym("a_{j}"), 1),
         Operator(:annihilate, Sym("a"), Sym("i"), Sym("a_{i}"), 1)]
    dstate_spin = double_excitation_state(orbital_class("i_{α}"), orbital_class("j_{β}"),
        orbital_class("a_{α}"), orbital_class("b_{β}"))
    @test dstate_spin.t_ai == Sym("t^a_{α}_i_{α}")
    @test dstate_spin.t_aibj == Sym("t^a_{α}b_{β}_i_{α}j_{β}")
    @test dstate_spin.op2 ==
        [Operator(:create, Sym("adag"), Sym("a_{α}"), Sym("a^{†}_{a_{α}}"), 0),
         Operator(:create, Sym("adag"), Sym("b_{β}"), Sym("a^{†}_{b_{β}}"), 0),
         Operator(:annihilate, Sym("a"), Sym("j_{β}"), Sym("a_{j_{β}}"), 1),
         Operator(:annihilate, Sym("a"), Sym("i_{α}"), Sym("a_{i_{α}}"), 1)]

    op_delta1 = creator(orbital_class("i")) * annihilator(orbital_class("j"))
    @test delta(op_delta1[1], op_delta1[2]).name == Sym("δ_{ij}")
    op_delta2 = creator(orbital_class("a")) * annihilator(orbital_class("b"))
    @test delta(op_delta2[1], op_delta2[2]).name == Sym("δ_{ab}")
    op_delta3 = creator(orbital_class("s")) * annihilator(orbital_class("t"))
    @test delta(op_delta3[1], op_delta3[2]).name == Sym("δ_{st}")
    @test orbital_class("s").type == :gen
    op_delta4 = creator(orbital_class("s")) * annihilator(orbital_class("i"))
    @test delta(op_delta4[1], op_delta4[2]).name == Sym("δ_{si}")
    op_delta5 = creator(orbital_class("s")) * annihilator(orbital_class("a"))
    @test delta(op_delta5[1], op_delta5[2]).name == Sym("δ_{sa}")
    op_delta6 = creator(orbital_class("i")) * annihilator(orbital_class("a"))
    @test delta(op_delta6[1], op_delta6[2]).name == Sym("δ_{ia}")
    @test delta(Sym("α"), Sym("β")).name == Sym(0)
    @test delta(Sym("α"), Sym("α")).name == Sym(1)

    op1 = annihilator(orbital_class("i_{α}"))
    op2 = creator(orbital_class("j_{α}"))
    @test contraction(op1, op2).name == Sym(0)

    op1b = creator(orbital_class("i_{α}"))
    op2b = annihilator(orbital_class("i_{α}"))
    @test contraction(op1b, op2b).name ==
        Sym("δ_{i_{α}i_{α}}")

    op3 = annihilator(orbital_class("s_{β}"))
    op4 = creator(orbital_class("t_{β}"))
    @test contraction(op3, op4).name == Sym("δ_{s_{β}t_{β}}")

    spin = spin_state_excitation(orbital_class("i"), orbital_class("a"))
    @test spin.s_zero isa OpSum
    @test spin.t_plus isa OpSum
    @test spin.t_minus isa OpSum
    @test spin.t_zero isa OpSum
    @test length(spin.s_zero.terms) == 2
    @test length(spin.t_plus.terms) == 1
    @test length(spin.t_minus.terms) == 1
    @test length(spin.t_zero.terms) == 2

    op_pair = creator(orbital_class("i")) * annihilator(orbital_class("a"))
    @test op_pair ==
        [Operator(:create, Sym("adag"), Sym("i"), Sym("a^{†}_{i}"), 1),
         Operator(:annihilate, Sym("a"), Sym("a"), Sym("a_{a}"), 0)]

    op_chain = op_pair * creator(orbital_class("b"))
    @test op_chain ==
        [Operator(:create, Sym("adag"), Sym("i"), Sym("a^{†}_{i}"), 1),
         Operator(:annihilate, Sym("a"), Sym("a"), Sym("a_{a}"), 0),
         Operator(:create, Sym("adag"), Sym("b"), Sym("a^{†}_{b}"), 0)]

    op_chain2 = annihilator(orbital_class("j")) * op_pair
    @test op_chain2 ==
        [Operator(:annihilate, Sym("a"), Sym("j"), Sym("a_{j}"), 1),
         Operator(:create, Sym("adag"), Sym("i"), Sym("a^{†}_{i}"), 1),
         Operator(:annihilate, Sym("a"), Sym("a"), Sym("a_{a}"), 0)]

    mb_two_flat = many_body_operator([orbital_class("p"), orbital_class("q")])
    mb_two = normal_order_op_pairs(mb_two_flat)
    @test mb_two_flat == one_body_operator(orbital_class("p"), orbital_class("q"))
    @test mb_two == [one_body_operator(orbital_class("p"), orbital_class("q"))]

    mb_four_flat = many_body_operator([
        orbital_class("p"), orbital_class("q"),
        orbital_class("r"), orbital_class("s")
    ])
    mb_four = normal_order_op_pairs(mb_four_flat)
    @test mb_four_flat == two_body_operator(
        orbital_class("p"), orbital_class("q"),
        orbital_class("r"), orbital_class("s")
    )
    @test mb_four == [
        one_body_operator(orbital_class("p"), orbital_class("s")),
        one_body_operator(orbital_class("p"), orbital_class("q")),
        one_body_operator(orbital_class("r"), orbital_class("s")),
        one_body_operator(orbital_class("r"), orbital_class("q"))
    ]

    h_pq = one_body_element(orbital_class("p"), orbital_class("q"); symbol="h")
    δap = delta(Sym("a"), Sym("p"))
    δbq = delta(Sym("b"), Sym("q"))
    expr1 = OpSum(IndicedObj[δap, δbq])
    @test SymPy.simplify(ops_product(contracted_terms(h_pq, expr1))) == Sym("h_{ab}")

    expr1x = Sym("x") * expr1
    @test SymPy.simplify(ops_product(contracted_terms(h_pq, expr1x))) ==
        Sym("h_{ab}") * Sym("x")

    δbp = delta(Sym("b"), Sym("p"))
    δaq = delta(Sym("a"), Sym("q"))
    exprdiff = OpSum(IndicedObj[δap, δbq]) - OpSum(IndicedObj[δbp, δaq])
    @test SymPy.simplify(ops_product(contracted_terms(h_pq, exprdiff))) ==
        Sym("h_{ab}") - Sym("h_{ba}")

    δaαp = delta(Sym("a_{α}"), Sym("p"))
    δbαq = delta(Sym("b_{α}"), Sym("q"))
    exprspin = OpSum(IndicedObj[δaαp, δbαq])
    @test SymPy.simplify(ops_product(contracted_terms(h_pq, exprspin))) ==
        Sym("h_{a_{α}b_{α}}")

    h_pq_spin = one_body_element(orbital_class("p_{α}"), orbital_class("q_{α}"); symbol="h")
    δap_spin = delta(Sym("a_{α}"), Sym("p_{α}"))
    δbq_spin = delta(Sym("b_{α}"), Sym("q_{α}"))
    exprspin2 = OpSum(IndicedObj[δap_spin, δbq_spin])
    @test SymPy.simplify(ops_product(contracted_terms(h_pq_spin, exprspin2))) ==
        Sym("h_{a_{α}b_{α}}")
end
