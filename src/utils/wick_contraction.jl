using SymPy

@doc raw"""
    isgreek(sigma)

Return true if `sigma` is a single lowercase Greek letter.
"""
function isgreek(sigma)
    s = sigma isa AbstractString ? sigma : string(sigma)
    return occursin(r"^\p{Greek}$", s) && islowercase(only(s))
end


@doc raw"""
    with_spin(label, sigma)

Attach a spin suffix to a label, e.g., "p" + "α" -> "p_{α}".
"""
function with_spin(label, sigma)
    if !isgreek(sigma)
        error("with_spin expects a single lowercase Greek letter for sigma.")
    end
    return Sym("$(string(label))_{$(string(sigma))}")
end


@doc raw"""
    split_orbital_spin(label)

Split a label like "p_{α}" (or "p_α") into ("p", "α"). If no spin suffix is
present, return (label, nothing).
"""
function split_orbital_spin(label)
    s = label isa AbstractString ? label : string(label)
    if s isa AbstractString
        m = match(r"^(.+)_\{(\p{Greek})\}$", s)
        if m !== nothing
            return Sym(m.captures[1]), Sym(m.captures[2])
        end
        parts = split(s, "_")
        if length(parts) == 2 && isgreek(parts[2])
            return Sym(parts[1]), Sym(parts[2])
        end
    end
    return label isa Sym ? label : Sym(s), nothing
end


@doc raw"""
    IndicedObj

Abstract supertype for indexed symbolic objects (e.g., orbitals, operators,
and matrix elements).
"""
abstract type IndicedObj end


@doc raw"""
    Orbital

Container for orbital metadata.
Fields: index, type, spin, label, nocc.
"""
struct Orbital <: IndicedObj
    type::Symbol
    index::Sym
    spin::Union{Sym, Nothing}
    label::Sym
    nocc::Union{Int, Nothing}
end


@doc raw"""
    MatrixElement

Container for integral metadata.
Fields: labels, str, name.
`labels` stores index symbols, `str` stores integral symbol labels (e.g., h/j/k),
and `name` is the full symbolic integral.
"""
struct MatrixElement <: IndicedObj
    labels::Vector{Sym}
    str::Vector{Sym}
    name::Sym
end


@doc raw"""
    orbital_class(label; orb_type=nothing)

Return an Orbital with index, type, and optional spin.
If `orb_type` is `nothing`, type is inferred:
:occ for i/j/k/l, :vir for a/b/c/d, else :gen (e.g., p/q/r/s).
`nocc` is 1 for occupied, 0 for virtual, and `nothing` for general.
General orbitals use no occupation factors in contractions.
"""
function orbital_class(label; orb_type=nothing)
    base, spin = split_orbital_spin(label)
    base_sym = base isa Sym ? base : Sym(string(base))
    full_label = spin === nothing ? base_sym : with_spin(base_sym, spin)
    if orb_type === nothing
        orb_type = :gen
        base_str = string(base_sym)
        if base_str in ("i", "j", "k", "l")
            orb_type = :occ
        elseif base_str in ("a", "b", "c", "d")
            orb_type = :vir
        end
    end
    nocc = orb_type == :occ ? 1 :
        (orb_type == :vir ? 0 : nothing)
    return Orbital(orb_type, base_sym, spin, full_label, nocc)
end


@doc raw"""
    Operator

Container for operator metadata.
Fields: type, str, label, name, nocc.
Note: for creators, `str` is "adag" while `name` is the symbolic a^{†}_{...}.
"""
struct Operator <: IndicedObj
    type::Symbol
    str::Sym
    label::Sym
    name::Sym
    nocc::Union{Int, Nothing}
end


@doc raw"""
    OpSum

Linear combination of operator strings, stored as coefficients and terms.
Each term is a `Vector{IndicedObj}` with an associated `Sym` coefficient.
"""
struct OpSum
    coeffs::Vector{Sym}
    terms::Vector{Vector{IndicedObj}}
end

OpSum() = OpSum([Sym(0)], [IndicedObj[]])
OpSum(ops::AbstractVector{<:IndicedObj}) = OpSum([Sym(1)], [IndicedObj[ops...]])


Base.:(==)(a::Operator, b::Operator) =
    a.type == b.type &&
    a.str == b.str &&
    a.label == b.label &&
    a.name == b.name &&
    a.nocc == b.nocc

Base.:+(a::AbstractVector{<:IndicedObj}, b::AbstractVector{<:IndicedObj}) =
    OpSum([Sym(1), Sym(1)], [IndicedObj[a...], IndicedObj[b...]])
Base.:+(a::OpSum, b::AbstractVector{<:IndicedObj}) =
    OpSum([a.coeffs; Sym(1)], [a.terms; [IndicedObj[b...]]])
Base.:+(a::AbstractVector{<:IndicedObj}, b::OpSum) =
    OpSum([Sym(1); b.coeffs], [IndicedObj[a...]; b.terms])
Base.:+(a::OpSum, b::OpSum) =
    OpSum([a.coeffs; b.coeffs], [a.terms; b.terms])

Base.:-(a::AbstractVector{<:IndicedObj}, b::AbstractVector{<:IndicedObj}) =
    OpSum([Sym(1), Sym(-1)], [IndicedObj[a...], IndicedObj[b...]])
Base.:-(a::OpSum, b::AbstractVector{<:IndicedObj}) =
    OpSum([a.coeffs; Sym(-1)], [a.terms; [IndicedObj[b...]]])
Base.:-(a::AbstractVector{<:IndicedObj}, b::OpSum) =
    OpSum([Sym(1); -b.coeffs], [IndicedObj[a...]; b.terms])
Base.:-(a::OpSum, b::OpSum) =
    OpSum([a.coeffs; -b.coeffs], [a.terms; b.terms])
Base.:-(a::OpSum) = Sym(-1) * a

Base.:*(a::IndicedObj, b::IndicedObj) = IndicedObj[a, b]
Base.:*(ops::AbstractVector{<:IndicedObj}, obj::IndicedObj) = [ops; obj]
Base.:*(obj::IndicedObj, ops::AbstractVector{<:IndicedObj}) = [obj; ops]
Base.:*(ops1::AbstractVector{<:IndicedObj}, ops2::AbstractVector{<:IndicedObj}) =
    [ops1; ops2]

Base.:*(c::Sym, a::OpSum) = OpSum(c .* a.coeffs, a.terms)
Base.:*(a::OpSum, c::Sym) = c * a
# Keep scalar factors from accidentally "infecting" every IndicedObj subtype;
# only Operators and MatrixElements participate directly in OpSum scaling.
Base.:*(c::Sym, obj::Operator) =
    (c == Sym(0) || obj.name == Sym(0)) ? OpSum() : OpSum([c], [IndicedObj[obj]])
Base.:*(c::Sym, obj::MatrixElement) =
    (c == Sym(0) || obj.name == Sym(0)) ? OpSum() : OpSum([c], [IndicedObj[obj]])
Base.:*(obj::Operator, c::Sym) = c * obj
Base.:*(obj::MatrixElement, c::Sym) = c * obj
Base.:*(c::Sym, ops::AbstractVector{<:IndicedObj}) =
    OpSum([c], [IndicedObj[ops...]])
Base.:*(ops::AbstractVector{<:IndicedObj}, c::Sym) = c * ops

Base.:*(a::OpSum, b::IndicedObj) =
    OpSum(a.coeffs, [vcat(t, IndicedObj[b]) for t in a.terms])
Base.:*(a::IndicedObj, b::OpSum) =
    OpSum(b.coeffs, [vcat(IndicedObj[a], t) for t in b.terms])
Base.:*(a::OpSum, b::AbstractVector{<:IndicedObj}) =
    OpSum(a.coeffs, [vcat(t, IndicedObj[b...]) for t in a.terms])
Base.:*(a::AbstractVector{<:IndicedObj}, b::OpSum) =
    OpSum(b.coeffs, [vcat(IndicedObj[a...], t) for t in b.terms])

function Base.:*(a::OpSum, b::OpSum)
    coeffs = Sym[]
    terms = Vector{Vector{IndicedObj}}()
    for (ca, ta) in zip(a.coeffs, a.terms)
        for (cb, tb) in zip(b.coeffs, b.terms)
            push!(coeffs, ca * cb)
            push!(terms, vcat(ta, tb))
        end
    end
    return OpSum(coeffs, terms)
end

function ops_product(os::OpSum)
    total = Sym(0)
    for (coeff, term) in zip(os.coeffs, os.terms)
        prod = Sym(1)
        for obj in term
            if obj isa Operator || obj isa MatrixElement
                prod *= obj.name
            else
                error("ops_product(OpSum) only supports Operator or MatrixElement.")
            end
        end
        total += coeff * prod
    end
    return total
end


@doc raw"""
    ops_product(ops)

Return the Sym product of a list of operators.
"""
function ops_product(ops::Vector{Operator})
    prod = Sym(1)
    for op in ops
        prod *= op.name
    end
    return prod
end


@doc raw"""
    creator(p)

Return a creation operator symbol a^{†}_{p,σ} (or a^{†}_p if spin is omitted).
The `str` field is "adag" while the symbolic `name` uses a^{†}_{...}.
`p` is an `Orbital` instance.
"""
function creator(p::Orbital)
    return Operator(
        :create, Sym("adag"), p.label, Sym("a^{†}_{$(string(p.label))}"), p.nocc
    )
end


@doc raw"""
    annihilator(p)

Return an annihilation operator symbol a_{p,σ} (or a_p if spin is omitted).
`p` is an `Orbital` instance.
"""
function annihilator(p::Orbital)
    return Operator(
        :annihilate, Sym("a"), p.label, Sym("a_{$(string(p.label))}"), p.nocc
    )
end


@doc raw"""
    is_creation(op)

Return true if the operator represents a creation operator.
"""
function is_creation(op::Operator)
    return op.type == :create
end


@doc raw"""
    is_annihilation(op)

Return true if the operator represents an annihilation operator.
"""
function is_annihilation(op::Operator)
    return op.type == :annihilate
end


@doc raw"""
    conjugate(op)

Return the conjugate operator by swapping creator and annihilator.
"""
function conjugate(op::Operator)
    if op.type == :create
        return Operator(
            :annihilate, Sym("a"), op.label, Sym("a_{$(string(op.label))}"), op.nocc
        )
    elseif op.type == :annihilate
        return Operator(
            :create, Sym("adag"), op.label, Sym("a^{†}_{$(string(op.label))}"), op.nocc
        )
    end
    error("Unknown operator type: $(op.type)")
end


@doc raw"""
    conjugate(ops::Vector{Operator})

Return the conjugate of each operator in the list.
"""
function conjugate(ops::Vector{Operator})
    return [conjugate(op) for op in reverse(ops)]
end


@doc raw"""
    fermion_op(kind, p)

Return a creation or annihilation operator by `kind` = :create or :annihilate.
`p` is an `Orbital` instance.
"""
function fermion_op(kind, p::Orbital)
    if kind == :create
        return creator(p)
    elseif kind == :annihilate
        return annihilator(p)
    end
    error("Unknown kind: $(kind). Use :create or :annihilate.")
end


@doc raw"""
    creators(ps)

Return a vector of creation operators for the given orbitals.
`ps` is a vector of `Orbital` instances.
"""
function creators(ps::Vector{Orbital})
    return [creator(p) for p in ps]
end


@doc raw"""
    annihilators(ps)

Return a vector of annihilation operators for the given orbitals.
`ps` is a vector of `Orbital` instances.
"""
function annihilators(ps::Vector{Orbital})
    return [annihilator(p) for p in ps]
end


@doc raw"""
    many_body_operator(ps)

Return a normal-ordered operator list with creators from odd positions
followed by annihilators from even positions in reverse order.
For `ps = [p,q,r,s,...]`, this returns
[a†_p, a†_r, ..., a_s, a_q].
"""
function many_body_operator(ps::Vector{Orbital})
    creators_list = [creator(p) for p in ps[1:2:end]]
    annihilators_list = [annihilator(q) for q in reverse(ps[2:2:end])]
    return vcat(creators_list, annihilators_list)
end


@doc raw"""
    normal_order_op_pairs(ops)

Return a vector of [a†, a] operator pairs for all combinations of creators
and annihilators in the given operator list.
"""
function normal_order_op_pairs(ops::AbstractVector{<:IndicedObj})
    ops = Operator[op::Operator for op in ops]
    creators_list = [op for op in ops if is_creation(op)]
    annihilators_list = [op for op in ops if is_annihilation(op)]
    ops = Vector{Vector{Operator}}()
    for c in creators_list
        for a in annihilators_list
            push!(ops, [c, a])
        end
    end
    return ops
end


@doc raw"""
    one_body_operator(p, q)

Return the one-electron operator list [a^{†}_p, a_q].
`p` and `q` are `Orbital` instances.
"""
function one_body_operator(p::Orbital, q::Orbital)
    return [creator(p), annihilator(q)]
end


@doc raw"""
    two_body_operator(p, q, r, s)

Return the two-electron operator list [a^{†}_p, a^{†}_r, a_s, a_q].
Here p,q label the first particle and r,s label the second particle.
`p,q,r,s` are `Orbital` instances.
"""
function two_body_operator(
    p::Orbital, q::Orbital, r::Orbital, s::Orbital
)
    return [creator(p), creator(r), annihilator(s), annihilator(q)]
end


@doc raw"""
    one_body_element(p, q; symbol="h")

Return h_pq, the one-body integral symbol.
`p` and `q` are `Orbital` instances.
"""
function one_body_element(p::Orbital, q::Orbital; symbol="h")
    labels = [p.label, q.label]
    str = [Sym(symbol)]
    name = Sym("$(symbol)_{$(string(p.label))$(string(q.label))}")
    return MatrixElement(labels, str, name)
end


@doc raw"""
    two_body_element(
        p, q, r, s; swap=false, symbol_j="j", symbol_k="k", c_ex=1
    )

Return g_pqrs = <pr||qs> = (pq|rs) - (ps|rq),
where j_pqrs = (pq|rs) and k_pqrs = (ps|rq).

If `swap=true`, return g_rspq = <rp||sq> (i.e., j_rspq - k_rspq).
The exchange term k_pqrs is scaled by `c_ex`.
Here p,q label the first particle and r,s label the second particle.
`p,q,r,s` are `Orbital` instances.
"""
function two_body_element(
    p::Orbital, q::Orbital, r::Orbital, s::Orbital;
    swap=false, symbol_j="j", symbol_k="k", c_ex=1
)
    if swap
        return two_body_element(
            r, s, p, q; swap=false, symbol_j=symbol_j, symbol_k=symbol_k,
            c_ex=c_ex
        )
    end
    labels = [p.label, q.label, r.label, s.label]
    str = [Sym(symbol_j), Sym(symbol_k)]
    j_pqrs = Sym(
        "$(symbol_j)_{$(string(p.label))$(string(q.label))$(string(r.label))$(string(s.label))}"
    )
    k_pqrs = Sym(
        "$(symbol_k)_{$(string(p.label))$(string(q.label))$(string(r.label))$(string(s.label))}"
    )
    name = j_pqrs - Sym(c_ex) * k_pqrs
    return MatrixElement(labels, str, name)
end


@doc raw"""
    one_body_term(p, q; symbol="h")

Return the normal-ordered one-electron Hamiltonian term: h_pq a^{†}_p a_q.
`p` and `q` are `Orbital` instances.
"""
function one_body_term(p::Orbital, q::Orbital; symbol="h")
    h_pq = one_body_element(p, q; symbol=symbol)
    return h_pq.name * ops_product(one_body_operator(p, q))
end


@doc raw"""
    two_body_term(
        p, q, r, s; swap=false, symmetric=false, symbol_j="j", symbol_k="k",
        c_ex=1
    )

Return the normal-ordered two-electron Hamiltonian term:
(1/4) g_pqrs a^{†}_p a^{†}_r a_s a_q, or (1/2) g_pqrs a^{†}_p a^{†}_r a_s a_q
if `symmetric=true`.
The exchange term k_pqrs in g_pqrs is scaled by `c_ex`.
Here p,q label the first particle and r,s label the second particle.
`p,q,r,s` are `Orbital` instances.
"""
function two_body_term(
    p::Orbital, q::Orbital, r::Orbital, s::Orbital;
    swap=false, symmetric=false, symbol_j="j", symbol_k="k", c_ex=1
)
    g_pqrs = two_body_element(
        p, q, r, s; swap=swap, symbol_j=symbol_j, symbol_k=symbol_k, c_ex=c_ex
    )
    factor = symmetric ? Sym(1//2) : Sym(1//4)
    return factor * g_pqrs.name * ops_product(two_body_operator(p, q, r, s))
end


@doc raw"""
    Hamiltonian

Container for a normal-ordered Hamiltonian with one- and two-electron terms.
Fields: h_pq, g_pqrs, c_ex, one_body, two_body.
"""
struct Hamiltonian
    h_pq::MatrixElement
    g_pqrs::MatrixElement
    c_ex::Sym
    one_body::Vector{Operator}
    two_body::Vector{Operator}
end


@doc raw"""
    hamiltonian(p, q, r, s; c_ex=1)

Build a Hamiltonian with h_pq, g_pqrs, exact exchange factor c_ex,
and operators a^{†}_p a_q, a^{†}_p a^{†}_r a_s a_q.
`p,q,r,s` are `Orbital` instances.
"""
function hamiltonian(
    p::Orbital, q::Orbital, r::Orbital, s::Orbital; c_ex=1
)
    h_pq = one_body_element(p, q)
    g_pqrs = two_body_element(p, q, r, s; c_ex=c_ex)
    one_body = one_body_operator(p, q)
    two_body = two_body_operator(p, q, r, s)
    return Hamiltonian(h_pq, g_pqrs, Sym(c_ex), one_body, two_body)
end


@doc raw"""
    one_body_excitation(i, a; reverse=false)

Return the one-electron excitation operator list [a^{†}_a, a_i].
If `reverse=true`, return [a^{†}_i, a_a] instead for de-excitation.
`i` and `a` are `Orbital` instances.
"""
function one_body_excitation(i::Orbital, a::Orbital; reverse=false)
    if i.type != :occ || a.type != :vir
        error("one_body_excitation expects i occupied and a virtual.")
    end
    if reverse
        return one_body_operator(i, a)
    end
    return one_body_operator(a, i)
end


@doc raw"""
    two_body_excitation(i, j, a, b)

Return the two-electron excitation operator list [a^{†}_a, a^{†}_b, a_j, a_i].
`i,j,a,b` are `Orbital` instances.
"""
function two_body_excitation(
    i::Orbital, j::Orbital, a::Orbital, b::Orbital
)
    if i.type != :occ || j.type != :occ || a.type != :vir || b.type != :vir
        error("two_body_excitation expects i,j occupied and a,b virtual.")
    end
    return two_body_operator(a, i, b, j)
end


@doc raw"""
    one_body_excitation_spin_conserve(i, a; sigma="α")

Return spin-conserving excitation E_ai(σ,σ) = a^{†}_{a,σ} a_{i,σ}.
`i` and `a` are `Orbital` instances.
"""
function one_body_excitation_spin_conserve(
    i::Orbital, a::Orbital; sigma="α"
)
    return one_body_excitation(orbital_class(with_spin(i.label, sigma)),
        orbital_class(with_spin(a.label, sigma)))
end


@doc raw"""
    one_body_excitation_spin_flip(i, a; sigma_from="α", sigma_to="β")

Return spin-flip excitation E_ai(σ_to,σ_from) = a^{†}_{a,σ_to} a_{i,σ_from}.
`i` and `a` are `Orbital` instances.
"""
function one_body_excitation_spin_flip(
    i::Orbital, a::Orbital; sigma_from="α", sigma_to="β"
)
    return one_body_excitation(orbital_class(with_spin(i.label, sigma_from)),
        orbital_class(with_spin(a.label, sigma_to)))
end


@doc raw"""
    spin_state_excitation(i, a)

Return the singlet and triplet components:
S0 = (E_ai(α,α) + E_ai(β,β)) / sqrt(2),
T+ = E_ai(α,β), T- = E_ai(β,α),
T0 = (E_ai(α,α) - E_ai(β,β)) / sqrt(2).
`i` and `a` are `Orbital` instances. Each component is an `OpSum`.
"""
function spin_state_excitation(i::Orbital, a::Orbital)
    coeff = Sym(1) / SymPy.sqrt(2)
    e_alpha = one_body_excitation_spin_conserve(i, a; sigma="α")
    e_beta = one_body_excitation_spin_conserve(i, a; sigma="β")
    s_zero = coeff * (OpSum(e_alpha) + OpSum(e_beta))
    t_plus = Sym(1) *
        one_body_excitation_spin_flip(i, a; sigma_from="β", sigma_to="α")
    t_minus = Sym(1) *
        one_body_excitation_spin_flip(i, a; sigma_from="α", sigma_to="β")
    t_zero = coeff * (OpSum(e_alpha) - OpSum(e_beta))
    return (s_zero=s_zero, t_plus=t_plus, t_minus=t_minus, t_zero=t_zero)
end


@doc raw"""
    ExcitedState

Parent type for excited-state containers.
"""
abstract type ExcitedState end


@doc raw"""
    SingleExcitationState

Container for a single-excitation state with amplitudes x_ai, y_ai and
operator a^{†}_a a_i.
Fields: x_ai, y_ai, t_ai, op.
"""
struct SingleExcitationState <: ExcitedState
    x_ai::Sym
    y_ai::Sym
    t_ai::Sym
    op1::Vector{Operator}
end


@doc raw"""
    DoubleExcitationState

Container for a double-excitation state with amplitudes t_ai, t_aibj and
operator a^{†}_a a^{†}_b a_j a_i.
Fields: t_ai, t_aibj, op.
"""
struct DoubleExcitationState <: ExcitedState
    t_ai::Sym
    t_aibj::Sym
    op1::Vector{Operator}
    op2::Vector{Operator}
end


@doc raw"""
    single_excitation_state(i, a)

Build a single-excitation state with amplitudes x_ai, y_ai and operator
a†_a a_i.
`i` and `a` are `Orbital` instances.
"""
function single_excitation_state(i::Orbital, a::Orbital)
    x_ai = Sym("x^$(string(a.label))_$(string(i.label))")
    y_ai = Sym("y^$(string(a.label))_$(string(i.label))")
    t_ai = x_ai
    op1 = one_body_excitation(i, a)
    return SingleExcitationState(x_ai, y_ai, t_ai, op1)
end


@doc raw"""
    double_excitation_state(i, j, a, b)

Build a double-excitation state with amplitude t_aibj and
operator a†_a a†_b a_j a_i.
`i,j,a,b` are `Orbital` instances.
"""
function double_excitation_state(
    i::Orbital, j::Orbital, a::Orbital, b::Orbital
)
    t_ai = Sym("t^$(string(a.label))_$(string(i.label))")
    t_aibj = Sym("t^$(string(a.label))$(string(b.label))_$(string(i.label))$(string(j.label))")
    op1 = one_body_excitation(i, a)
    op2 = two_body_excitation(i, j, a, b)
    return DoubleExcitationState(t_ai, t_aibj, op1, op2)
end


@doc raw"""
    delta(p::Sym, q::Sym)

Return the delta matrix element δ_pq. `p` and `q` may be orbital labels, spin
labels, or spin-orbital labels.
"""
function delta(p::Sym, q::Sym)
    ps = string(p)
    qs = string(q)
    δsym = Sym("δ")
    if ps == "α" || ps == "β" || qs == "α" || qs == "β"
        if ps == qs
            return MatrixElement(Sym[], Sym[δsym], Sym(1))
        end
        return MatrixElement(Sym[], Sym[δsym], Sym(0))
    end

    p_spin = endswith(ps, "_{α}") ? "α" : (endswith(ps, "_{β}") ? "β" : nothing)
    q_spin = endswith(qs, "_{α}") ? "α" : (endswith(qs, "_{β}") ? "β" : nothing)
    if p_spin !== nothing && q_spin !== nothing && p_spin != q_spin
        return MatrixElement(Sym[], Sym[δsym], Sym(0))
    end

    name = Sym("δ_{$ps$qs}")
    return MatrixElement(Sym[p, q], Sym[δsym], name)
end


@doc raw"""
    delta(op1::Operator, op2::Operator)

Return the Kronecker delta matrix element δ_pq for operator labels.
If the pair is not a creation/annihilation (either order), return 0.
"""
function delta(op1::Operator, op2::Operator)
    if !((is_creation(op1) && is_annihilation(op2)) ||
        (is_annihilation(op1) && is_creation(op2)))
        return MatrixElement(Sym[], Sym[Sym("δ")], Sym(0))
    end
    return delta(op1.label, op2.label)
end


@doc raw"""
    contraction(op1, op2)

Return the fermionic contraction with occupation factors:
if the creator is left of annihilator, use n_p;
if the creator is right of annihilator, use (1 - n_p).
`op1` and `op2` are `Operator` instances.
General orbitals (nocc = nothing) use no occupation factors.
"""
function contraction(op1::Operator, op2::Operator)
    if op1.name == Sym(0) || op2.name == Sym(0)
        return MatrixElement(Sym[], Sym[Sym("δ")], Sym(0))
    end
    n1 = op1.nocc === nothing ? Sym(1) :
        (is_creation(op1) ? Sym(op1.nocc) : Sym(1) - Sym(op1.nocc))
    n2 = op2.nocc === nothing ? Sym(1) :
        (is_annihilation(op2) ? Sym(op2.nocc) : Sym(1) - Sym(op2.nocc))
    n = n1 * n2
    if n == Sym(0)
        return MatrixElement(Sym[], Sym[Sym("δ")], Sym(0))
    end
    d = delta(op1, op2)
    return MatrixElement(d.labels, d.str, n * d.name)
end


@doc raw"""
    contraction_sign(ops, pairs)

Compute fermionic sign by swapping paired operators together and removing them.
"""
function contraction_sign(ops::Vector{Operator}, pairs::Vector{Tuple{Int, Int}})
    sign = Sym(1)
    positions = collect(1:length(ops))
    for (i, j) in pairs
        i_pos = findfirst(==(i), positions)
        j_pos = findfirst(==(j), positions)
        if j_pos < i_pos
            i_pos, j_pos = j_pos, i_pos
        end
        swaps = j_pos - i_pos - 1
        if isodd(swaps)
            sign = -sign
        end
        deleteat!(positions, j_pos)
        deleteat!(positions, i_pos)
    end
    return sign
end


@doc raw"""
    normal_order(ops; sign=Sym(1))

Return the normal-ordered operator list and sign by moving creation operators
left.
"""
function normal_order(ops::Vector{Operator}; sign=Sym(1))
    ops_copy = copy(ops)
    for i in 1:length(ops_copy)
        for j in i:-1:2
            if is_creation(ops_copy[j]) && is_annihilation(ops_copy[j - 1])
                ops_copy[j - 1], ops_copy[j] = ops_copy[j], ops_copy[j - 1]
                sign = -sign
            end
        end
    end
    return sign * ops_copy
end


@doc raw"""
    wick_pairs(ops; forbidden_pairs=Vector{Vector{Operator}}())

Return all creator/annihilator index pairs (i,j) with i < j.
"""
function wick_pairs(
    ops::Vector{Operator};
    forbidden_pairs::Vector{Vector{Operator}}=Vector{Vector{Operator}}()
)
    pairs = Tuple{Int, Int}[]
    n = length(ops)
    for i in 1:n
        for j in i + 1:n
            if (is_creation(ops[i]) && is_annihilation(ops[j])) ||
                (is_annihilation(ops[i]) && is_creation(ops[j]))
                skip = false
                for pair in forbidden_pairs
                    if length(pair) == 2 &&
                        ((ops[i] == pair[1] && ops[j] == pair[2]) ||
                        (ops[i] == pair[2] && ops[j] == pair[1]))
                        skip = true
                        break
                    end
                end
                if skip
                    continue
                end
                push!(pairs, (i, j))
            end
        end
    end
    return pairs
end


@doc raw"""
    pairings_for_order(pairs, k)

Return all k-pair matchings from a list of pairs, with no shared indices.
"""
function pairings_for_order(pairs::Vector{Tuple{Int, Int}}, k::Int)
    if k == 0
        return [Tuple{Int, Int}[]]
    end
    if length(pairs) < k
        return Vector{Vector{Tuple{Int, Int}}}()
    end
    results = Vector{Vector{Tuple{Int, Int}}}()
    first = pairs[1]
    rest = pairs[2:end]
    (i1, j1) = first
    remaining = [p for p in rest if !(i1 in p) && !(j1 in p)]
    for tail in pairings_for_order(remaining, k - 1)
        push!(results, [first; tail])
    end
    for tail in pairings_for_order(rest, k)
        push!(results, tail)
    end
    return results
end


@doc raw"""
    wick_contractions_order(ops, k)

Return the sum of all k-contraction terms with the remaining operators
normal-ordered.
"""
function wick_contractions_order(
    ops::Vector{Operator}, k::Int, pairs::Vector{Tuple{Int, Int}}
)
    total = OpSum()
    n = length(ops)
    if k == 0
        return normal_order(ops)
    end
    if n < 2k
        return total
    end
    for pairs_k in pairings_for_order(pairs, k)
        c = OpSum(IndicedObj[])  # identity OpSum
        zero = false
        for (i, j) in pairs_k
            cij = contraction(ops[i], ops[j])
            if cij.name == Sym(0)
                zero = true
                break
            end
            c *= cij
        end
        if zero
            continue
        end
        remaining = Operator[]
        contracted = Set(vcat([p[1] for p in pairs_k], [p[2] for p in pairs_k]))
        for idx in 1:n
            if !(idx in contracted)
                push!(remaining, ops[idx])
            end
        end
        sign = contraction_sign(ops, pairs_k)
        rem_ordered = normal_order(remaining)
        total += (sign * c) * rem_ordered
    end
    return total
end


@doc raw"""
    wick_expand(ops)

Return the normal-ordered product plus contractions of all orders
up to min(number of creators, number of annihilators).
"""
function wick_expand(
    ops::Vector{Operator};
    forbidden_pairs::Vector{Vector{Operator}}=Vector{Vector{Operator}}()
)
    n_cre = count(is_creation, ops)
    n_ann = count(is_annihilation, ops)
    max_pairs = min(n_cre, n_ann)
    pairs = wick_pairs(ops; forbidden_pairs=forbidden_pairs)
    total = OpSum()
    for k in 0:max_pairs
        total += wick_contractions_order(ops, k, pairs)
    end
    return total
end


@doc raw"""
    wick_full_contraction(ops)

Return the fully contracted term when the number of creators equals the
number of annihilators. If they differ, return 0.
"""
function wick_full_contraction(
    ops::Vector{Operator};
    forbidden_pairs::Vector{Vector{Operator}}=Vector{Vector{Operator}}()
)
    n_cre = count(is_creation, ops)
    n_ann = count(is_annihilation, ops)
    if n_cre != n_ann
        return OpSum()
    end
    pairs = wick_pairs(ops; forbidden_pairs=forbidden_pairs)
    return wick_contractions_order(ops, n_cre, pairs)
end


@doc raw"""
    contracted_terms(integral, expr)

Given a matrix element (e.g., h_{pq} or g_{pqrs}) and a Wick contraction
expression, replace integral indices using Kronecker deltas and remove the
delta factors. If a delta introduces multiple replacements for the same index,
an error is raised. Returns an `OpSum` (use `ops_product` to obtain a `Sym`).
"""
function contracted_terms(integral::MatrixElement, expr::OpSum)
    if isempty(integral.str)
        error("contracted_terms expects a matrix element with a symbol label.")
    end
    prefix = integral.str[1]
    indices = integral.labels
    index_set = Set(indices)
    δsym = Sym("δ")

    function collect_replacements(term::Vector{IndicedObj})
        replacements = Dict{Sym, Sym}()
        indexed = Set{Int}()
        for (j, obj) in enumerate(term)
            if obj isa MatrixElement && !isempty(obj.str) && obj.str[1] == δsym
                if length(obj.labels) != 2
                    error("Delta must have exactly two indices.")
                end
                a, b = obj.labels
                touched = false
                if a in index_set
                    if haskey(replacements, a) && replacements[a] != b
                        error("Multiple replacements for index $(a).")
                    end
                    replacements[a] = b
                    touched = true
                end
                if b in index_set
                    if haskey(replacements, b) && replacements[b] != a
                        error("Multiple replacements for index $(b).")
                    end
                    replacements[b] = a
                    touched = true
                end
                if touched
                    push!(indexed, j)
                end
            end
        end

        others = OpSum(IndicedObj[])  # identity OpSum
        for (j, obj) in enumerate(term)
            if !(j in indexed)
                others *= obj
            end
        end
        return replacements, others
    end

    total = OpSum()
    for (coeff, term) in zip(expr.coeffs, expr.terms)
        replacements, others = collect_replacements(term)
        replaced = Sym[]
        for idx in indices
            if haskey(replacements, idx)
                push!(replaced, replacements[idx])
            else
                push!(replaced, idx)
            end
        end
        new_name = string(prefix, "_{", join(string.(replaced), ""), "}")
        new_integral = MatrixElement(replaced, integral.str, Sym(new_name))
        total += coeff * new_integral * others
    end
    return total
end



@doc raw"""
    commutator(A, B)

Return the commutator [A, B] = A*B - B*A.
"""
function commutator(A, B)
    return A * B - B * A
end


@doc raw"""
    double_commutator(A, B, C)

Return the average of [[A, B], C] and [A, [B, C]].
"""
function double_commutator(A, B, C)
    term1 = commutator(commutator(A, B), C)
    term2 = commutator(A, commutator(B, C))
    return Sym(1//2) * (term1 + term2)
end


if abspath(PROGRAM_FILE) == @__FILE__
    println("Running wick_contraction.jl example")

    i, j = orbital_class("i_{α}"), orbital_class("j_{α}")
    a, b = orbital_class("a_{α}"), orbital_class("b_{α}")
    p, q = orbital_class("p"), orbital_class("q")
    r, s = orbital_class("r"), orbital_class("s")

    t_ai = one_body_excitation(i, a)
    t_bj = one_body_excitation(j, b)
    h_pq = one_body_element(p, q; symbol="h")
    h_op = one_body_operator(p, q)
    g_pqrs = two_body_element(p, q, r, s; symbol_j="g", c_ex=0)
    g_op = many_body_operator([p, q, r, s])
    filter = normal_order_op_pairs(g_op)

    ops = conjugate(t_ai) * h_op * t_bj
    wick_term = wick_full_contraction(ops; forbidden_pairs=filter)
    wick_term_sym = ops_product(wick_term)
    wick_simplified = SymPy.simplify(wick_term_sym)
    wick_latex = SymPy.latex(wick_simplified)
    contracted_h = contracted_terms(h_pq, wick_term)
    contracted_h_simplified = SymPy.simplify(ops_product(contracted_h))
    contracted_h_latex = SymPy.latex(contracted_h_simplified)

    ops2 = conjugate(t_ai) * g_op * t_bj
    wick_term2 = wick_full_contraction(ops2; forbidden_pairs=filter)
    wick_term2_sym = ops_product(wick_term2)
    wick_simplified2 = SymPy.simplify(wick_term2_sym)
    wick_latex2 = SymPy.latex(wick_simplified2)
    contracted_g = contracted_terms(g_pqrs, wick_term2)
    contracted_g_simplified = SymPy.simplify(ops_product(contracted_g))
    contracted_g_latex = SymPy.latex(contracted_g_simplified)

    t_ai_sym = SymPy.simplify(ops_product(t_ai))
    t_bj_sym = SymPy.simplify(ops_product(t_bj))
    h_op_sym = SymPy.simplify(ops_product(h_op))

    println("t_ai: ", t_ai_sym)
    println("t_bj: ", t_bj_sym)
    println("h_pq a^{†}_p a_q: ", h_pq.name, " * ", h_op_sym)
    println("Wick full contraction: ", wick_term_sym)
    println("Contracted h_pq (simplified): ", contracted_h_simplified)
    println("Contracted h_pq (LaTeX): ", contracted_h_latex)
    println("Wick full contraction (simplified): ", wick_simplified)
    println("Wick full contraction (LaTeX): ", wick_latex)
    println("g_pqrs a^{†}_p a^{†}_r a_s a_q: ", g_pqrs.name, " * ",
        ops_product(g_op))
    println("Two-electron full contraction: ", wick_term2_sym)
    println("Contracted g_pqrs (simplified): ", contracted_g_simplified)
    println("Contracted g_pqrs (LaTeX): ", contracted_g_latex)
    println("Two-electron full contraction (simplified): ", wick_simplified2)
    println("Two-electron full contraction (LaTeX): ", wick_latex2)
end
