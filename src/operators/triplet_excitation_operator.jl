export τ

"""
    SingletExcitationOperator

The T_pq type operator.
"""
struct TripletExcitationOperator <: Operator
    p::Int
    q::Int
end

function Base.show(io::IO,
    (
        e, constraints, translation
    )::Tuple{TripletExcitationOperator,Constraints,IndexTranslation})
    print(io, "T_")
    print_mo_index(io, constraints, translation, e.p, e.q)
end

function print_latex(io::IO,
    (
        e, constraints, translation
    )::Tuple{TripletExcitationOperator,Constraints,IndexTranslation})
    print(io, "T_{")
    print_latex_mo_index(io, constraints, translation, e.p, e.q)
    print(io, "}")
end

function exchange_indices(e::TripletExcitationOperator, mapping)
    TripletExcitationOperator(
        exchange_index(e.p, mapping),
        exchange_index(e.q, mapping)
    )
end

function get_all_indices(e::TripletExcitationOperator)
    (e.p, e.q)
end

function Base.:(==)(a::TripletExcitationOperator, b::TripletExcitationOperator)
    (a.p, a.q) == (b.p, b.q)
end

function Base.isless(a::TripletExcitationOperator, b::TripletExcitationOperator)
    (a.p, a.q) < (b.p, b.q)
end

"""
    τ(p, q)

Constructs an expression containing a single triplet excitation operator.
"""
τ(p, q) = Expression(TripletExcitationOperator(p, q))

function convert_to_elementary_operators(o::TripletExcitationOperator)
    Expression(
        [
            (fermiondag(o.p, α)*fermion(o.q, α))[1]
            -(fermiondag(o.p, β)*fermion(o.q, β))[1]
        ]
    )
end

function act_on_ket(op::TripletExcitationOperator)
    p = op.p
    q = op.q
    τ(p, q) * virtual(p) * occupied(q)
end

function Base.adjoint(op::TripletExcitationOperator)
    TripletExcitationOperator(op.q, op.p)
end
