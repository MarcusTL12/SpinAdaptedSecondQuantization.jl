export boson, bosondag

"""
    Boson Operators

The basic boson type operator.
"""
struct BosonOperator <: Operator
    dag::Bool
end

function Base.show(io::IO,
    (b, _, _)::Tuple{BosonOperator,Constraints,IndexTranslation})
    dag = b.dag ? '†' : '⁻'
    print(io, 'b', dag)
end

function exchange_indices(b::BosonOperator, mapping)
    b
end

function get_all_indices(::BosonOperator)
    ()
end

function Base.isless(a::BosonOperator, b::BosonOperator)
    a.dag < b.dag
end

function Base.:(==)(a::BosonOperator, b::BosonOperator)
    a.dag == b.dag
end

"""
    boson()

Constructs a boson annihilation operator.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> boson()
b⁻
```
"""
boson() = Expression(BosonOperator(false))

"""
    bosondag()

Constructs a boson creation operator.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> bosondag()
b†
julia> boson()'
b†
```
"""
bosondag() = Expression(BosonOperator(true))

function act_on_ket(op::BosonOperator)
    (op.dag) * Expression(op)
end

function Base.adjoint(op::BosonOperator)
    BosonOperator(!op.dag)
end
