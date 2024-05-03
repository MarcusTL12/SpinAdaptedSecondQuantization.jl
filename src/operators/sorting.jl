# Implement ordering of new operator types here:

# current sorting is: T, E, a, Q, b

function Base.isless(::A, ::B) where {A<:Operator,B<:Operator}
    A < B
end

function Base.isless(::Type{FermionOperator}, ::Type{SingletExcitationOperator})
    false
end

function Base.isless(::Type{BosonOperator}, ::Type{SingletExcitationOperator})
    false
end

function Base.isless(::Type{BosonOperator}, ::Type{FermionOperator})
    false
end

function Base.isless(::Type{SingletExcitationOperator}, ::Type{TripletExcitationOperator})
    false
end

function Base.isless(::Type{FermionOperator}, ::Type{TripletExcitationOperator})
    false
end

function Base.isless(::Type{BosonOperator}, ::Type{TripletExcitationOperator})
    false
end

function Base.isless(::Type{PairFermionOperator}, ::Type{TripletExcitationOperator})
    false
end

function Base.isless(::Type{PairFermionOperator}, ::Type{SingletExcitationOperator})
    false
end

function Base.isless(::Type{PairFermionOperator}, ::Type{FermionOperator})
    false
end

function Base.isless(::Type{BosonOperator}, ::Type{PairFermionOperator})
    false
end
