# Implement ordering of new operator types here:

# current sorting is: T, E, e, a, b

function Base.isless(::A, ::B) where {A<:Operator,B<:Operator}
    A < B
end

Base.isless(::Type{BosonOperator}, ::Type{FermionOperator}) = false
Base.isless(::Type{BosonOperator}, ::Type{SingletDoubleExcitationOperator}) = false
Base.isless(::Type{BosonOperator}, ::Type{SingletExcitationOperator}) = false
Base.isless(::Type{BosonOperator}, ::Type{TripletExcitationOperator}) = false

Base.isless(::Type{FermionOperator}, ::Type{SingletDoubleExcitationOperator}) = false
Base.isless(::Type{FermionOperator}, ::Type{SingletExcitationOperator}) = false
Base.isless(::Type{FermionOperator}, ::Type{TripletExcitationOperator}) = false

Base.isless(::Type{SingletDoubleExcitationOperator}, ::Type{SingletExcitationOperator}) = false
Base.isless(::Type{SingletDoubleExcitationOperator}, ::Type{TripletExcitationOperator}) = false

Base.isless(::Type{SingletExcitationOperator}, ::Type{TripletExcitationOperator}) = false

