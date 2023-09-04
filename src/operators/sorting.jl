# Implement ordering of new operator types here:

# current sorting is: T, E, a, b

function Base.isless(::FermionOperator, ::SingletExcitationOperator)
    false
end

function Base.isless(::BosonOperator, ::SingletExcitationOperator)
    false
end

function Base.isless(::BosonOperator, ::FermionOperator)
    false
end

function Base.isless(::SingletExcitationOperator, ::TripletExcitationOperator)
    false
end

function Base.isless(::FermionOperator, ::TripletExcitationOperator)
    false
end

function Base.isless(::BosonOperator, ::TripletExcitationOperator)
    false
end
