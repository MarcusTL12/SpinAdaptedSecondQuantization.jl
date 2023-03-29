# Implement ordering of new operator types here:

function Base.isless(::FermionOperator, ::SingletExcitationOperator)
    false
end

function Base.isless(::BosonOperator, ::SingletExcitationOperator)
    false
end

function Base.isless(::BosonOperator, ::FermionOperator)
    false
end