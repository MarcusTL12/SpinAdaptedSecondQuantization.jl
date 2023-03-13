# Implement ordering of new operator types here:

function Base.isless(::FermionOperator, ::SingletExcitationOperator)
    false
end
