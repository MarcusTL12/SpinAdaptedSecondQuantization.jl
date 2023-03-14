using Documenter
using SpinAdaptedSecondQuantization

makedocs(
    sitename = "SpinAdaptedSecondQuantization.jl",
    format = Documenter.HTML(),
    modules = [SpinAdaptedSecondQuantization]
)

deploydocs(repo = "github.com/MarcusTL12/SpinAdaptedSecondQuantization.jl.git")
