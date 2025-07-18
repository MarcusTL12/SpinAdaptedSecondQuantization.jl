using Documenter
using SpinAdaptedSecondQuantization

makedocs(
    sitename="SpinAdaptedSecondQuantization.jl",
    authors="Marcus T. Lexander, Tor S. Haugland, Alexander C. Paul",
    format=Documenter.HTML(),
    modules=[SpinAdaptedSecondQuantization],
    pages=[
        "Home" => "index.md",
        "CCSD example" => "ccsd_example.md",
        "Reference/API" => "reference.md",
        "Internals" => "internals.md",
    ]
)

deploydocs(
    repo="github.com/MarcusTL12/SpinAdaptedSecondQuantization.jl.git",
    push_preview=true
)
