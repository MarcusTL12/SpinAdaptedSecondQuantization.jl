using Documenter
using SpinAdaptedSecondQuantization

makedocs(
    sitename="SpinAdaptedSecondQuantization.jl",
    authors="Marcus T. Lexander, Tor S. Haugland, Federico Rossi, Alexander C. Paul",
    format=Documenter.HTML(),
    modules=[SpinAdaptedSecondQuantization],
    pages=[
        "Getting Started" => "index.md",
        "Examples" => [
            "ccsd_example.md",
            "qed_ccsd_example.md",
            "box13.2.md",
            "code_generation_example.md",
        ],
        "Developer Tutorials" => [
            "new_type_tutorial.md",
            "new_indices_tutorial.md",
        ],
        "Reference/API" => [
            "reference.md",
            "internals.md",
        ]
    ]
)

deploydocs(
    repo="github.com/MarcusTL12/SpinAdaptedSecondQuantization.jl.git",
    push_preview=true
)
