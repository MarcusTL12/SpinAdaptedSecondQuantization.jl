# SpinAdaptedSecondQuantization.jl

[SpinAdaptedSecondQuantization.jl]
(https://github.com/MarcusTL12/SpinAdaptedSecondQuantization.jl)
is a Julia package for doing symbolic second quantization mainly targeted at
quantum chemistry methods, such as post Hartree-Fock methods.

To install the package, run

```julia
(v1.8) pkg> add https://github.com/MarcusTL12/SpinAdaptedSecondQuantization.jl
```

Julia version 1.8 or higher is required.

The package can then be loaded

```julia
using SpinAdaptedSecondQuantization
```

```@setup 1
using SpinAdaptedSecondQuantization
```

```@meta
DocTestSetup = quote
  using SpinAdaptedSecondQuantization
  SASQ.disable_color()
end
```

## Quick Start

Example:

```@example 1
E(1, 2) * occupied(2)
```

Another example:

```@example 1
E(1, 2) * virtual(1)
```
