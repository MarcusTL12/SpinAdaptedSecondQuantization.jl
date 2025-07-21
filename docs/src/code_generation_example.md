
# Generating runnable code from symbolic expressions

Here we show a small example of how to translate a derived expression into
runnable code. For this we have the `print_julia_function` function which
translates an expression into a julia function that can be run using the
[`TensorOperations.jl`](https://github.com/Jutho/TensorOperations.jl)
package.

As an example we load the `omega_aibj_r` expression from the CCSD example.

```@repl 1
using SpinAdaptedSecondQuantization
using Serialization
omega_aibj_r = deserialize("omega_aibj_r")
```

Then we can turn this into a function:

```@repl 1
print_julia_function("omega_doubles_r", omega_aibj_r);
println(ans)
```

By default the generated function will take in tensors with all indices
assumed to run over all molecular orbitals and make views into these to get
sub-blocks where some indices run over only occupied or virtual indices. For
integral tensors this can be nice if many different blocks are used, then
the caller of the function does not have to think about which blocks to give
to the function. However, for tensors such as the T amplitudes where only one
specific block exists (the ovov...) block, it is better to only have to give
this, rather than embedding it in a much larger tensor. To achieve this one
can specify the optional keyword argument `explicit_tensor_blocks` which
makes all sub-blocks used of these tensors be specified as individual
inputs to the generated function.

```@repl 1
print_julia_function("omega_doubles_r", omega_aibj_r;
    explicit_tensor_blocks = ["t", "u"]);
println(ans)
```

After generating the function, it can be copied into a new julia session

```@repl 2
using TensorOperations
function omega_doubles_r(no, nv, F, t_vovo, g, u_vovo)
    nao = no + nv
    o = 1:no
    v = no+1:nao
    
    X = zeros(nv, no, nv, no)
    g_vvoo = @view g[v,v,o,o]
    g_voov = @view g[v,o,o,v]
    g_ovov = @view g[o,v,o,v]
    F_vv = @view F[v,v]
    F_oo = @view F[o,o]
    @tensoropt (a=>10χ,b=>10χ,c=>10χ,d=>10χ,i=>χ,j=>χ,k=>χ,l=>χ) begin
        X[a,i,b,j] += F_vv[a,c]*t_vovo[b,j,c,i]
        X[a,i,b,j] += -F_oo[k,i]*t_vovo[a,k,b,j]
        X[a,i,b,j] += -g_vvoo[a,c,k,i]*t_vovo[b,j,c,k]
        X[a,i,b,j] += -g_vvoo[a,c,k,j]*t_vovo[b,k,c,i]
        X[a,i,b,j] += g_voov[a,i,k,c]*u_vovo[b,j,c,k]
        X[a,i,b,j] += -g_ovov[k,c,l,d]*t_vovo[a,i,b,k]*u_vovo[c,j,d,l]
        X[a,i,b,j] += -g_ovov[k,c,l,d]*t_vovo[a,i,c,j]*u_vovo[b,k,d,l]
    end
    X
end
```

Then to test it we can setup some test tensors.

```@repl 2
no = 4 # Number of occupied orbitals
nv = 10 # Number of virtual orbitals
nmo = no + nv # Number of molecular orbitals
F = randn(nmo, nmo) # Fock matrix populated by random numbers
g = randn(nmo, nmo, nmo, nmo);
t = randn(nv, no, nv, no);
u = 2 * t - PermutedDimsArray(t, (1, 4, 3, 2));
```

Then we can call the function and post-process the output by symmetrizing and
scaling the diagonal

```@repl 2
Ω_aibj = @time omega_doubles_r(no, nv, F, t, g, u);
Ω_aibj += PermutedDimsArray(t, (3, 4, 1, 2)); # Symmetrize
for a in 1:nv, i in 1:no
    Ω_aibj[a, i, a, i] *= 0.5 # Scale diagonal by 0.5
end
```
