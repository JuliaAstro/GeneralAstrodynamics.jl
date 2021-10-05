# N-body Problem Dynamics
_Also known as NBP dynamics!_

```@setup main
using AstrodynamicalModels
using ModelingToolkit
using Latexify
Latexify.auto_display(true)
```

## Overview

In an astrodynamical context, the N-body problem assumes 
$N$ celestial bodies which move with respect to some common 
origin. A body $i$ moves due to the cumulative gravity
of every other body in the system. This problem is notoriously
difficult because it cannot be solved analytically for $N>3$!

## Examples

All `NBP` calls require the number of bodies to be specified as 
the first argument, like so. As always, use the `stm` and 
`structural_simplify` arguments at your leisure.

```@repl main
model = NBP(6; stm=true, structural_simplify=false)
```

Like other models, we can compute the Jacobian for these dynamics.

```@repl main
using SparseArrays
J = sparse(calculate_jacobian(NBP(30)))
```

Finally, let's construct a Julia function which implements these dynamics!

```@repl main
f = NBPFunction(2)
let u = randn(12), m = randn(2), G = rand(), t = 0
    f(u, [G, m...], t)
end
```

Note that if you'd like to utilize state transition dynamics, 
or calculate the system's Jacobian,you may benefit from using 
sparcity to describe your dynamics more efficiently.

```@repl main
f = NBPFunction(6; stm=true, jac=true, sparse=true)
let u = randn(6 * 2 * 6), m = randn(20), G = rand(), t = 0
    f(Val{:jac}, u, [G, m...], t)
end
```