# N-body Problem Dynamics

_Also known as NBP dynamics!_

```@setup main
using AstrodynamicalModels
using ModelingToolkit
using Latexify
Latexify.auto_display(true)
```

## Overview

In an astrodynamical context, the N-body problem assumes $N$ celestial bodies
which move with respect to some common origin. A body $i$ moves due to the
cumulative gravity of every other body in the system. This problem is
notoriously difficult because it cannot be solved analytically for $N\geq3$!

## Examples

All `NBSystem` calls require the number of bodies to be specified as the first
argument, like so. As always, use the `stm` argument at your leisure. Beware,
though! using `stm=true` for N-body systems with more than 5 bodies may cause
`NBSystem` to compute for a _really, really_ long time! True story – `v1.0.1` of this
package had a version of these docs which tried to compute `NBSystem(30; stm=true)`.
That resulted in GitHub and JuliaHub failing the job on a timeout after _several
hours_ – at first I was surprised, until I realized that appending state
transition matrix dynamics to a 30-body system results in a total of 32580
states!

If you're curious, the number of states of an N-body system with state
transition matrix dynamics appended is equivalent to `N*6 + (N*6)^2`.

```@repl main
model = NBSystem(2; stm=true)
```

Like other models, we can compute the Jacobian for these dynamics.

```@repl main
using SparseArrays
J = sparse(calculate_jacobian(NBSystem(4)))
```

Finally, let's construct a Julia function which implements these dynamics!

```@repl main
f = NBFunction(2)
let u = randn(12), m = randn(2), G = rand(), t = 0
    f(u, [G, m...], t)
end
```
