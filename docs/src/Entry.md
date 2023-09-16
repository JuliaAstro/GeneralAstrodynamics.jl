# Planar Entry Dynamics

_Also known as canonical entry dynamics!_

```@setup main
using AstrodynamicalModels
using ModelingToolkit
using Latexify
Latexify.auto_display(true)
```

## Overview

The Planar Entry model assumes a spacecraft moving in an exponential atmosphere
about a spherical planet. Acceleration due to gravity is ignored. The equations
of motion are shown below.

$$\begin{aligned}
  \dot{\gamma} &= \frac{1}{v} \left( L_m - (1 - \frac{v^2}{v_c^2}) g \cos{\gamma} \right) \\
  \dot{v} &= -D_m - g \sin{\gamma} \\
  \dot{r} &= v \sin{\gamma} \\
  \dot{\theta} &= \frac{v}{r} \cos{\gamma} \\
\end{aligned}$$

## Examples

```@repl main
model = PlanarEntry()
```

Let's compute the Jacobian for these dynamics.

```@repl main
J = calculate_jacobian(PlanarEntry())
```

Finally, let's construct a Julia function which implements these dynamics!

```@repl main
f = PlanarEntryFunction()
let u = abs.(randn(4)), p = abs.(randn(7)), t = 0
    f(u, p, t)
end
```
