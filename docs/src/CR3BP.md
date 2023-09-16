# Circular Restricted Three-body Dynamics
_Also known as CR3BP dynamics!_

```@setup main
using AstrodynamicalModels
using ModelingToolkit
using Latexify
Latexify.auto_display(true)
```

## Overview

The Circular Restricted Three-body Problem (CR3BP) assumes a massless spacecraft which moves
due to the gravity of __two__ celestial bodies which orbit their common center of mass.
This may seem like an arbitrary model, but it's actually a pretty decent approximation for how 
a spacecraft moves nearby the Earth and the Sun, the Earth and the Moon, the Sun and 
Jupiter, and other systems in our solar system! The equations of motion 
are provided below.

$\begin{aligned}
\frac{dx(t)}{dt} =& ẋ\left( t \right) \\
\frac{dy(t)}{dt} =& ẏ\left( t \right) \\
\frac{dz(t)}{dt} =& ż\left( t \right) \\
\frac{dẋ(t)}{dt} =& 2 ẏ\left( t \right) - \left( \frac{1}{\sqrt{\left( \mu + x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}}} \right)^{3} \left( 1 - \mu \right) \left( \mu + x\left( t \right) \right) - \left( \frac{1}{\sqrt{\left( -1 + \mu + x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}}} \right)^{3} \mu \left( -1 + \mu + x\left( t \right) \right) + x\left( t \right) \\
\frac{dẏ(t)}{dt} =&  - 2 ẋ\left( t \right) - \left( \left( \frac{1}{\sqrt{\left( \mu + x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}}} \right)^{3} \left( 1 - \mu \right) + \left( \frac{1}{\sqrt{\left( -1 + \mu + x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}}} \right)^{3} \mu \right) y\left( t \right) + y\left( t \right) \\
\frac{dż(t)}{dt} =& \left(  - \left( \frac{1}{\sqrt{\left( \mu + x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}}} \right)^{3} \left( 1 - \mu \right) - \left( \frac{1}{\sqrt{\left( -1 + \mu + x\left( t \right) \right)^{2} + \left( y\left( t \right) \right)^{2} + \left( z\left( t \right) \right)^{2}}} \right)^{3} \mu \right) z\left( t \right)
\end{aligned}$

## Examples

State transition dynamics are particularly valuable for CR3BP models.
Recall that the state transition matrix is simply the local linearization
of a spacecraft within CR3BP dynamics. Let's look at the Jacobian (another
word for "local linearization") below, evaluated at some random state.

```@repl main
f = CR3BPFunction(; jac=true)
let x = randn(6), p = rand((0.0, 0.5)), t = 0
    f.jac(x, p, t)
end
```

The Jacobian will always have this form (zeros in the top-left,
the identity matrix in the top-right, a dense matrix in the 
bottom-left, and the same sparse "-2, 2" matrix in the bottom-right).
We can include the state transition dynamics in our model with 
`stm=true`, initialize the state transition matrix states to the 
identity matrix, and propagate our spacecraft for one periodic orbit:
the result is known as the Monodromy Matrix! The Monodromy Matrix
provides stability characteristics for the entire periodic orbit.

```@repl main
model = CR3BP(; stm=true)
```

Note that periodic orbits are not easy to find within CR3BP dynamics.
Various algorithms have been developed to analytically approximate, 
and numerically refine, periodic CR3BP orbits. Some of those 
algorithms have already been implemented in Julia! See 
[`OrbitalTrajectories`](https://github.com/dpad/OrbitalTrajectories.jl)
and [`GeneralAstrodynamics`](https://github.com/cadojo/GeneralAstrodynamics.jl).