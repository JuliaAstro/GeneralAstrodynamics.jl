# Restricted Two-body Dynamics

_Also known as R2BP dynamics!_

```@setup main
using AstrodynamicalModels
using ModelingToolkit
using Latexify
Latexify.auto_display(true)
```

## Overview

The Restricted Two-body Problem (R2BP) assumes a massless spacecraft which moves
due to the gravity of **one** celestial body: one star, or one planet, or one
moon, or one asteroid. The equations of motion for R2BP dynamics are shown
below.

$\begin{align*}
\frac{dx(t)}{dt} =& ẋ\left( t \right) \\
\frac{dy(t)}{dt} =& ẏ\left( t \right) \\
\frac{dz(t)}{dt} =& ż\left( t \right) \\
\frac{dẋ(t)}{dt} =& \frac{ - \mu x\left( t \right)}{\left( \sqrt{x^2\left(t\right) + y^2\left(t\right) + z^2\left(t\right)} \right)^{3}} \\
\frac{dẏ(t)}{dt} =& \frac{ - \mu y\left( t \right)}{\left( \sqrt{x^2\left(t\right) + y^2\left(t\right) + z^2\left(t\right)} \right)^{3}} \\
\frac{dż(t)}{dt} =& \frac{ - \mu z\left( t \right)}{\left( \sqrt{x^2\left(t\right) + y^2\left(t\right) + z^2\left(t\right)} \right)^{3}}
\end{align*}$

## Examples

```@repl main
model = R2BP()
```

Every model also offers _optional_ state transition matrix dynamics. Use
`stm=true` to append the state transition matrix dynamics to your model's
equations of motion. State transition dynamics can also be thought of the
model's _local linearization_.

!!! note The state transition dynamics for `R2BP` systems are not _nearly_ as
useful as the state transition dynamics within [`CR3BP`](CR3BP.md) models.
Within CR3BP dynamics, a spacecraft's local linearization offers stability
characteristics for periodic orbits, and provides stable and unstable directions
(in state-space) for invariant manifolds about periodic orbits and Lagrange
points.

```@repl main
model = R2BP(; stm=true)
```

Let's compute the Jacobian for these dynamics.

```@repl main
J = calculate_jacobian(R2BP())
```

Finally, let's construct a Julia function which implements these dynamics!

```@repl main
f = R2BPFunction()
let u = randn(6), p = [3e6], t = 0
    f(u, p, t)
end
```
