# Attitude Dynamics

_Quaternion kinematics, and dynamics!_

```@setup main
using AstrodynamicalModels
using ModelingToolkit
using Latexify
Latexify.auto_display(true)
```

## Overview

The Attitude model assumes a spacecraft with some orientation described by a
**scalar-last** quaternion, and body rates which are small enough such that they
appear constant for small numerical integration tolerance values.

!!! danger 
    You should normalize the quaternion vector at each time step using a
    `ManifoldCallback` or `DiscreteCallback` when simulating this model!
    Without normalizing, the solution will drift such that the quaternion
    state vector is no longer a unit quaternion. The dynamics in this 
    model *assume* a unit quaternion norm!

$\begin{aligned}
    \dot{q} &= \frac{1}{2} \begin{bmatrix}
        0 && \omega_3 && -\omega_2 && \omega_1 \\
        -\omega_3 && 0 && \omega_1 && \omega_2 \\
        \omega_2 && -\omega_1 && 0 && \omega_3 \\
        -\omega_1 && -\omega_2 && -\omega_3 && 0
    \end{bmatrix} q \\
    \dot{\omega} &= -J^{-1} (\omega\times) J \omega + J^{-1} L + u \\
\end{aligned}$

## Examples

```@repl main
model = AttitudeSystem()
```

Let's compute the Jacobian for these dynamics.

```@repl main
J = calculate_jacobian(AttitudeSystem())
```

Finally, let's construct a Julia function which implements these dynamics!

```@repl main
f = AttitudeFunction()
let u = randn(7), p = randn(15), t = 0
    f(u, p, t)
end
```
