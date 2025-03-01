var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"AstrodynamicalSolvers","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers","text":"All exported names.","category":"page"},{"location":"reference/","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers","text":"Modules = [\n    AstrodynamicalSolvers,\n    AstrodynamicalSolvers.Propagation,\n]\nOrder = [:module, :type, :function, :constant]","category":"page"},{"location":"reference/#AstrodynamicalSolvers.AstrodynamicalSolvers","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers.AstrodynamicalSolvers","text":"Provides astrodynamical solvers, including Lyapunov and halo orbit correctors.\n\nExtended help\n\nLicense\n\nMIT License\n\nCopyright (c) 2023 Joseph D Carpinelli\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n\nExports\n\nCR3BSolvers\nPropagation\nconvergent_manifold\ndivergent_manifold\nhalo\nlyapunov\nmonodromy\npropagate\npropagate!\n\nImports\n\nAstrodynamicalSolvers.CR3BSolvers\nAstrodynamicalSolvers.Propagation\nBase\nCore\nDocStringExtensions\nReexport\n\n\n\n\n\n","category":"module"},{"location":"reference/#AstrodynamicalSolvers.Propagation","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers.Propagation","text":"Wrappers around SciML differential equation solvers for fast and convenient  orbit propagation.\n\nExtended Help\n\nExports\n\nconvergent_manifold\ndivergent_manifold\nmonodromy\npropagate\npropagate!\n\nImports\n\nAstrodynamicalCalculations\nAstrodynamicalModels\nBase\nCore\nDocStringExtensions\nModelingToolkit\nOrdinaryDiffEqVerner\nSciMLBase\nStaticArrays\n\n\n\n\n\n","category":"module"},{"location":"reference/#AstrodynamicalSolvers.Propagation.convergent_manifold-Tuple{Any, Any, Any}","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers.Propagation.convergent_manifold","text":"convergent_manifold(u, μ, Δt; eps, trajectories, kwargs...)\n\n\nReturn a vector of orbits along the manifold which converges to the provided  halo orbit.\n\n\n\n\n\n","category":"method"},{"location":"reference/#AstrodynamicalSolvers.Propagation.divergent_manifold-Tuple{Any, Any, Any}","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers.Propagation.divergent_manifold","text":"divergent_manifold(u, μ, Δt; eps, trajectories, kwargs...)\n\n\nReturn a vector of orbits along the manifold which diverges from the provided  halo orbit.\n\n\n\n\n\n","category":"method"},{"location":"reference/#AstrodynamicalSolvers.Propagation.monodromy-Tuple{AbstractVector, Any, Any, Function}","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers.Propagation.monodromy","text":"monodromy(\n    u,\n    μ,\n    T,\n    f;\n    algorithm,\n    reltol,\n    abstol,\n    save_everystep,\n    kwargs...\n)\n\n\nSolve for the monodromy matrix of the periodic orbit.\n\n\n\n\n\n","category":"method"},{"location":"reference/#AstrodynamicalSolvers.Propagation.monodromy-Tuple{AstrodynamicalModels.AstrodynamicalOrbit, Any}","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers.Propagation.monodromy","text":"monodromy(orbit, Δt; algorithm, reltol, abstol, kwargs...)\n\n\nCompute the monodromy matrix for any periodic orbit.\n\n\n\n\n\n","category":"method"},{"location":"reference/#AstrodynamicalSolvers.Propagation.propagate!-Tuple{AstrodynamicalModels.AstrodynamicalOrbit, Any}","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers.Propagation.propagate!","text":"propagate!(\n    orbit,\n    Δt;\n    stm,\n    algorithm,\n    reltol,\n    abstol,\n    kwargs...\n)\n\n\nNumerically integrate the orbit forward (or backward) in time, modifying the  state vector in-place within the AstrodynamicalOrbit instance.\n\n\n\n\n\n","category":"method"},{"location":"reference/#AstrodynamicalSolvers.Propagation.propagate-Tuple{AstrodynamicalModels.AstrodynamicalOrbit, Any}","page":"AstrodynamicalSolvers","title":"AstrodynamicalSolvers.Propagation.propagate","text":"propagate(\n    orbit,\n    Δt;\n    stm,\n    algorithm,\n    reltol,\n    abstol,\n    kwargs...\n)\n\n\nNumerically integrate the orbit forward (or backward) in time, and return a new  AstrodynamicalOrbit instance with identical parameters to the provided orbit.\n\n\n\n\n\n","category":"method"},{"location":"#AstrodynamicalSolvers.jl","page":"Getting Started","title":"AstrodynamicalSolvers.jl","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Common solvers within orbital mechanics and astrodynamics.","category":"page"},{"location":"#Installation","page":"Getting Started","title":"Installation","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"pkg> add AstrodynamicalSolvers","category":"page"},{"location":"#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"This package currently provides periodic orbit, and manifold computations within  Circular Restricted Three Body Problem dynamics.","category":"page"},{"location":"#Periodic-Orbits","page":"Getting Started","title":"Periodic Orbits","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"This package contains differential correctors, and helpful wrapper functions, for  finding periodic orbits within Circular Restricted Three Body Problem dynamics.","category":"page"},{"location":"#Plots.jl","page":"Getting Started","title":"Plots.jl","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"using AstrodynamicalSolvers\nusing AstrodynamicalModels\nusing OrdinaryDiffEq\nusing Plots\n\nμ = 0.012150584395829193\n\nplanar = let\n    ic = halo(μ, 1) # lyapunov (planar) orbit\n    u = Vector(CartesianState(ic))\n    problem = ODEProblem(CR3BFunction(), u, (0, ic.Δt), (μ,))\n    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14)\n    plot(solution, idxs=(:x,:y,:z), title = \"Lyapunov Orbit\", label=:none, size=(1600,900), dpi=400, aspect_ratio=1)\nend\n\nextraplanar = let\n    ic = halo(μ, 2; amplitude=0.01) # halo (non-planar) orbit\n    u = Vector(CartesianState(ic))\n    problem = ODEProblem(CR3BFunction(), u, (0, ic.Δt), (μ,))\n    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14)\n    plot(solution, idxs=(:x,:y,:z), title = \"Halo Orbit\", label=:none, size=(1600,900), dpi=400, aspect_ratio=1)\nend\n\nplot(planar, extraplanar, layout=(1,2))","category":"page"},{"location":"#Makie.jl","page":"Getting Started","title":"Makie.jl","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"using AstrodynamicalSolvers\nusing AstrodynamicalModels\nusing OrdinaryDiffEq\nusing CairoMakie\n\nμ = 0.012150584395829193\n\nsol_planar = let\n    ic = halo(μ, 1) # lyapunov (planar) orbit\n    u = Vector(CartesianState(ic))\n    problem = ODEProblem(CR3BFunction(), u, (0, ic.Δt), (μ,))\n    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14)\nend\n\nsol_extraplanar = let\n    ic = halo(μ, 2; amplitude=0.01) # halo (non-planar) orbit\n    u = Vector(CartesianState(ic))\n    problem = ODEProblem(CR3BFunction(), u, (0, ic.Δt), (μ,))\n    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14)\nend\n\nfig = Figure(size=(800, 400); fontsize=11)\n\nax_kwargs_common = (; aspect=:equal, azimuth=-π/3)\n\nax_left = Axis3(fig[1, 1];\n    title = \"Lyapunov Orbit\",\n    limits = (0.78, 0.90, -0.09, 0.09, -0.02, 1.04),\n    ax_kwargs_common...,\n)\nax_right = Axis3(fig[1, 2];\n    title = \"Halo Orbit\",\n    limits = (1.05, 1.26, -0.1, 0.1, -0.02, 0.01),\n    protrusions = (30, 100, 0, 0),\n    ax_kwargs_common...,\n)\n\nplot!(ax_left, sol_planar; idxs=(:x, :y, :z))\nplot!(ax_right, sol_extraplanar; idxs=(:x, :y, :z))\n\nfig","category":"page"},{"location":"#Manifold-Computations","page":"Getting Started","title":"Manifold Computations","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Manifold computations, provided by AstrodynamicalCalculations.jl, can perturb  halo orbits onto their unstable or stable manifolds.","category":"page"},{"location":"#Plots.jl-2","page":"Getting Started","title":"Plots.jl","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"using AstrodynamicalSolvers\nusing AstrodynamicalCalculations\nusing AstrodynamicalModels\nusing OrdinaryDiffEq\nusing LinearAlgebra\nusing Plots\n\nμ = 0.012150584395829193\n\nunstable = let\n    ic = halo(μ, 1; amplitude=0.005)\n\n    u = CartesianState(ic)\n    Φ = monodromy(u, μ, ic.Δt, CR3BFunction(stm=true))\n\n    ics = let\n        problem = ODEProblem(CR3BFunction(stm=true), vcat(u, vec(I(6))), (0, ic.Δt), (μ,))\n        solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=(ic.Δt / 10))\n\n        solution.u\n    end\n\n    perturbations = [\n        diverge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=-1e-7)\n        for ic in ics\n    ]\n\n    problem = EnsembleProblem(\n        ODEProblem(CR3BFunction(), u, (0.0, 2 * ic.Δt), (μ,)),\n        prob_func=(prob, i, repeat) -> remake(prob; u0=perturbations[i]),\n    )\n\n    solution = solve(problem, Vern9(), trajectories=length(perturbations), reltol=1e-14, abstol=1e-14)\nend\n\nstable = let\n    ic = halo(μ, 2; amplitude=0.005)\n\n    u = CartesianState(ic)\n    Φ = monodromy(u, μ, ic.Δt, CR3BFunction(stm=true))\n\n    ics = let\n        problem = ODEProblem(CR3BFunction(stm=true), vcat(u, vec(I(6))), (0, ic.Δt), (μ,))\n        solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=(ic.Δt / 10))\n\n        solution.u\n    end\n    \n    perturbations = [\n        converge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=1e-7)\n        for ic in ics\n    ]\n\n    problem = EnsembleProblem(\n        ODEProblem(CR3BFunction(), u, (0.0, -2.1 * ic.Δt), (μ,)),\n        prob_func=(prob, i, repeat) -> remake(prob; u0=perturbations[i]),\n    )\n\n    solution = solve(problem, Vern9(), trajectories=length(perturbations), reltol=1e-14, abstol=1e-14)\nend\n\nfigure = plot(; \n    aspect_ratio = 1.0,\n    background = :transparent,\n    grid = true,\n    title = \"Unstable and Stable Invariant Manifolds\",\n)\n\nplot!(figure, unstable, idxs=(:x, :y), aspect_ratio=1, label=:none, palette=:blues)\nplot!(figure, stable, idxs=(:x, :y), aspect_ratio=1, label=:none, palette=:blues)\nscatter!(figure, [1-μ], [0], label=\"Moon\", xlabel=\"X (Earth-Moon Distance)\", ylabel=\"Y (Earth-Moon Distance)\", marker=:x, color=:black, markersize=10,)\n\nfigure # hide","category":"page"},{"location":"#Makie.jl-2","page":"Getting Started","title":"Makie.jl","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"using AstrodynamicalSolvers\nusing AstrodynamicalCalculations\nusing AstrodynamicalModels\nusing OrdinaryDiffEq\nusing LinearAlgebra\nusing CairoMakie\n\nμ = 0.012150584395829193\n\nunstable = let\n    ic = halo(μ, 1; amplitude=0.005)\n\n    u = CartesianState(ic)\n    Φ = monodromy(u, μ, ic.Δt, CR3BFunction(stm=true))\n\n    ics = let\n        problem = ODEProblem(CR3BFunction(stm=true), vcat(u, vec(I(6))), (0, ic.Δt), (μ,))\n        solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=(ic.Δt / 10))\n\n        solution.u\n    end\n\n    perturbations = [\n        diverge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=-1e-7)\n        for ic in ics\n    ]\n\n    problem = EnsembleProblem(\n        ODEProblem(CR3BFunction(), u, (0.0, 2 * ic.Δt), (μ,)),\n        prob_func=(prob, i, repeat) -> remake(prob; u0=perturbations[i]),\n    )\n\n    solution = solve(problem, Vern9(), trajectories=length(perturbations), reltol=1e-14, abstol=1e-14)\nend\n\nstable = let\n    ic = halo(μ, 2; amplitude=0.005)\n\n    u = CartesianState(ic)\n    Φ = monodromy(u, μ, ic.Δt, CR3BFunction(stm=true))\n\n    ics = let\n        problem = ODEProblem(CR3BFunction(stm=true), vcat(u, vec(I(6))), (0, ic.Δt), (μ,))\n        solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=(ic.Δt / 10))\n\n        solution.u\n    end\n\n    perturbations = [\n        converge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=1e-7)\n        for ic in ics\n    ]\n\n    problem = EnsembleProblem(\n        ODEProblem(CR3BFunction(), u, (0.0, -2.1 * ic.Δt), (μ,)),\n        prob_func=(prob, i, repeat) -> remake(prob; u0=perturbations[i]),\n    )\n\n    solution = solve(problem, Vern9(), trajectories=length(perturbations), reltol=1e-14, abstol=1e-14)\nend\n\nfig = Figure(size=(800, 400), fontsize=20)\n\nax = Axis(fig[1, 1];\n    xreversed = true,\n    xticks = LinearTicks(5),\n    yticks = LinearTicks(5),\n    aspect = DataAspect(),\n    xlabel = \"X (Earth-Moon Distance)\",\n    ylabel = \"Y (Earth-Moon Distance)\",\n    title = \"Unstable and Stable Invariant Manifolds\",\n    titlesize = 24,\n)\n\nidxs = (:x, :y)\n\n# TODO: replace this manual workaround when\n# https://github.com/SciML/SciMLBase.jl/issues/697#issuecomment-2135801331\n# is addressed\nfor (traj, color) in zip(unstable, resample_cmap(:blues, length(unstable)))\n    plot!(ax, traj; idxs, color)\nend\n\nfor (traj, color) in zip(stable, resample_cmap(:blues, length(stable)))\n    plot!(ax, traj; idxs, color)\nend\n\nscatter!(ax, [1-μ], [0]; marker='⨯', color=:black, markersize=50, label=\"Moon\")\n\nfig","category":"page"},{"location":"cr3bp/#Circular-Restricted-Three-Body-Solvers","page":"CR3BSolvers","title":"Circular Restricted Three Body Solvers","text":"","category":"section"},{"location":"cr3bp/","page":"CR3BSolvers","title":"CR3BSolvers","text":"All three-body solvers!","category":"page"},{"location":"cr3bp/","page":"CR3BSolvers","title":"CR3BSolvers","text":"Modules = [\n    AstrodynamicalSolvers.CR3BSolvers,\n]\nOrder = [:module, :type, :function, :constant]","category":"page"},{"location":"cr3bp/#AstrodynamicalSolvers.CR3BSolvers","page":"CR3BSolvers","title":"AstrodynamicalSolvers.CR3BSolvers","text":"Solvers specific to the Circular Restricted Three Body Problem.\n\nExtended Help\n\nExports\n\nhalo\nlyapunov\n\nImports\n\nAstrodynamicalCalculations\nAstrodynamicalModels\nBase\nCore\nDocStringExtensions\nLinearAlgebra\nModelingToolkit\nOrdinaryDiffEqVerner\nStaticArrays\n\n\n\n\n\n","category":"module"},{"location":"cr3bp/#AstrodynamicalSolvers.CR3BSolvers.extraplanar_differential-Tuple{AbstractVector, Any}","page":"CR3BSolvers","title":"AstrodynamicalSolvers.CR3BSolvers.extraplanar_differential","text":"extraplanar_differential(state, μ)\n\n\nwarning: CR3BP Dynamics\nThis computation is valid for Circular Restricted Three Body Problem dynamics.\n\nGiven a full state vector for CR3BP dynamics, including vertically concatenated columns of the state transition matrix, return the differential correction term for a periodic orbit.\n\n\n\n\n\n","category":"method"},{"location":"cr3bp/#AstrodynamicalSolvers.CR3BSolvers.halo-NTuple{5, Any}","page":"CR3BSolvers","title":"AstrodynamicalSolvers.CR3BSolvers.halo","text":"halo(x, z, ẏ, μ, T; reltol, abstol, maxiters)\n\n\nwarning: CR3BP Dynamics\nThis computation is valid for Circular Restricted Three Body Problem dynamics.\n\nIterate on an initial guess for halo orbit conditions.\n\n\n\n\n\n","category":"method"},{"location":"cr3bp/#AstrodynamicalSolvers.CR3BSolvers.halo-Tuple{Any, Int64}","page":"CR3BSolvers","title":"AstrodynamicalSolvers.CR3BSolvers.halo","text":"halo(μ, lagrange; amplitude, phase, hemisphere, kwargs...)\n\n\nwarning: CR3BP Dynamics\nThis computation is valid for Circular Restricted Three Body Problem dynamics.\n\nGiven a nondimensional mass parameter μ, and orbit characteristics, construct  an initial guess using Richardson's analytical solution, and iterate on that guess using a differential corrector. \n\n\n\n\n\n","category":"method"},{"location":"cr3bp/#AstrodynamicalSolvers.CR3BSolvers.lyapunov-NTuple{4, Any}","page":"CR3BSolvers","title":"AstrodynamicalSolvers.CR3BSolvers.lyapunov","text":"lyapunov(x, ẏ, μ, T; reltol, abstol, maxiters)\n\n\nwarning: CR3BP Dynamics\nThis computation is valid for Circular Restricted Three Body Problem dynamics.\n\nIterate on an initial guess for Lyapunov orbit conditions.\n\n\n\n\n\n","category":"method"},{"location":"cr3bp/#AstrodynamicalSolvers.CR3BSolvers.planar_differential-Tuple{AbstractVector, Any}","page":"CR3BSolvers","title":"AstrodynamicalSolvers.CR3BSolvers.planar_differential","text":"planar_differential(state, μ)\n\n\nwarning: CR3BP Dynamics\nThis computation is valid for Circular Restricted Three Body Problem dynamics.\n\nGiven a full state vector for CR3BP dynamics, including vertically concatenated columns of the state transition matrix, return the differential correction term for a planar periodic orbit.\n\n\n\n\n\n","category":"method"}]
}
