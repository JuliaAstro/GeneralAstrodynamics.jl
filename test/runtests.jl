using Pkg, Test

const GROUP = get(ENV, "GROUP", "ALL")

# 
# Copied from https://github.com/SciML/Optimization.jl/blob/master/test/runtests.jl
#
# Copyright (c) 2022 Julia Computing

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

function activate_subpkg_env(subpkg)
    subpkg_path = joinpath(dirname(@__DIR__), "lib", subpkg)
    Pkg.activate(subpkg_path)
    return Pkg.instantiate()
end

if uppercase(GROUP) in ("ALL", "GENERALASTRODYNAMICS.JL")
    for package in readdir(joinpath(@__DIR__, "..", "lib"))
        activate_subpkg_env(package)
        Pkg.test()
    end
else
    activate_subpkg_env(GROUP)
    Pkg.test()
end