"""
# Compressed sensing JLab project code

This package provides three ways to form "single-pixel" images.

- `GhostImage`: Computational ghost imaging.
- `ExplicitCS`: Explicit compressive sensing (BPDN) by convex optimization.
- `ImplicitCS`: Implicit compressive sensing (BPDN) via the GenSPGL package.
"""
module CSJLab

using LinearAlgebra
using LinearAlgebra.BLAS
using DelimitedFiles
using Images.ImageCore
using FFTW
using JOLI
using GenSPGL

include("Utils.jl")
include("GhostImage.jl")
include("ImplicitCS.jl")

end # module

