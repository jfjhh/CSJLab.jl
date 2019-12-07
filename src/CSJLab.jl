__precompile__()

"""
# Compressed sensing JLab project code

This package provides three ways to form "single-pixel" images.

- `GhostImage`: Computational ghost imaging.
- `ExplicitCS`: Explicit compressive sensing (BPDN) by convex optimization.
- `ImplicitCS`: Implicit compressive sensing (BPDN) via the GenSPGL package.
"""
module CSJLab

export CSImage, loadimage, imvec, vecim, randmask, measure, SNR, renormalize, binvalue, contrast

include("Utils.jl")
include("GhostImage.jl")
include("ExplicitCS.jl")
include("ImplicitCS.jl")

end # module

