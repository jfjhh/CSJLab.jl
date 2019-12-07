__precompile__()

"""
## Computational Ghost Imaging
"""
module GhostImage

export ghostimage, simulate_ghostimage

using LinearAlgebra
using Images

include("Utils.jl")

function ghostimage(cs)
    N   = length(cs.measurements)
    Σs  = sum(cs.measurements)
    ΣsI = sum(cs.patterns .* cs.measurements, dims=1)
    ΣI  = sum(cs.patterns, dims=1)

    renormalize((ΣsI / N) - (Σs * ΣI / (N*N)))
end

simulate_ghostimage(img, m) = ghostimage(measure(img, m))

end

