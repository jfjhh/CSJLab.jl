__precompile__()

"""
Implicit compressed sensing via the GenSPGL package, which allows an implicit
DCT-II transform that significantly speeds things up.
"""
module ImplicitCS

using LinearAlgebra
using LinearAlgebra.BLAS
using Images
using FFTW
using JOLI
using GenSPGL

include("Utils.jl")

export csimage, simulate_csimage

const dσ = 1e-1
const doptions = (verbosity=0,)

function csimage(cs; σ = dσ, options = doptions)
    BLAS.set_num_threads(16)
    s = spgl1(joMatrix(cs.patterns) * joDCT(size(cs.patterns)[2]),
              cs.measurements,
              tau     = 0.0,
              sigma   = σ,
              options = spgOptions(;stepMin=10.0*eps(), iterations=5000,
                                   options...))[1]
    renormalize(dct(s))
end

"Simulate compressed sensing of `img` with `m` measurements."
function simulate_csimage(img, m, σ = dσ, options = doptions)
    spgl_sense(measure(imvec(img), m), σ, options)
end

end

