# Explicit compressed sensing, solving BPDN with explicit convex optimization
__precompile__()

module ExplicitCS

using Images
using FFTW
using Convex
using SCS

include("Utils.jl")

export convex_sense

# Bases Ψ
stdbasis(n) = Matrix(1.0I, n, n)
dctbasis(n) = dct(stdbasis(n), 1)

coeffs(basis, imgdata) = (imgdata' * basis)'

# Θ = ΦΨ, where Ψ is the DCT-II basis.
Θdct(Φ, s) = Φ * dct(s)

# Basis persuit denoising (BPDN)
#
# Finds (sparse) s that minimizes ∥Θs - y∥² + λ∥s∥₁, where
# y is a vector of measurements,
# Θ is as above, and
# λ is the regularization parameter that determines the balance between
# data misfit and the penalty from the 1-norm of the coefficients s (sparsity).
#
function sparse_coeffs(Θ, y, λ)
    s = Variable(size(Θ)[2])
    problem = minimize(sumsquares(Θ*s - y) + λ*norm(s, 1))
    solve!(problem, SCSSolver())
    s.value, problem.status
end

# Simulate compressed sensing
function convex_sense(img, filename, M, λ, basis=dctbasis)
    N = length(img)
    Φ = randmask(M, N)
    y = Φ * vec(float.(channelview(img)))
    Ψ = basis(N)
    Θ = Φ * Ψ
    (s, status) = sparse_coeffs(Θ, y, λ)
    if status == :Optimal
        x = Ψ * s
        print("Successful compressed sensing.")
        out = reshape(Gray.(renormalize(x)), size(img))
        save(filename, out)
        out
    else
        error("Could not find sparse coefficients.")
    end
end

end

