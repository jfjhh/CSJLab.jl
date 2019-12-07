# Utility functions for compressed sensing
__precompile__()

using LinearAlgebra
using Images
using DelimitedFiles

"Compressed sensing image data"
struct CSImage{T}
    "The measurements of light intensity"
    measurements :: Array{T,1}
    "The patterns that the light was transmitted through, each as a row"
    patterns     :: Array{T,2}
end

"""
Read image data CSV file with `filename` as `(measurements, patterns)`.
"""
function loadimage(filename)
    d = readdlm(filename, ',')
    CSImage(d[:,1], d[:,2:end])
end

"Turn an image into a vector suitable for processing."
imvec(img, conversion=float) = vec(channelview(conversion.(Gray.(img))))

"Turn a processed vector `vec` into a grayscale image with given `size`."
vecim(vec, size = (Int.(sqrt(length(vec))), :)) = Gray.(reshape(vec, size))

"""
Create a random mask array of zero or one entries.

### Examples

```
julia> randmask(2, 3)
2×3 Array{Float64,2}:
 1.0  0.0  1.0
 1.0  0.0  0.0
```
"""
randmask(dims...) = (sign.(randn(dims...)) .+ 1.0) ./ 2.0

"Simulates measurements of `img` through `pattern` or through random masks."
function measure(img, m, patterns = randmask(m, length(img)))
    CSImage(patterns * img, patterns)
end

"""
Calculate the signal-to-noise ratio of `signal` (\$x\$) to `actual` (\$y\$)
using

\$\\textrm{SNR} = -20 \\log\\left(\\frac{\\|x - y\\|}{\\|y\\|}\\right).\$
"""
SNR(signal, actual) = -20.0 * log10(norm(signal - actual) / norm(actual))

"Rescale `data` to the range \$[0, 1]\$, then multiply by `scale` if given."
function renormalize(data, scale=oneunit(eltype(data)))
    (a, b) = extrema(data)
    if a == b
        fill(scale / 2, size(data))
    else
        @. (data - a) * (scale / (b - a))
    end
end

"""
Creates a function that bins the values of an image according to `values`.

The partition `values` is of the form `[a_0, a_1, ..., a_n]`, where
\$0 = a_0 < a_1 < \\cdots < a_n = 1.\$
By default, the value for a given bin is the midpoint of the values at the
endpoints of the bin, but this may be changed by specifying `bins`.

### Examples

```
julia> binvalue([0.0, 0.2, 0.8, 1.0]).([0.01, 0.3, 0.9])
3-element Array{Float64,1}:
 0.1
 0.5
 0.9

julia> binvalue([0.0, 0.5, 1.0], [0.0, 0.8]).([0.2, 0.6])
2-element Array{Float64,1}:
 0.0
 0.8
```
"""
function binvalue(values, bins = (values[2:end] + values[1:end-1]) ./ 2.0)
    v -> if v ≤ values[1]
        bins[1]
    else
        bins[(1:length(values))[v .≤ values][1] - 1]
    end
end

"Change contrast and brightness of `img` by a clamped linear transformation."
contrast(a, b, img) = Gray.(clamp.(a .* img .+ b, 0.0, 1.0))

