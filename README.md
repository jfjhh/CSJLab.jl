# CSJLab

This is a small Julia package that provides *computational ghost imaging* and
*compressed sensing* (basis pursuit denoising) functionality for constructing
images using a "single-pixel camera" proof-of-concept. Three submodules are
provided: `GhostImaging`, `ImplicitCS`, and `ExplicitCS`.

## The setup

We have set up a photodiode with a square centimeter of sensitive area to read
the light intensity from a computer monitor configured with *LabView* to output
random black and white pixel patterns in a square grid. We place an opaque
object like a strip of paper or a transparency sheet between the displayed
pattern and the photodiode and record the pattern used and the intensity
measured. From these data, the code in this package can be used to reconstruct
the transmission function (an image) of the object.

## Submodules

The computational ghost imaging (`GhostImage` module) is simple enough, and just
uses standard linear algebra.

The basis pursuit denoising is implemented using implicit operators via
[GenSPGL.jl](https://github.com/slimgroup/GenSPGL.jl) package and is made
available in the `ImplicitCS` module.

Explicit solution capability by convex optimization is also provided by the
`ExplicitCS` module, which uses
[Convex.jl](https://github.com/JuliaOpt/Convex.jl), but this quickly becomes
intractable as the size of the image grows.

