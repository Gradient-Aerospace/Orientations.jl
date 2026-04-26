"A package for representing and operating on orientation-related types."
module GradientOrientations

using LinearAlgebra
using Random
using StaticArrays

"""
All subtypes of AbstractOrientation are expected to support:

* reframe(orientation, vector)
* compose(a, b)
* difference(a, b)
* distance(orientation)
* distance(a, b)
* interpolate(a, b, fraction)
* Base.one
* Base.inv
* Random.rand

To enable a type to support numerical integration, implement `rate` for the type.

Where it makes sense, types may implement:

* other(orientation)
* smallest(orientation)

It's helpful to write conversions for each abstract type, especially to and from
EulerRodriguesParameters. That provides a way to convert to all other types. Where more
direct conversions exist, those are helpful to implement as well.
"""
abstract type AbstractOrientation{T} end

# Necessary things
export reframe, compose, difference, distance, interpolate, identity_orientation

# May not exist for all orientation types
export other, smallest, rate

# Users can `convert` from one type to another, but this lets them specify their intention
# more clearly and is easier to write.
export aa2dcm, aa2erp, aa2rpy, aa2rv
export dcm2aa, dcm2erp, dcm2rpy, dcm2rv
export erp2aa, erp2dcm, erp2rpy, erp2rv
export rv2aa, rv2dcm, rv2erp, rv2rpy
export rpy2aa, rpy2dcm, rpy2erp, rpy2rv

# Composition operator
export ⊗

include("utilities.jl")

# The primary types
include("erp.jl")
include("aa.jl")
include("dcm.jl")
include("rpy.jl")
include("rv.jl")
include("conversions.jl")

# Implement some fallback methods for AbstractOrientation that just convert to ERP, do the
# operation, and then convert back.

"""
    reframe(B_wrt_A::AbstractOrientation, v_A)

If `B_wrt_A` is the orientation of frame B wrt frame A and `v_A` is a vector expressed in
frame A, this returns the vector expressed in frame B. (Shuster eq. 158)
"""
function reframe(a::AbstractOrientation, v)
    return reframe(convert(ERP, a), v)
end

"""
Composes two orientations s.t. if the first argument is the orientation of frame C wrt frame
B, and the second is the orientation of B wrt A, the result is the orientation of frame C
wrt frame A. (Shuster eq. 173)
"""
function compose(a::T, b::T) where {T <: AbstractOrientation}
    return convert(T, compose(convert(ERP, a), convert(ERP, b)))
end

"""
If `a` is the orientation of frame A wrt C and `b` is the orientation of frame B wrt C, then
this function returns the orientation of A wrt B.
"""
function difference(a::T, b::T) where {T <: AbstractOrientation}
    return convert(T, difference(convert(ERP, a), convert(ERP, b)))
end

"""
Returns the angle that the given orientation is rotated wrt its reference (shortest way
around).
"""
function distance(a::AbstractOrientation{T}) where {T}
    return distance(convert(ERP{T}, a))
end

"""
If `a` is the orientation of frame A wrt C and `b` is the orientation of frame
B wrt C, then this function returns the angle by which A is rotated from C (shortest way
around).
"""
function distance(a::T, b::T) where {T <: AbstractOrientation}
    return distance(convert(ERP, a), convert(ERP, b))
end

"""
Compares two orientations by angular distance in radians.
"""
function Base.isapprox(a::AbstractOrientation, b::AbstractOrientation; atol = eps(4π), rtol = zero(atol), kwargs...)
    return distance(convert(ERP, a), convert(ERP, b)) <= max(atol, rtol * max(distance(a), distance(b)))
end

"""
    interpolate(a, b, t; shortest_path = true)

Spherical linear interpolation between two orientations.

The parameter `f` can range from 0 to 1 to specify where along the sphere to interpolate.
When f = 0, the result is `a` and when f = 1, the result is `b`.

Setting `shortest_path` to false will allow the interpolation to take the "long way around".

References
- https://en.wikipedia.org/wiki/Slerp
- https://blog.magnum.graphics/backstage/the-unnecessarily-short-ways-to-do-a-quaternion-slerp/
"""
function interpolate(a::T, b::T, f; kwargs...) where {T <: AbstractOrientation}
    return convert(T, interpolate(convert(ERP, a), convert(ERP, b), f; kwargs...))
end

"Returns the element type used by this orientation."
Base.eltype(::AbstractOrientation{T}) where {T} = T
Base.eltype(::Type{OT}) where {T, OT <: AbstractOrientation{T}} = T

"""
Returns an orientation that is not rotated from its reference, using the same type as the input.
"""
Base.one(x::AbstractOrientation) = one(typeof(x))

"""
Returns an orientation that is not rotated from its reference, using the given type.
"""
Base.one(type::Type{OT}) where {T, OT <: AbstractOrientation{T}} = convert(OT, one(ERP{T}))

"""
    identity_orientation(x::AbstractOrientation)
    identity_orientation(::Type{<:AbstractOrientation})

Returns an orientation that is not rotated from its reference.
"""
identity_orientation(x::AbstractOrientation) = one(x)
identity_orientation(type::Type{<:AbstractOrientation}) = one(type)

"""
Inverts the orientation such that if the input is the orientation of B wrt A, this returns
the orientation of A wrt B.
"""
Base.inv(x::AbstractOrientation) = convert(typeof(x), inv(convert(ERP, x)))

"""
Returns a uniform random orientation with the given type.
"""
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{OT}) where {OT <: AbstractOrientation{T}} where {T}
    return convert(OT, rand(rng, ERP{T}))
end

"Composition operator, with the same interface as `compose`."
⊗(a, b) = compose(a, b)

# Allow a user to compose multiple orientations.
compose(a, b, c, args...) = compose(compose(a, b), c, args...)

end # module GradientOrientations
