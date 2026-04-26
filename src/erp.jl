export EulerRodriguesParameters, ERP, ERP_F64
export erpx, erpy, erpz

"""
Represents an orientation using Euler-Rodrigues Parameters, as described by Shuster in "A
Survey of Attitude Representations".

This is the same type as a quaternion or rotation where the scalar is last. It follows the
JPL convention for composition rules and conversion to other types.
"""
@kwdef struct EulerRodriguesParameters{T} <: AbstractOrientation{T}
    x::T
    y::T
    z::T
    s::T
end
const ERP = EulerRodriguesParameters
const ERP_F64 = EulerRodriguesParameters{Float64}

################
# Constructors #
################

"Constructs a EulerRodriguesParameters from a 4-element vector."
function EulerRodriguesParameters(v::AbstractVector)
    @assert length(v) == 4 "Cannot construct EulerRodriguesParameters from a $(length(v))-element vector."
    return EulerRodriguesParameters(v...)
end

"Converts a 4-element vector to a EulerRodriguesParameters."
function Base.convert(type::Type{<:EulerRodriguesParameters}, v::AbstractVector)
    @assert length(v) == 4 "Cannot convert a $(length(v))-element vector to EulerRodriguesParameters."
    return type(v...)
end

"Creates a EulerRodriguesParameters representing no rotation."
Base.zero(::ERP{T}) where {T} = ERP{T}(zero(T), zero(T), zero(T), one(T))

"Creates a EulerRodriguesParameters representing no rotation."
Base.zero(::Type{ERP{T}}) where {T} = ERP{T}(zero(T), zero(T), zero(T), one(T))

"Create a random EulerRodriguesParameters following a uniform distribution on SO(3)."
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{ERP{T}}) where {T}
    return normalize(ERP(randn(rng, T), randn(rng, T), randn(rng, T), randn(rng, T)))
end

"Returns the EulerRodriguesParameters corresponding to a rotation of `θ` about x."
function erpx(θ)
    s, c = sincos(θ/2)
    n = zero(s)
    return EulerRodriguesParameters(s, n, n, c)
end

"Returns the EulerRodriguesParameters corresponding to a rotation of `θ` about y."
function erpy(θ)
    s, c = sincos(θ/2)
    n = zero(s)
    return EulerRodriguesParameters(n, s, n, c)
end

"Returns the EulerRodriguesParameters corresponding to a rotation of `θ` about z."
function erpz(θ)
    s, c = sincos(θ/2)
    n = zero(s)
    return EulerRodriguesParameters(n, n, s, c)
end

##############
# Operations #
##############

"Returns the opposite orientation of the input. (Shuster eq. 177)"
function Base.inv(erp::EulerRodriguesParameters)
    return typeof(erp)(-erp.x, -erp.y, -erp.z, erp.s)
end

"""
Returns the EulerRodriguesParameters corresponding to "the other way around the rotation
axis".
"""
function other(erp::EulerRodriguesParameters)
    return typeof(erp)(-erp.x, -erp.y, -erp.z, -erp.s)
end

"""
Returns the smaller of the two rotations these EulerRodriguesParameters represent (positive
scalar) -- "shortest way around the rotation axis."
"""
function smallest(erp::EulerRodriguesParameters)
    if erp.s >= 0
        return erp
    else
        return other(erp)
    end
end

"""
    reframe(erp_BA::EulerRodriguesParameters, v_A::NTuple)

If `erp_BA` is the orientation of frame B wrt frame A and `v_A` is a vector expressed in
frame A, this returns the vector expressed in frame B. (Shuster eq. 158)
"""
function reframe(e::EulerRodriguesParameters, v::NTuple{3, T}) where {T}
    v1, v2, v3 = v
    x = e.x
    y = e.y
    z = e.z
    s = e.s
    x2 = x^2
    y2 = y^2
    z2 = z^2
    s2 = s^2
    return (
            (x2 - y2 - z2 + s2) * v1
        + 2 * (x * y + s * z) * v2
        + 2 * (x * z - s * y) * v3,
            2 * (y * x - s * z) * v1
        + (-x2 + y2 - z2 + s2) * v2
        + 2 * (y * z + s * x) * v3,
            2 * (z * x + s * y) * v1
        + 2 * (z * y - s * x) * v2
        + (-x2 - y2 + z2 + s2) * v3,
    )
end

"""
    reframe(erp_BA::EulerRodriguesParameters, v_A::Vector)

If `erp_BA` is the orientation of frame B wrt frame A and `v_A` is a vector expressed in
frame A, this returns the vector expressed in frame B. (Shuster eq. 158)
"""
function reframe(e::EulerRodriguesParameters, v::Vector{T}) where {T}
    @assert length(v) == 3 "reframe expects a 3-element vector."
    return T[reframe(e, (v[1], v[2], v[3]))...]
end

"""
    reframe(erp_BA::EulerRodriguesParameters, v_A::AbstractVector)

If `erp_BA` is the orientation of frame B wrt frame A and `v_A` is a vector expressed in
frame A, this returns the vector expressed in frame B. (Shuster eq. 158)
"""
function reframe(e::EulerRodriguesParameters, v::T) where {T <: AbstractVector}
    @assert length(v) == 3 "reframe expects 3-element vectors."
    return T(reframe(e, (v[1], v[2], v[3]))...)
end

"""
Composes two ERP s.t. if the first argument is the rotation of frame C wrt frame B, and
the second is the rotation of B wrt A, the result is the rotation of frame C wrt
frame A. (Shuster eq. 173)
"""
function compose(a::EulerRodriguesParameters, b::EulerRodriguesParameters)
    return normalize(
        EulerRodriguesParameters(
            a.s * b.x + a.z * b.y - a.y * b.z + a.x * b.s,
            -a.z * b.x + a.s * b.y + a.x * b.z + a.y * b.s,
            a.y * b.x - a.x * b.y + a.s * b.z + a.z * b.s,
            -a.x * b.x - a.y * b.y - a.z * b.z + a.s * b.s,
        )
    )
end

"""
If `a` is the rotation of frame A wrt C and `b` is the rotation of frame
B wrt C, then this function returns the rotation of A wrt B.
"""
function difference(a::EulerRodriguesParameters, b::EulerRodriguesParameters)
    return compose(a, inv(b))
end

"""
If `a` is the rotation of frame A wrt C and `b` is the rotation of frame
B wrt C, then this function returns the angle by which A is rotated from C (smallest angle).
"""
function distance(a::EulerRodriguesParameters, b::EulerRodriguesParameters)
    return distance(difference(a, b))
end

"""
Returns the "distance" (non-negative angle, the "short way around") of the ERP wrt the
reference orientation.
"""
function distance(a::EulerRodriguesParameters{T}) where {T}
    return 2 * atan(sqrt(a.x^2 + a.y^2 + a.z^2), abs(a.s))
end

"""
    rate(erp::EulerRodriguesParameters, ω; k = 1/2)

Returns the rate of change of the ERP of B wrt A (`erp`) given the rotation rate of B
wrt A (`ω`), expressed in either B or A (doesn't matter which). (Shuster eq. 306)

This function includes a correction factor (`k`) to push the norm towards unity. This is
useful for numerical integration and should have no affect if the norm is already unity. In
general, choosing `k` such that `k * Dt < 1`, where `Dt` is the overall timestep
used for integration, provides good stability.

(Note that we're using the ERP type as a container, and the thing being contained is the
rate of change of the ERP, not valid ERP themselves. E.g., the output shouldn't have a
unit norm.)
"""
function rate(erp::EulerRodriguesParameters{T}, ω, k = one(T)/2) where {T}
    hω1 = ω[1] / 2
    hω2 = ω[2] / 2
    hω3 = ω[3] / 2
    correction = k * (one(T) - squared_norm(erp))
    return EulerRodriguesParameters(
         erp.s * hω1 - erp.z * hω2 + erp.y * hω3 + correction * erp.x,
         erp.z * hω1 + erp.s * hω2 - erp.x * hω3 + correction * erp.y,
        -erp.y * hω1 + erp.x * hω2 + erp.s * hω3 + correction * erp.z,
        -erp.x * hω1 - erp.y * hω2 - erp.z * hω3 + correction * erp.s,
    )
end

"""
    interpolate(a, b, f; shortest_path = true)

Spherical linear interpolation between two EulerRodriguesParameters.

The parameter `f` can range from 0 to 1 to specify where along the sphere to interpolate.
When f = 0, the result is `a` and when f = 1, the result is `b`.

Setting `shortest_path` to false will allow the interpolation to take the "long way around".

References
- https://en.wikipedia.org/wiki/Slerp
- https://blog.magnum.graphics/backstage/the-unnecessarily-short-ways-to-do-a-quaternion-slerp/
"""
function interpolate(ep_start::EulerRodriguesParameters{T}, ep_end::EulerRodriguesParameters{T}, f; shortest_path::Bool = true) where {T}

    if f < 0 || f > 1
        throw(DomainError(f, "f was $f but must be in the range [0, 1]."))
    end
    if iszero(f)
        return ep_start
    end
    if isone(f)
        return ep_end
    end

    # This is the intuitive implementation, but what's below is faster.
    #
    #   aa = erp2aa(smallest(difference(ep_end, ep_start)))
    #   return normalize(compose(aa2erp(AxisAngle(aa.axis, aa.angle * f)), ep_start))
    #

    d = dot(ep_start, ep_end)
    if d <= -one(T) && !shortest_path
        return ep_start
    end

    if shortest_path
        # Flip `ep_start` if there is a shorter path.
        ep_startp = d >= 0 ? ep_start : other(ep_start)
    else
        ep_startp = ep_start
    end

    θ = acos(clamp(abs(d), zero(T), one(T)))
    if θ <= sqrt(eps(float(T)))
        return normalize((1 - f) * ep_startp + f * ep_end)
    end

    # If -1 < cos(θ) < 1, then sin(θ) should be > 0.
    c1 = sin((1 - f) * θ) / sin(θ)
    c2 = sin(f * θ) / sin(θ)
    return normalize(c1 * ep_startp + c2 * ep_end)

end

#############
# Iteration #
#############

# Helpful functions for iterating over EulerRodriguesParameters.
Base.length(::EulerRodriguesParameters) = 4
Base.eltype(::EulerRodriguesParameters{T}) where {T} = T
Base.size(::EulerRodriguesParameters) = (4,)

# Allow a user to iterate over the elements of ERP, e.g. for splatting.
function Base.iterate(erp::EulerRodriguesParameters, state = 1)
    state == 1 && return (erp.x, state + 1)
    state == 2 && return (erp.y, state + 1)
    state == 3 && return (erp.z, state + 1)
    state == 4 && return (erp.s, state + 1)
    return nothing
end

# Provide linear indexing behavior.
function Base.getindex(erp::EulerRodriguesParameters, k)
    k == 1 && return erp.x
    k == 2 && return erp.y
    k == 3 && return erp.z
    k == 4 && return erp.s
    throw(BoundsError(erp, k))
end
Base.firstindex(erp::EulerRodriguesParameters) = 1
Base.lastindex(erp::EulerRodriguesParameters) = 4

# Let's tell the user they can't do this.
function Base.setindex!(erp::EulerRodriguesParameters, value, k)
    error("EulerRodriguesParameters are immutable and cannot support setindex!.")
end

##############
# Arithmetic #
##############

"""
Returns the squared norm of the ERP.
"""
function squared_norm(erp::EulerRodriguesParameters)
    return erp.x * erp.x + erp.y * erp.y + erp.z * erp.z + erp.s * erp.s
end

"""
Returns the 2-norm of the ERP.
"""
function LinearAlgebra.norm(erp::EulerRodriguesParameters)
    return sqrt(squared_norm(erp))
end

"""
Returns the EulerRodriguesParameters divided by its 2-norm.
"""
function LinearAlgebra.normalize(erp::EulerRodriguesParameters)
    m = norm(erp)
    @assert m > 0 "The EulerRodriguesParameters had a non-positive magnitude."
    return typeof(erp)(erp.x/m, erp.y/m, erp.z/m, erp.s/m)
end

"Adds the elements of two EulerRodriguesParameters."
function Base.:+(a::T, b::T) where {T <: EulerRodriguesParameters}
    return T(a.x + b.x, a.y + b.y, a.z + b.z, a.s + b.s)
end

"Subtracts the elements of second EulerRodriguesParameters from the first."
function Base.:-(a::T, b::T) where {T <: EulerRodriguesParameters}
    return T(a.x - b.x, a.y - b.y, a.z - b.z, a.s - b.s)
end

"Multiplies the elements of a EulerRodriguesParameters by a scalar (first argument)."
function Base.:*(s::ST, a::T) where {ST <: Number, T <: EulerRodriguesParameters}
    return T(s * a.x, s * a.y, s * a.z, s * a.s)
end

"Multiplies the elements of a EulerRodriguesParameters by a scalar (second argument)."
function Base.:*(a::T, s::ST) where {T <: EulerRodriguesParameters, ST <: Number}
    return T(a.x * s, a.y * s, a.z * s, a.s * s)
end


#################
# Miscellaneous #
#################

"""
A convenient method for printing EP that highlights the vector and scalar parts.
(We don't show i,j,k; those are for quaternions and not part of the EP notation.)
"""
function Base.show(io::IO, erp::EulerRodriguesParameters)
    print(io, "EulerRodriguesParameters(x = ", erp.x, ", y = ", erp.y, ", z = ", erp.z, ", s = ", erp.s, ")")
end

########################################
# Conversions to Non-Orientation Types #
########################################

# Get a Tuple of the internal values.
function Base.convert(type::Type{NTuple{4, T}}, erp::EulerRodriguesParameters) where {T}
    return type((erp.x, erp.y, erp.z, erp.s))
end

# Get a Vector of the internal values.
function Base.convert(::Type{Vector{T}}, erp::EulerRodriguesParameters) where {T}
    return T[erp.x, erp.y, erp.z, erp.s]
end

# Get a Dict of the internal values.
function Base.convert(type::Type{<:AbstractDict{KT, VT}}, erp::EulerRodriguesParameters) where {KT <: AbstractString, VT}
    return type("x" => erp.x, "y" => erp.y, "z" => erp.z, "s" => erp.s)
end
