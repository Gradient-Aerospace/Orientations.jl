export AxisAngle, AA, AA_F64, AADeg, AADeg_F64
export aax, aay, aaz

"""
Represents an orientation that is rotated a given angle (rad) about the given axis from a
reference. The fields are `axis` and `angle`.

The constructor stores the axis exactly as provided; it does not normalize the axis. Use
`LinearAlgebra.normalize(aa)` to return an equivalent `AxisAngle` with a unit axis.
"""
@kwdef struct AxisAngle{T} <: AbstractOrientation{T}
    axis::SVector{3, T}
    angle::T
    AxisAngle{T}(axis::SVector{3, T}, angle::T) where {T} = new{T}(axis, angle)
end
const AA = AxisAngle
const AA_F64 = AxisAngle{Float64}
AxisAngle(axis::SVector{3, T}, angle::T) where {T} = AxisAngle{T}(axis, angle)

"""
Represents an orientation that is rotated a given angle (deg) about the given axis from a
reference. The fields are `axis` and `angle`.

The constructor stores the axis exactly as provided; it does not normalize the axis. Use
`LinearAlgebra.normalize(aa_deg)` to return an equivalent `AxisAngleDeg` with a unit axis.

This type can be converted to `AxisAngle` via the `deg2rad` function.
"""
@kwdef struct AxisAngleDeg{T} <: AbstractOrientation{T}
    axis::SVector{3, T}
    angle::T
    AxisAngleDeg{T}(axis::SVector{3, T}, angle::T) where {T} = new{T}(axis, angle)
end
const AADeg = AxisAngleDeg
const AADeg_F64 = AxisAngleDeg{Float64}
AxisAngleDeg(axis::SVector{3, T}, angle::T) where {T} = AxisAngleDeg{T}(axis, angle)

Base.rad2deg(aa::AxisAngle) = AADeg(aa.axis, rad2deg(aa.angle))
Base.deg2rad(aa_deg::AxisAngleDeg) = AA(aa_deg.axis, deg2rad(aa_deg.angle))
Base.convert(::Type{AxisAngle}, aa_deg::AxisAngleDeg) = deg2rad(aa_deg)
Base.convert(::Type{AxisAngle{T}}, aa_deg::AxisAngleDeg) where {T} = convert(AxisAngle{T}, deg2rad(aa_deg))
Base.convert(::Type{AxisAngleDeg}, aa::AxisAngle) = rad2deg(aa)
Base.convert(::Type{AxisAngleDeg{T}}, aa::AxisAngle) where {T} = convert(AxisAngleDeg{T}, rad2deg(aa))

################
# Constructors #
################

"Constructs an AxisAngle orientation from a 3-element axis vector and an angle."
function AxisAngle(axis::AbstractVector, angle)
    @assert length(axis) == 3 "Cannot construct AxisAngle from a $(length(axis))-element axis."
    T = promote_type(eltype(axis), typeof(angle))
    return AxisAngle{T}(SVector{3, T}(axis), convert(T, angle))
end

"Constructs an AxisAngle orientation from a 3-element axis vector and an angle."
function AxisAngle{T}(axis::AbstractVector, angle) where {T}
    @assert length(axis) == 3 "Cannot construct AxisAngle from a $(length(axis))-element axis."
    return AxisAngle{T}(SVector{3, T}(axis), convert(T, angle))
end

"Constructs an AxisAngleDeg orientation from a 3-element axis vector and an angle."
function AxisAngleDeg(axis::AbstractVector, angle)
    @assert length(axis) == 3 "Cannot construct AxisAngleDeg from a $(length(axis))-element axis."
    T = promote_type(eltype(axis), typeof(angle))
    return AxisAngleDeg{T}(SVector{3, T}(axis), convert(T, angle))
end

"Constructs an AxisAngleDeg orientation from a 3-element axis vector and an angle."
function AxisAngleDeg{T}(axis::AbstractVector, angle) where {T}
    @assert length(axis) == 3 "Cannot construct AxisAngleDeg from a $(length(axis))-element axis."
    return AxisAngleDeg{T}(SVector{3, T}(axis), convert(T, angle))
end

"Constructs an AxisAngle orientation that is rotated by the given angle about the x axis."
aax(θ::T) where {T} = AA{T}(SA[one(T), zero(T), zero(T)], θ)

"Constructs an AxisAngle orientation that is rotated by the given angle about the y axis."
aay(θ::T) where {T} = AA{T}(SA[zero(T), one(T), zero(T)], θ)

"Constructs an AxisAngle orientation that is rotated by the given angle about the z axis."
aaz(θ::T) where {T} = AA{T}(SA[zero(T), zero(T), one(T)], θ)

Base.one(::Type{AA}) = one(AA_F64)
Base.one(::Type{AA{T}}) where {T} = AA{T}(SA[one(T), zero(T), zero(T)], zero(T))

Base.one(::Type{AADeg}) = one(AADeg_F64)
Base.one(::Type{AADeg{T}}) where {T} = AADeg{T}(SA[one(T), zero(T), zero(T)], zero(T))

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{AADeg{T}}) where {T}
    return rad2deg(rand(rng, AA{T}))
end

##############
# Operations #
##############

# We could swap the angle or axis. We choose the axis.
Base.inv(aa::AA) = typeof(aa)(aa.axis, -aa.angle)

"""
Returns an equivalent AxisAngle with a unit axis.
"""
function LinearAlgebra.normalize(aa::AxisAngle)
    axis_norm = norm(aa.axis)
    @assert axis_norm > zero(axis_norm) "Cannot normalize an AxisAngle with a zero axis."
    axis = aa.axis ./ axis_norm
    T = promote_type(eltype(axis), typeof(aa.angle))
    return AxisAngle{T}(convert(SVector{3, T}, axis), convert(T, aa.angle))
end

"""
Returns an equivalent AxisAngleDeg with a unit axis.
"""
function LinearAlgebra.normalize(aa_deg::AxisAngleDeg)
    axis_norm = norm(aa_deg.axis)
    @assert axis_norm > zero(axis_norm) "Cannot normalize an AxisAngleDeg with a zero axis."
    axis = aa_deg.axis ./ axis_norm
    T = promote_type(eltype(axis), typeof(aa_deg.angle))
    return AxisAngleDeg{T}(convert(SVector{3, T}, axis), convert(T, aa_deg.angle))
end

# Going to ERPs is actually the best way to do this calculation.
compose(a::AA, b::AA) = erp2aa(compose(aa2erp(a), aa2erp(b)))

# We can reframe by going to DCM, but actually it's slightly faster without
# going to DCM.
function reframe(aa::AA, v)
    aa = normalize(aa)
    T = eltype(aa.axis)
    s, c = sincos(aa.angle)
    r = aa.axis
    return c .* v + ((one(T) - c) * (r ⋅ v)) .* r - s .* cross(r, v)
end

# Like composition, going to ERPs is the best way to do this.
difference(a::AA, b::AA) = erp2aa(difference(aa2erp(a), aa2erp(b)))

# Getting the distance is trivial.
function distance(aa::AA)
    angle = mod(aa.angle, 2π) # Now angle is in [0, 2pi).
    if angle > π
        return 2π - angle
    else
        return angle
    end
end

function smallest(aa::AA)
    angle = mod(aa.angle, 2π) # Now angle is in [0, 2pi).
    if angle > π
        return AA(-aa.axis, 2π - angle)
    else
        return AA(aa.axis, angle)
    end
end

other(aa::AA) = AA(-aa.axis, -aa.angle)

# This is also best handled with ERP interpolation.
interpolate(a::AA, b::AA, f) = erp2aa(interpolate(aa2erp(a), aa2erp(b), f))

#############
# Iteration #
#############

# Helpful functions for iterating over AxisAngle.
Base.length(::AA) = 2
Base.eltype(::AA{T}) where {T} = Union{T, SVector{3, T}}
Base.size(::AA) = (2,)

# Allow a user to iterate over the elements of AA, e.g. for splatting.
function Base.iterate(aa::AA, state = 1)
    state == 1 && return (aa.axis, state + 1)
    state == 2 && return (aa.angle, state + 1)
    return nothing
end

# Provide linear indexing behavior.
function Base.getindex(aa::AA, k)
    k == 1 && return aa.axis
    k == 2 && return aa.angle
    throw(BoundsError(aa, k))
end
Base.firstindex(aa::AA) = 1
Base.lastindex(aa::AA) = 2

# Let's tell the user they can't do this.
function Base.setindex!(aa::AA, value, k)
    error("AxisAngle is immutable and cannot support setindex!.")
end

#################
# Miscellaneous #
#################

function Base.show(io::IO, aa::AxisAngle)
    print(io, "AxisAngle(axis = $(aa.axis), angle = $(aa.angle))")
end

########################################
# Conversions to Non-Orientation Types #
########################################

# What would be a normal type to convert to? A tuple?
