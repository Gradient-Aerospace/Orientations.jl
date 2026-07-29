export RollPitchYaw, RPY, RPY_F64
export RollPitchYawDeg, RPYDeg, RPYDeg_F64

"""
Represents an orientation as a roll, pitch, and yaw (rad) from a reference. Yaw is the first
rotation, followed by pitch, followed by roll. That is, this describes a frame oriented
`roll` around the x axis of a frame that's oriented `pitch` around the y axis of a frame
that's oriented `yaw` around the z axis of the reference frame.
"""
@kwdef struct RollPitchYaw{T} <: AbstractOrientation{T}
    roll::T
    pitch::T
    yaw::T
end
const RPY = RollPitchYaw
const RPY_F64 = RollPitchYaw{Float64}

"""
Represents an orientation as a roll, pitch, and yaw (deg) from a reference. Yaw is the first
rotation, followed by pitch, followed by roll. That is, this describes a frame oriented
`roll` around the x axis of a frame that's oriented `pitch` around the y axis of a frame
that's oriented `yaw` around the z axis of the reference frame.

This type can be converted to RollPitchYaw via `deg2rad`.
"""
@kwdef struct RollPitchYawDeg{T} <: AbstractOrientation{T}
    roll::T
    pitch::T
    yaw::T
end
const RPYDeg = RollPitchYawDeg
const RPYDeg_F64 = RollPitchYawDeg{Float64}

##############################
# Conversion To/From Degrees #
##############################

Base.deg2rad(rpy_deg::RPYDeg{T}) where {T} = RPY{T}(
    deg2rad(rpy_deg.roll),
    deg2rad(rpy_deg.pitch),
    deg2rad(rpy_deg.yaw),
)

Base.rad2deg(rpy::RPY{T}) where {T} = RPYDeg{T}(
    rad2deg(rpy.roll),
    rad2deg(rpy.pitch),
    rad2deg(rpy.yaw),
)

Base.convert(::Type{RPY}, rpy_deg::RPYDeg) = deg2rad(rpy_deg)
Base.convert(::Type{RPY{T}}, rpy_deg::RPYDeg) where {T} = convert(RPY{T}, deg2rad(rpy_deg))
Base.convert(::Type{RPYDeg}, rpy::RPY) = rad2deg(rpy)
Base.convert(::Type{RPYDeg{T}}, rpy::RPY) where {T} = convert(RPYDeg{T}, rad2deg(rpy))

################
# Constructors #
################

Base.one(::Type{RPY}) = one(RPY_F64)
Base.one(::Type{RPY{T}}) where {T} = RPY{T}(zero(T), zero(T), zero(T))

Base.one(::Type{RPYDeg}) = one(RPYDeg_F64)
Base.one(::Type{RPYDeg{T}}) where {T} = RPYDeg{T}(zero(T), zero(T), zero(T))

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{RPYDeg{T}}) where {T}
    return rad2deg(rand(rng, RPY{T}))
end

##############
# Operations #
##############

# We use the abstract fallback for all operations. RPYs aren't really meant for operations.

##########################
# Indexing and Iteration #
##########################

# RollPitchYaw and RollPitchYawDeg have the same indexing behavior.
const RollPitchYawTypes = Union{RollPitchYaw, RollPitchYawDeg}

Base.length(::RollPitchYawTypes) = 3
Base.size(::RollPitchYawTypes) = (3,)
Base.axes(::RollPitchYawTypes) = (Base.OneTo(3),)
function Base.axes(rpy::RollPitchYawTypes, k)
    if k != 1
        throw(BoundsError(rpy, k))
    end
    return Base.OneTo(3)
end
Base.keys(::RollPitchYawTypes) = Base.OneTo(3)
Base.firstindex(::RollPitchYawTypes) = 1
Base.lastindex(::RollPitchYawTypes) = 3

# Allow users to iterate over the angles, e.g. for splatting.
function Base.iterate(rpy::RollPitchYawTypes, state = 1)
    state == 1 && return (rpy.roll, state + 1)
    state == 2 && return (rpy.pitch, state + 1)
    state == 3 && return (rpy.yaw, state + 1)
    return nothing
end

# Provide linear indexing behavior.
function Base.getindex(rpy::RollPitchYawTypes, k)
    k == 1 && return rpy.roll
    k == 2 && return rpy.pitch
    k == 3 && return rpy.yaw
    throw(BoundsError(rpy, k))
end

# Let's tell the user they can't do this.
function Base.setindex!(rpy::RollPitchYawTypes, value, k)
    error("$(nameof(typeof(rpy))) is immutable and cannot support setindex!.")
end
