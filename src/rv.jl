export RotationVector, RV, RV_F64

"""
Represents an orientation using a rotation vector -- essentially the angle of rotation times
the axis of rotation between the two frames. It has a single field, `vector`.
"""
@kwdef struct RotationVector{T} <: AbstractOrientation{T}
    vector::SVector{3, T}
    RotationVector{T}(vector::SVector{3, T}) where {T} = new{T}(vector)
end
const RV = RotationVector
const RV_F64 = RotationVector{Float64}
RotationVector(vector::SVector{3, T}) where {T} = RotationVector{T}(vector)

################
# Constructors #
################

"Constructs a RotationVector from a 3-element vector."
function RotationVector(v::AbstractVector)
    @assert length(v) == 3 "Cannot construct RotationVector from a $(length(v))-element vector."
    T = eltype(v)
    return RotationVector{T}(SVector{3, T}(v))
end

"Constructs a RotationVector from a 3-element vector."
function RotationVector{T}(v::AbstractVector) where {T}
    @assert length(v) == 3 "Cannot construct RotationVector from a $(length(v))-element vector."
    return RotationVector{T}(SVector{3, T}(v))
end

Base.convert(type::Type{<:RV}, v::AbstractVector) = type(v)

Base.one(::Type{RV}) = one(RV_F64)
Base.one(::Type{RV{T}}) where {T} = RV(zero(SVector{3, T}))

# TODO: We could do the deg2rad/rad2deg thing here too.

##############
# Operations #
##############

Base.inv(rv::RV) = RV(-rv.vector)
compose(a::RV, b::RV) = erp2rv(compose(rv2erp(a), rv2erp(b)))
reframe(rv::RV, v) = reframe(rv2aa(rv), v)
difference(a::RV, b::RV) = erp2rv(difference(rv2erp(a), rv2erp(b)))

function distance(rv::RV)
    angle = norm(rv.vector)
    angle = mod(angle, 2π)
    if angle > π
        return 2π - angle
    else
        return angle
    end
end

interpolate(a::RV, b::RV, f) = erp2rv(interpolate(rv2erp(a), rv2erp(b), f))

#############
# Iteration #
#############

# TODO? Do we want to be able to splat the elements of an RV?

#################
# Miscellaneous #
#################

# The default `show` method is fine.

########################################
# Conversions to Non-Orientation Types #
########################################

# Get a Vector of the internal values.
function Base.convert(type::Type{<:AbstractVector}, rv::RotationVector)
    return convert(type, rv.vector)
end
