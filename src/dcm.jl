export DirectionCosineMatrix, DCM, DCM_F64
export Rx, Ry, Rz

"""
Represents an orientation as a direction cosine matrix, with a single field called `matrix`.
"""
@kwdef struct DirectionCosineMatrix{T} <: AbstractOrientation{T} # Not an AbstractMatrix!
    matrix::SMatrix{3, 3, T, 9}
end
const DCM = DirectionCosineMatrix
const DCM_F64 = DirectionCosineMatrix{Float64}

################
# Constructors #
################

"""
Creates a direction cosine matrix such that if frame B is rotated θ about x from frame A,
then a vector expressed in frame B can be determined from a vector in frame A as:

```
v_B = Rx(θ) * v_A
```
"""
function Rx(θ)
    s, c = sincos(θ)
    u = one(typeof(s))
    z = zero(typeof(s))
    return DCM(
        @SMatrix [
            u  z  z;
            z  c  s;
            z -s  c;
        ]
    )
end

"""
Creates a direction cosine matrix such that if frame B is rotated θ about y from frame A,
then a vector expressed in frame B can be determined from a vector in frame A as:

```
v_B = Ry(θ) * v_A
```
"""
function Ry(θ)
    s, c = sincos(θ)
    u = one(typeof(s))
    z = zero(typeof(s))
    return DCM(
        @SMatrix [
            c  z -s;
            z  u  z;
            s  z  c;
        ]
    )
end

"""
Creates a direction cosine matrix such that if frame B is rotated θ about z from frame A,
then a vector expressed in frame B can be determined from a vector in frame A as:

```
v_B = Rz(θ) * v_A
```
"""
function Rz(θ)
    s, c = sincos(θ)
    u = one(typeof(s))
    z = zero(typeof(s))
    return DCM(
        @SMatrix [
            c  s  z;
           -s  c  z;
            z  z  u;
        ]
    )
end

# Note that zero means "no rotation". You might be tempted to try `one` because you're
# thinking of an identity matrix, but that's not the right concept here! A "one" orientation
# would be an orientation displaced by "one unit" (and about what??) from a reference
# orientation.
Base.zero(dcm::DCM{T}) where {T} = zero(typeof(dcm))
Base.zero(::Type{DCM}) = zero(DCM_F64)
Base.zero(::Type{DCM{T}}) where {T} = DCM(one(SMatrix{3, 3, T, 9})) # <- A bit ironic...

"Converts a matrix to a DirectionCosineMatrix."
function Base.convert(type::Type{<:DCM}, m::AbstractMatrix)
    @assert size(m) == (3,3) "Cannot convert a $(size(m)) matrix to a DirectionCosineMatrix."
    return type(convert(SMatrix, m))
end

"Converts a StaticMatrix to a DirectionCosineMatrix."
function Base.convert(type::Type{<:DCM}, m::SMatrix)
    @assert size(m) == (3,3) "Cannot convert a $(size(m)) matrix to a DirectionCosineMatrix."
    return type(m)
end

##############
# Operations #
##############

Base.inv(dcm::DCM) = DCM(dcm.matrix')
compose(a::DCM, b::DCM) = DCM(a.matrix * b.matrix)
reframe(dcm::DCM, v) = dcm.matrix * v
difference(a::DCM, b::DCM) = DCM(a.matrix * b.matrix')
distance(dcm::DCM) = distance(dcm2erp(dcm))
interpolate(a::DCM, b::DCM, f) = erp2dcm(interpolate(dcm2erp(a), dcm2erp(b), f))

###################
# Matrix Behavior #
###################

# Let this behave like a regular matix, insofar as that makes sense for a DCM.
@inline Base.:*(a::DCM, b::DCM) = DCM(a.matrix * b.matrix)
@inline Base.:*(dcm::DCM, rhs) = dcm.matrix * rhs
@inline Base.:*(lhs, dcm::DCM) = lhs * dcm.matrix
@inline Base.getindex(dcm::DCM, args...) = getindex(dcm.matrix, args...)
@inline LinearAlgebra.tr(dcm::DCM) = tr(dcm.matrix)
# Note: We don't support broadcasting or scalar multiplication because that probably means
# the resulting type isn't a DCM any more. If we supported `rate`, that would make sense
# though.

#################
# Miscellaneous #
#################

# The default `show` implementation is fine here; we don't need anything special.

########################################
# Conversions to Non-Orientation Types #
########################################

Base.convert(t::Type{<:Matrix}, dcm::DCM) = convert(t, dcm.matrix)
Base.convert(::Type{<:SMatrix}, dcm::DCM) = dcm.matrix
