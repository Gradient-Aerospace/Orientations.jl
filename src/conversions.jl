############
# AA to... #
############

"Returns the DirectionCosineMatrix for the given AxisAngle."
function aa2dcm(aa::AA)
    aa = normalize(aa)
    T = eltype(aa.axis)
    s, c = sincos(aa.angle)
    r = aa.axis
    R = diagm(SVector{3, T}(c, c, c)) + (one(T) - c) .* (r * r') - s .* crs3(r)
    return DCM{T}(R)
end
Base.convert(::Type{DCM}, aa::AA) = aa2dcm(aa) # Type is not specified on the LHS.
Base.convert(::Type{DCM{T}}, aa::AA{T}) where {T} = aa2dcm(aa) # Type is not specified on the LHS.

"Returns the EulerRodriguesParameters for the given AxisAngle."
function aa2erp(aa::AA)
    aa = normalize(aa)
    T = eltype(aa.axis)
    r = aa.axis
    s, c = sincos(aa.angle/2)
    return ERP{T}(s * r[1], s * r[2], s * r[3], c)
end
Base.convert(::Type{ERP}, aa::AA) = aa2erp(aa) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, aa::AA{T}) where {T} = aa2erp(aa) # Type is not specified on the LHS.
Base.convert(::Type{ERP}, aa_deg::AADeg) = aa2erp(deg2rad(aa_deg)) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, aa_deg::AADeg{T}) where {T} = aa2erp(deg2rad(aa_deg))

"Returns the RotationVector for the given AxisAngle."
function aa2rv(aa::AA)
    aa = normalize(aa)
    T = eltype(aa.axis)
    return RV{T}(aa.angle .* aa.axis)
end
Base.convert(::Type{RV}, aa::AA) = aa2rv(aa) # Type is not specified on the LHS.
Base.convert(::Type{RV{T}}, aa::AA{T}) where {T} = aa2rv(aa) # Type is not specified on the LHS.

"Returns the RollPitchYaw for the given AxisAngle."
function aa2rpy(aa::AA)
    return erp2rpy(aa2erp(aa))
end
Base.convert(::Type{RPY}, aa::AA) = aa2rpy(aa) # Type is not specified on the LHS.
Base.convert(::Type{RPY{T}}, aa::AA{T}) where {T} = aa2rpy(aa) # Type is not specified on the LHS.

#############
# DCM to... #
#############

"Returns the AxisAngle for the given DirectionCosineMatrix."
function dcm2aa(dcm::DCM{T}) where {T}

    R = dcm.matrix
    vx = R[2, 3] - R[3, 2]
    vy = R[3, 1] - R[1, 3]
    vz = R[1, 2] - R[2, 1]

    # Calculate angle from atan for best numerical stability.
    s = sqrt(vx^2 + vy^2 + vz^2)
    c = tr(R) - one(T)
    angle = atan(s, c)

    # Identity rotation (or sufficiently close to it).
    if iszero(angle) || angle <= sqrt(eps(T))
        return AA{T}(SA[one(T), zero(T), zero(T)], zero(T))
    end

    # Near pi, the skew-symmetric terms become too small for robust axis extraction.
    sin_angle = s/T(2)
    if sin_angle <= sqrt(eps(T))

        ax = sqrt(max(zero(T), (R[1, 1] + one(T)) / T(2)))
        ay = sqrt(max(zero(T), (R[2, 2] + one(T)) / T(2)))
        az = sqrt(max(zero(T), (R[3, 3] + one(T)) / T(2)))

        if ax >= ay && ax >= az
            ay = copysign(ay, R[1, 2] + R[2, 1])
            az = copysign(az, R[1, 3] + R[3, 1])
        elseif ay >= ax && ay >= az
            ax = copysign(ax, R[1, 2] + R[2, 1])
            az = copysign(az, R[2, 3] + R[3, 2])
        else
            ax = copysign(ax, R[1, 3] + R[3, 1])
            ay = copysign(ay, R[2, 3] + R[3, 2])
        end

        axis = normalize(SVector{3, T}(ax, ay, az))

        return AA{T}(axis, angle)

    end

    axis = SVector{3, T}(vx/s, vy/s, vz/s)

    return AA{T}(axis, angle)

end
Base.convert(::Type{AA}, dcm::DCM) = dcm2aa(dcm) # Type is not specified on the LHS.
Base.convert(::Type{AA{T}}, dcm::DCM{T}) where {T} = dcm2aa(dcm) # Type is not specified on the LHS.

"Returns the EulerRodriguesParameters for the given DirectionCosineMatrix."
function dcm2erp(dcm::DCM{T}) where {T}

    R = dcm.matrix

    # Split the conversion so as to divide by the largest possible number.
    if R[1,1] + R[2,2] + R[3,3] >= 0
        η4 = 0.5 * √(1. + R[1,1] + R[2,2] + R[3,3])
        α  = 0.25 / η4 # Necessarily safe from above
        return EulerRodriguesParameters{T}(
            α * (R[2,3] - R[3,2]),
            α * (R[3,1] - R[1,3]),
            α * (R[1,2] - R[2,1]),
            η4,
        )
    elseif R[1,1] - R[2,2] - R[3,3] >= 0
        η1 = 0.5 * √(1. + R[1,1] - R[2,2] - R[3,3])
        α  = 0.25 / η1 # Necessarily safe from above
        return EulerRodriguesParameters{T}(
            η1,
            α * (R[1,2] + R[2,1]),
            α * (R[3,1] + R[1,3]),
            α * (R[2,3] - R[3,2]),
        )
    elseif - R[1,1] + R[2,2] - R[3,3] >= 0
        η2 = 0.5 * √(1. - R[1,1] + R[2,2] - R[3,3])
        α  = 0.25 / η2 # Necessarily safe from above
        return EulerRodriguesParameters{T}(
            α * (R[1,2] + R[2,1]),
            η2,
            α * (R[3,2] + R[2,3]),
            α * (R[3,1] - R[1,3]),
        )
    else
        η3 = 0.5 * √(1. - R[1,1] - R[2,2] + R[3,3])
        α  = 0.25 / η3 # Safe if R is a DCM
        return EulerRodriguesParameters{T}(
            α * (R[1,3] + R[3,1]),
            α * (R[3,2] + R[2,3]),
            η3,
            α * (R[1,2] - R[2,1]),
        )
    end

end
Base.convert(::Type{ERP}, dcm::DCM) = dcm2erp(dcm) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, dcm::DCM{T}) where {T} = dcm2erp(dcm) # Type is not specified on the LHS.

"Returns the RotationVector for the given DirectionCosineMatrix."
function dcm2rv(dcm::DCM)
    return aa2rv(dcm2aa(dcm))
end
Base.convert(::Type{RV}, dcm::DCM) = dcm2rv(dcm) # Type is not specified on the LHS.
Base.convert(::Type{RV{T}}, dcm::DCM{T}) where {T} = dcm2rv(dcm) # Type is not specified on the LHS.

"Returns the RollPitchYaw for the given DirectionCosineMatrix."
function dcm2rpy(dcm::DCM{T}) where {T}
    # Zipfel, 75
    sin_pitch = -dcm.matrix[1, 3]
    if sin_pitch <= -one(T) + eps(one(T))
        pitch = -π/2
        roll = zero(T)
        yaw = atan(-dcm.matrix[2, 1], dcm.matrix[2, 2])
    elseif sin_pitch >= one(T) - eps(one(T))
        pitch = π/2
        roll = zero(T)
        yaw = atan(-dcm.matrix[2, 1], dcm.matrix[2, 2])
    else
        pitch = asin(clamp(sin_pitch, -one(T), one(T)))
        yaw = atan(dcm.matrix[1, 2], dcm.matrix[1, 1])
        roll = atan(dcm.matrix[2, 3], dcm.matrix[3, 3])
    end
    return RPY{T}(roll, pitch, yaw)
end
Base.convert(::Type{RPY}, dcm::DCM) = dcm2rpy(dcm) # Type is not specified on the LHS.
Base.convert(::Type{RPY{T}}, dcm::DCM{T}) where {T} = dcm2rpy(dcm) # Type is not specified on the LHS.

#############
# ERP to... #
#############

"Returns the AxisAngle for the given EulerRodriguesParameters."
function erp2aa(erp::EulerRodriguesParameters{T}) where {T}
    m = sqrt(erp.x^2 + erp.y^2 + erp.z^2)
    θ = 2 * atan(m, erp.s) # Better numerical stability than θ = 2 * acos(erp.s)
    if iszero(m)
        r = SVector{3,T}(one(T), zero(T), zero(T)) # Arbitarily choose x.
    else
        r = SVector{3,T}(erp.x/m, erp.y/m, erp.z/m)
    end
    return AA{T}(r, θ)
end
Base.convert(::Type{AA}, erp::ERP) = erp2aa(erp) # Type is not specified on the LHS.
Base.convert(::Type{AA{T}}, erp::ERP{T}) where {T} = erp2aa(erp) # Types are the same.
# Base.convert(::Type{AA{T}}, erp::ERP) = AA{T}(erp2aa(erp)) # Different types? Weird.
Base.convert(::Type{AADeg}, erp::ERP) = rad2deg(erp2aa(erp))
Base.convert(::Type{AADeg{T}}, erp::ERP{T}) where {T} = rad2deg(erp2aa(erp))

"Returns the DirectionCosineMatrix for the given EulerRodriguesParameters."
function erp2dcm(e::EulerRodriguesParameters{T}) where {T}
    x2 = e.x^2
    y2 = e.y^2
    z2 = e.z^2
    s2 = e.s^2
    m = @SMatrix [
        (x2 - y2 - z2 + s2)         2 * (e.x * e.y + e.s * e.z) 2 * (e.x * e.z - e.s * e.y);
        2 * (e.y * e.x - e.s * e.z)        (-x2 + y2 - z2 + s2) 2 * (e.y * e.z + e.s * e.x);
        2 * (e.z * e.x + e.s * e.y) 2 * (e.z * e.y - e.s * e.x)        (-x2 - y2 + z2 + s2);
    ]
    return DCM{T}(m)
end
Base.convert(::Type{DCM}, erp::ERP) = erp2dcm(erp) # Type is not specified on the LHS.
Base.convert(::Type{DCM{T}}, erp::ERP{T}) where {T} = erp2dcm(erp) # Types are the same.

"Returns the RotationVector for the given EulerRodriguesParameters."
function erp2rv(erp::EulerRodriguesParameters{T}) where {T}
    aa = erp2aa(erp)
    return RV{T}(aa.angle .* aa.axis)
end
Base.convert(::Type{RV}, erp::ERP) = erp2rv(erp) # Type is not specified on the LHS.
Base.convert(::Type{RV{T}}, erp::ERP{T}) where {T} = erp2rv(erp) # Types are the same.

"Returns the RollPitchYaw for the given EulerRodriguesParameters."
function erp2rpy(erp::EulerRodriguesParameters{T}) where {T}

    # Zipfel, 127
    q0 = erp.s
    q1 = erp.x
    q2 = erp.y
    q3 = erp.z
    tol = 1e-12 # TODO: Make this type-aware. This is for Float64.
    sin_pitch = 2 * (q0 * q2 - q1 * q3)
    if sin_pitch >= one(T) - tol
        pitch = π/2
        roll = zero(T) # Pitch is indistinguishable from yaw here.
        yaw = 2 * atan(-q1, q2) # From inspection of eq. 4.78
    elseif sin_pitch < -one(T) + tol
        pitch = -π/2
        roll = zero(T) # Pitch is indistinguishable from yaw here.
        yaw = 2 * atan(-q1, q2) # From inspection of eq. 4.78
    else
        pitch = asin(sin_pitch)
        yaw = atan(2 * (q1 * q2 + q0 * q3), q0^2 + q1^2 - q2^2 - q3^2)
        roll = atan(2 * (q2 * q3 + q0 * q1), q0^2 - q1^2 - q2^2 + q3^2)
    end

    return RPY{T}(roll, pitch, yaw)

end
Base.convert(::Type{RPY}, erp::ERP) = erp2rpy(erp) # Type is not specified on the LHS.
Base.convert(::Type{RPY{T}}, erp::ERP{T}) where {T} = erp2rpy(erp) # Types are the same.
Base.convert(::Type{RPYDeg}, erp::ERP) = rad2deg(erp2rpy(erp))
Base.convert(::Type{RPYDeg{T}}, erp::ERP{T}) where {T} = rad2deg(erp2rpy(erp))

############
# RV to... #
############

"Returns the AxisAngle for the given RotationVector."
function rv2aa(rv::RV{T}) where {T}
    θ = norm(rv.vector)
    if iszero(θ) # If no rotation...
        return zero(AA{T})
    end
    return AA{T}(rv.vector ./ θ, θ)
end
Base.convert(::Type{AA}, rv::RV) = rv2aa(rv) # Type is not specified on the LHS.
Base.convert(::Type{AA{T}}, rv::RV{T}) where {T} = rv2aa(rv) # Type is not specified on the LHS.

"Returns the DirectionCosineMatrix for the given RotationVector."
function rv2dcm(rv::RV)
    return aa2dcm(rv2aa(rv))
end
Base.convert(::Type{DCM}, rv::RV) = rv2dcm(rv) # Type is not specified on the LHS.
Base.convert(::Type{DCM{T}}, rv::RV{T}) where {T} = rv2dcm(rv) # Type is not specified on the LHS.

"Returns the EulerRodriguesParameters for the given RotationVector."
function rv2erp(rv::RV)
    return aa2erp(rv2aa(rv)) # This is already efficient.
end
Base.convert(::Type{ERP}, rv::RV) = rv2erp(rv) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, rv::RV{T}) where {T} = rv2erp(rv) # Type is not specified on the LHS.

"Returns the RollPitchYaw for the given RotationVector."
function rv2rpy(rv::RV)
    return erp2rpy(rv2erp(rv))
end
Base.convert(::Type{RPY}, rv::RV) = rv2rpy(rv) # Type is not specified on the LHS.
Base.convert(::Type{RPY{T}}, rv::RV{T}) where {T} = rv2rpy(rv) # Type is not specified on the LHS.

#############
# RPY to... #
#############

"Returns the AxisAngle for the given RollPitchYaw."
function rpy2aa(rpy::RPY)
    return erp2aa(rpy2erp(rpy))
end
Base.convert(::Type{AA}, rpy::RPY) = rpy2aa(rpy) # Type is not specified on the LHS.
Base.convert(::Type{AA{T}}, rpy::RPY{T}) where {T} = rpy2aa(rpy) # Type is not specified on the LHS.

"Returns the DirectionCosineMatrix for the given RollPitchYaw."
function rpy2dcm(rpy::RPY)
    sy, cy = sincos(rpy.yaw)
    sp, cp = sincos(rpy.pitch)
    sr, cr = sincos(rpy.roll)
    return DCM(
        @SMatrix [ # Zipfel 75, eq. 3.14
            cy * cp                    sy * cp                      -sp;
            cy * sp * sr - sy * cr     sy * sp * sr + cy * cr       cp * sr;
            cy * sp * cr + sy * sr     sy * sp * cr - cy * sr       cp * cr;
        ]
    )
end
Base.convert(::Type{DCM}, rpy::RPY) = rpy2dcm(rpy) # Type is not specified on the LHS.
Base.convert(::Type{DCM{T}}, rpy::RPY{T}) where {T} = rpy2dcm(rpy) # Type is not specified on the LHS.

"Returns the EulerRodriguesParameters for the given RollPitchYaw."
function rpy2erp(rpy::RPY)
    sy, cy = sincos(rpy.yaw / 2)
    sp, cp = sincos(rpy.pitch / 2)
    sr, cr = sincos(rpy.roll / 2)
    return ERP( # Zipfel 126, eq. 4.78
        cy * cp * sr - sy * sp * cr,
        cy * sp * cr + sy * cp * sr,
        sy * cp * cr - cy * sp * sr,
        cy * cp * cr + sy * sp * sr,
    )
end
Base.convert(::Type{ERP}, rpy::RPY) = rpy2erp(rpy) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, rpy::RPY{T}) where {T} = rpy2erp(rpy)
Base.convert(::Type{ERP}, rpy_deg::RPYDeg) = rpy2erp(deg2rad(rpy_deg)) # Type is not specified on the LHS.
Base.convert(::Type{ERP{T}}, rpy_deg::RPYDeg{T}) where {T} = rpy2erp(deg2rad(rpy_deg))

"Returns the RotationVector for the given RollPitchYaw."
function rpy2rv(rpy::RPY)
    return erp2rv(rpy2erp(rpy))
end
Base.convert(::Type{RV}, rpy::RPY) = rpy2rv(rpy) # Type is not specified on the LHS.
Base.convert(::Type{RV{T}}, rpy::RPY{T}) where {T} = rpy2rv(rpy) # Type is not specified on the LHS.

####################################
# Fallback Orientation Conversions #
####################################

Base.convert(::Type{AA{T}}, aa::AA) where {T} = AA{T}(SVector{3, T}(aa.axis), convert(T, aa.angle))
Base.convert(::Type{AADeg{T}}, aa::AADeg) where {T} = AADeg{T}(SVector{3, T}(aa.axis), convert(T, aa.angle))
Base.convert(::Type{DCM{T}}, dcm::DCM) where {T} = DCM{T}(SMatrix{3, 3, T, 9}(dcm.matrix))
Base.convert(::Type{ERP{T}}, erp::ERP) where {T} = ERP{T}(convert(T, erp.x), convert(T, erp.y), convert(T, erp.z), convert(T, erp.s))
Base.convert(::Type{RPY{T}}, rpy::RPY) where {T} = RPY{T}(convert(T, rpy.roll), convert(T, rpy.pitch), convert(T, rpy.yaw))
Base.convert(::Type{RPYDeg{T}}, rpy::RPYDeg) where {T} = RPYDeg{T}(convert(T, rpy.roll), convert(T, rpy.pitch), convert(T, rpy.yaw))
Base.convert(::Type{RV{T}}, rv::RV) where {T} = RV{T}(SVector{3, T}(rv.vector))

Base.convert(::Type{AA}, a::AbstractOrientation{T}) where {T} = convert(AA{T}, a)
Base.convert(::Type{AADeg}, a::AbstractOrientation{T}) where {T} = convert(AADeg{T}, a)
Base.convert(::Type{DCM}, a::AbstractOrientation{T}) where {T} = convert(DCM{T}, a)
Base.convert(::Type{ERP}, a::AbstractOrientation{T}) where {T} = convert(ERP{T}, a)
Base.convert(::Type{RPY}, a::AbstractOrientation{T}) where {T} = convert(RPY{T}, a)
Base.convert(::Type{RPYDeg}, a::AbstractOrientation{T}) where {T} = convert(RPYDeg{T}, a)
Base.convert(::Type{RV}, a::AbstractOrientation{T}) where {T} = convert(RV{T}, a)

function Base.convert(::Type{OT}, a::AbstractOrientation) where {T, OT <: AbstractOrientation{T}}
    erp = a isa ERP ? a : convert(ERP, a)
    return convert(OT, convert(ERP{T}, erp))
end
