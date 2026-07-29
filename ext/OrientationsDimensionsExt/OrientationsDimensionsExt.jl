module OrientationsDimensionsExt

import Dimensions
using Dimensions: RealVectorDimensionStyle, StructDimensionStyle
using Orientations: AxisAngle, AxisAngleDeg, DirectionCosineMatrix,
    EulerRodriguesParameters, RollPitchYaw, RollPitchYawDeg, RotationVector

# Here, we expose the scalar dimensions of each orientation type. The simplest approach is
# to treat each stored scalar as a dimension. In this sense, RotationVector is
# 3-dimensional, EulerRodriguesParameters is 4-dimensional, and DirectionCosineMatrix is
# 9-dimensional. Another approach would be to reduce all of these to just 3 dimensions (the
# RotationVector type) or to reduce them to 4 dimensions (EulerRodriguesParameters), since
# these are the only true degrees of freedom (depending on exactly what you want from the
# notion of a "dimension" in this context).
#
# Which should we choose? The purpose of Dimensions isn't to reduce everything to a minimal
# form. Rather, it's to make it easy to access the nth dimension of a type. In that context,
# it seems to make the most sense to expose the underlying data as the dimensions rather
# than the underlying meaning of an orientation. Then, if you want a reduced form, simply
# convert to RotationVector or EulerRodriguesParameters first.

# AxisAngle dimensions are the three scalar axis components followed by the angle. Its
# indexing behavior is intentionally different: index 1 is the complete axis.
Dimensions.dimstyle(::Type{<:AxisAngle}) = StructDimensionStyle()
Dimensions.numdims(::AxisAngle) = 4
Dimensions.numdims_for_type(::Type{<:AxisAngle}) = 4
Dimensions.getdim(aa::AxisAngle, d) = getindex((aa.axis..., aa.angle), d)
Dimensions.eachdim(aa::AxisAngle) = (aa.axis..., aa.angle)

# AxisAngleDeg has the same scalar layout as AxisAngle.
Dimensions.dimstyle(::Type{<:AxisAngleDeg}) = StructDimensionStyle()
Dimensions.numdims(::AxisAngleDeg) = 4
Dimensions.numdims_for_type(::Type{<:AxisAngleDeg}) = 4
Dimensions.getdim(aa::AxisAngleDeg, d) = getindex((aa.axis..., aa.angle), d)
Dimensions.eachdim(aa::AxisAngleDeg) = (aa.axis..., aa.angle)

# DirectionCosineMatrix dimensions follow its column-major linear indexing.
Dimensions.dimstyle(::Type{<:DirectionCosineMatrix}) = RealVectorDimensionStyle()
Dimensions.numdims(::DirectionCosineMatrix) = 9
Dimensions.numdims_for_type(::Type{<:DirectionCosineMatrix}) = 9
Dimensions.getdim(dcm::DirectionCosineMatrix, d) = getindex(dcm, d)
Dimensions.eachdim(dcm::DirectionCosineMatrix) = ntuple(d -> getindex(dcm, d), 9)

# EulerRodriguesParameters dimensions follow the x, y, z, scalar indexing order.
Dimensions.dimstyle(::Type{<:EulerRodriguesParameters}) = RealVectorDimensionStyle()
Dimensions.numdims(::EulerRodriguesParameters) = 4
Dimensions.numdims_for_type(::Type{<:EulerRodriguesParameters}) = 4
Dimensions.getdim(erp::EulerRodriguesParameters, d) = getindex(erp, d)
Dimensions.eachdim(erp::EulerRodriguesParameters) = (erp.x, erp.y, erp.z, erp.s)

# RotationVector dimensions follow its x, y, z indexing order.
Dimensions.dimstyle(::Type{<:RotationVector}) = RealVectorDimensionStyle()
Dimensions.numdims(::RotationVector) = 3
Dimensions.numdims_for_type(::Type{<:RotationVector}) = 3
Dimensions.getdim(rv::RotationVector, d) = getindex(rv, d)
Dimensions.eachdim(rv::RotationVector) = Tuple(rv.vector)

# RollPitchYaw dimensions follow its roll, pitch, yaw indexing order.
Dimensions.dimstyle(::Type{<:RollPitchYaw}) = StructDimensionStyle()
Dimensions.numdims(::RollPitchYaw) = 3
Dimensions.numdims_for_type(::Type{<:RollPitchYaw}) = 3
Dimensions.getdim(rpy::RollPitchYaw, d) = getindex(rpy, d)
Dimensions.eachdim(rpy::RollPitchYaw) = (rpy.roll, rpy.pitch, rpy.yaw)

# RollPitchYawDeg has the same scalar layout as RollPitchYaw.
Dimensions.dimstyle(::Type{<:RollPitchYawDeg}) = StructDimensionStyle()
Dimensions.numdims(::RollPitchYawDeg) = 3
Dimensions.numdims_for_type(::Type{<:RollPitchYawDeg}) = 3
Dimensions.getdim(rpy::RollPitchYawDeg, d) = getindex(rpy, d)
Dimensions.eachdim(rpy::RollPitchYawDeg) = (rpy.roll, rpy.pitch, rpy.yaw)

end
