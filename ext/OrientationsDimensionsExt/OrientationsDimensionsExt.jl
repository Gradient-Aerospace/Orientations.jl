module OrientationsDimensionsExt

import Dimensions
using Dimensions: RealVectorDimensionStyle, StructDimensionStyle
using Orientations: AxisAngle, AxisAngleDeg, DirectionCosineMatrix,
    EulerRodriguesParameters, RollPitchYaw, RollPitchYawDeg, RotationVector

# Here, we expose the scalar dimensions of each of the orientation types. The simplest
# approach is to simply treat each field of each type as a dimension. In this sense, the
# RotationVector type would be 3-dimensional, the EulerRodriguesParameters type would be
# 4-dimensional, and the DirectionCosineMatrix type would be 9-dimensional. Another approach
# would be to reduce all of these to just 3 dimensions (the RotationVector type) or to
# reduce them to 4 dimensions (EulerRodriguesParameters), since these are the only true
# degrees of freedom (depending on exactly what you want from the notion of a "dimension" in
# this context).
#
# Which should we chose? The notion of the Dimensions package isn't to reduce everything to
# a minimal form. Rather, it's just trying to make it easy to access the nth dimension of a
# type. In that context, it seems to make the most sense to expose the underlying data as
# the dimensions, not the underlying meaning of an orientation. That is the desing choice
# we've made here. Then, if you want a reduced form, simply convert to RotationVector or
# EulerRodriguesParameters first before asking for dimensions.

# These behave like structs. To get each scalar, we need to go through each field, in order.
Dimensions.dimstyle(::Type{<:AxisAngle}) = StructDimensionStyle()
Dimensions.dimstyle(::Type{<:AxisAngleDeg}) = StructDimensionStyle()
Dimensions.dimstyle(::Type{<:RollPitchYaw}) = StructDimensionStyle()
Dimensions.dimstyle(::Type{<:RollPitchYawDeg}) = StructDimensionStyle()
# TODO: This could be made more efficient.

# These behave more like vectors because they support getindex.
Dimensions.dimstyle(::Type{<:EulerRodriguesParameters}) = RealVectorDimensionStyle()
Dimensions.dimstyle(::Type{<:RotationVector}) = RealVectorDimensionStyle()
Dimensions.dimstyle(::Type{<:DirectionCosineMatrix}) = RealVectorDimensionStyle()

# RealVectorDimensionStyle normally uses `axes(x, 1)` to support vectors with non-traditional
# indices. These orientation types use traditional one-based indexing but aren't
# AbstractArrays, so access their dimensions directly rather than adding `axes` methods just
# for Dimensions.
Dimensions.getdim(erp::EulerRodriguesParameters, d) = getindex(erp, d)
Dimensions.getdim(rv::RotationVector, d) = getindex(rv, d)
Dimensions.getdim(dcm::DirectionCosineMatrix, d) = getindex(dcm, d)

end
