module TestDimensionsExt

using Test
using StaticArrays
using Orientations
using Dimensions: Dimensions, dimstyle, numdims, getdim, eachdim, numdims_for_type

@testset "Dimensions interface for all types" begin

    aa = AxisAngle(SA[1., 0., 0.], 2.)
    aa_deg = AxisAngleDeg(SA[0., 1., 0.], 3.)
    dcm = DirectionCosineMatrix(SA[1. 0. 0.; 0. 1. 0.; 0. 0. 1.])
    erp = EulerRodriguesParameters(0., 0., 0., 1.)
    rv = RotationVector(SA[1., 0., 0.])
    rpy = RollPitchYaw(1., 2., 3.)
    rpy_deg = RollPitchYawDeg(1., 2., 3.)

    @test dimstyle(AxisAngle) == Dimensions.StructDimensionStyle()
    @test dimstyle(AxisAngleDeg) == Dimensions.StructDimensionStyle()
    @test dimstyle(DirectionCosineMatrix) == Dimensions.RealVectorDimensionStyle()
    @test dimstyle(EulerRodriguesParameters) == Dimensions.RealVectorDimensionStyle()
    @test dimstyle(RotationVector) == Dimensions.RealVectorDimensionStyle()
    @test dimstyle(RollPitchYaw) == Dimensions.StructDimensionStyle()
    @test dimstyle(RollPitchYawDeg) == Dimensions.StructDimensionStyle()

    @test numdims(aa) == 4
    @test numdims(aa_deg) == 4
    @test numdims(dcm) == 9
    @test numdims(erp) == 4
    @test numdims(rv) == 3
    @test numdims(rpy) == 3
    @test numdims(rpy_deg) == 3

    # The dimensionality is fixed and available from parameterized or unparameterized types.
    @test numdims_for_type(AxisAngle) == 4
    @test numdims_for_type(AxisAngleDeg) == 4
    @test numdims_for_type(DirectionCosineMatrix) == 9
    @test numdims_for_type(EulerRodriguesParameters) == 4
    @test numdims_for_type(RotationVector) == 3
    @test numdims_for_type(RollPitchYaw) == 3
    @test numdims_for_type(RollPitchYawDeg) == 3

    @test numdims_for_type(typeof(aa)) == 4
    @test numdims_for_type(typeof(aa_deg)) == 4
    @test numdims_for_type(typeof(dcm)) == 9
    @test numdims_for_type(typeof(erp)) == 4
    @test numdims_for_type(typeof(rv)) == 3
    @test numdims_for_type(typeof(rpy)) == 3
    @test numdims_for_type(typeof(rpy_deg)) == 3

    @test getdim(aa, 1) == 1.
    @test getdim(aa, 2) == 0.
    @test getdim(aa, 3) == 0.
    @test getdim(aa, 4) == 2.

    @test getdim(aa_deg, 1) == 0.
    @test getdim(aa_deg, 2) == 1.
    @test getdim(aa_deg, 3) == 0.
    @test getdim(aa_deg, 4) == 3.

    @test getdim(dcm, 1) == 1.
    @test getdim(dcm, 2) == 0.
    @test getdim(dcm, 3) == 0.
    @test getdim(dcm, 4) == 0.
    @test getdim(dcm, 5) == 1.
    @test getdim(dcm, 6) == 0.
    @test getdim(dcm, 7) == 0.
    @test getdim(dcm, 8) == 0.
    @test getdim(dcm, 9) == 1.

    @test getdim(erp, 1) == 0.
    @test getdim(erp, 2) == 0.
    @test getdim(erp, 3) == 0.
    @test getdim(erp, 4) == 1.

    @test getdim(rv, 1) == 1.
    @test getdim(rv, 2) == 0.
    @test getdim(rv, 3) == 0.

    @test getdim(rpy, 1) == 1.
    @test getdim(rpy, 2) == 2.
    @test getdim(rpy, 3) == 3.

    @test getdim(rpy_deg, 1) == 1.
    @test getdim(rpy_deg, 2) == 2.
    @test getdim(rpy_deg, 3) == 3.

    @test_throws BoundsError getdim(aa, 5)
    @test_throws BoundsError getdim(aa_deg, 5)
    @test_throws BoundsError getdim(dcm, 10)
    @test_throws BoundsError getdim(erp, 5)
    @test_throws BoundsError getdim(rv, 4)
    @test_throws BoundsError getdim(rpy, 4)
    @test_throws BoundsError getdim(rpy_deg, 4)

    # Each direct iterator follows the same scalar order as getdim.
    @test collect(eachdim(dcm)) == collect(dcm)
    @test collect(eachdim(erp)) == collect(erp)
    @test collect(eachdim(rv)) == collect(rv)

    @test collect(eachdim(aa)) == [1., 0., 0., 2.]
    @test collect(eachdim(aa_deg)) == [0., 1., 0., 3.]
    @test collect(eachdim(rpy)) == [1., 2., 3.]
    @test collect(eachdim(rpy_deg)) == [1., 2., 3.]

end

end
