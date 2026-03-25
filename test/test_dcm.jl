using GradientOrientations
using Test
using StaticArrays
using LinearAlgebra

@testset "DCM" begin

    a = rv2dcm(RV(SA[1., 2., 3.]))
    b = rv2dcm(RV(SA[1., 0., -1.]))
    v = SA[5., 2., -4.]
    f = 0.1

    @test dcm2erp(zero(DCM_F64)) ≈ zero(ERP_F64)

    # Triggy conversions have much worse tolerances than the non-trig types.
    dcm_tol = 1e-7

    @test reframe(a, v) ≈ reframe(dcm2erp(a), v)
    @test a ⊗ b ≈ erp2dcm(dcm2erp(a) ⊗ dcm2erp(b)) atol = dcm_tol
    @test difference(a, b) ≈ erp2dcm(difference(dcm2erp(a), dcm2erp(b)))
    @test distance(a, b) ≈ distance(dcm2erp(a), dcm2erp(b))
    @test interpolate(a, b, f) ≈ erp2dcm(interpolate(dcm2erp(a), dcm2erp(b), f)) atol = dcm_tol
    @test inv(a) ≈ erp2dcm(inv(dcm2erp(a))) atol = dcm_tol

    rand(DCM_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test dcm2aa(a) ≈ a atol = tol
    @test dcm2erp(a) ≈ a atol = tol
    @test dcm2rpy(a) ≈ a atol = tol
    @test dcm2rv(a) ≈ a atol = tol

    dcm = Rx(0.1)
    @test distance(dcm) ≈ 0.1
    dcm = Ry(-0.2)
    @test distance(dcm) ≈ 0.2
    dcm = Rz(3π/2)
    @test distance(dcm) ≈ π/2

end

@testset "dcm2aa" begin

    # Test 1: Identity matrix should give zero angle
    R_identity = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    dcm_identity = DCM(R_identity)
    aa_identity = dcm2aa(dcm_identity)
    @test aa_identity.angle ≈ 0.0 atol=1e-10
    @test norm(aa_identity.axis) == 1.

    # Test 2: 90 degree rotation about x-axis
    aa_90x = AxisAngle(SA[1.0, 0.0, 0.0], π/2)
    dcm_90x = aa2dcm(aa_90x)
    aa_from_dcm = dcm2aa(dcm_90x)
    @test aa_from_dcm.angle ≈ π/2 atol=1e-10
    @test norm(aa_from_dcm.axis .- SA[1.0, 0.0, 0.0]) < 1e-10

    # Test 3: 90 degree rotation about y-axis
    aa_90y = AxisAngle(SA[0.0, 1.0, 0.0], π/2)
    dcm_90y = aa2dcm(aa_90y)
    aa_from_dcm_y = dcm2aa(dcm_90y)
    @test aa_from_dcm_y.angle ≈ π/2 atol=1e-10
    @test norm(aa_from_dcm_y.axis .- SA[0.0, 1.0, 0.0]) < 1e-10

    # Test 4: 90 degree rotation about z-axis
    aa_90z = AxisAngle(SA[0.0, 0.0, 1.0], π/2)
    dcm_90z = aa2dcm(aa_90z)
    aa_from_dcm_z = dcm2aa(dcm_90z)
    @test aa_from_dcm_z.angle ≈ π/2 atol=1e-10
    @test norm(aa_from_dcm_z.axis .- SA[0.0, 0.0, 1.0]) < 1e-10

    # Test 5: Arbitrary rotation (45 degrees about [1,1,1] axis)
    axis_arbitrary = normalize(SA[1.0, 1.0, 1.0])
    angle_arbitrary = π/4
    aa_arbitrary = AxisAngle(axis_arbitrary, angle_arbitrary)
    dcm_arbitrary = aa2dcm(aa_arbitrary)
    aa_from_dcm_arb = dcm2aa(dcm_arbitrary)
    @test aa_from_dcm_arb.angle ≈ angle_arbitrary atol=1e-10
    @test norm(aa_from_dcm_arb.axis .- axis_arbitrary) < 1e-10

    # Test 6: Large angle (180 degrees), and mostly rotated about x
    axis_180 = normalize(SA[2.0, 0.0, 1.0])
    angle_180 = 1π
    aa_180 = AxisAngle(axis_180, angle_180)
    dcm_180 = aa2dcm(aa_180)
    aa_from_dcm_180 = dcm2aa(dcm_180)
    # This rotation could go either way.
    if aa_from_dcm_180.angle < 0
        @test aa_from_dcm_180.angle ≈ angle_180 atol=1e-10
        @test norm(aa_from_dcm_180.axis .- -axis_180) < 1e-10
    else
        @test aa_from_dcm_180.angle ≈ angle_180 atol=1e-10
        @test norm(aa_from_dcm_180.axis .- axis_180) < 1e-10
    end

    # Test 6(b): Large angle (not really less than 180 degrees), and mostly rotated about y
    axis_180 = normalize(SA[2.0, -3.0, -1.0])
    angle_180 = 1π - eps(1π) / 2
    aa_180 = AxisAngle(axis_180, angle_180)
    dcm_180 = aa2dcm(aa_180)
    aa_from_dcm_180 = dcm2aa(dcm_180)
    # This rotation could go either way.
    if aa_from_dcm_180.angle < 0
        @test aa_from_dcm_180.angle ≈ angle_180 atol=1e-10
        @test norm(aa_from_dcm_180.axis .- -axis_180) < 1e-10
    else
        @test aa_from_dcm_180.angle ≈ angle_180 atol=1e-10
        @test norm(aa_from_dcm_180.axis .- axis_180) < 1e-10
    end

    # Test 6(c): Large angle (just less than 180 degrees), and mostly rotated about z
    axis_180 = normalize(SA[1.0, 2.0, 3.0])
    angle_180 = 1π - 2 * eps(1π)
    aa_180 = AxisAngle(axis_180, angle_180)
    dcm_180 = aa2dcm(aa_180)
    aa_from_dcm_180 = dcm2aa(dcm_180)
    @test aa_from_dcm_180.angle ≈ angle_180 atol=1e-10
    @test norm(aa_from_dcm_180.axis .- axis_180) < 1e-10

    # Test 7: Small angle (very small rotation)
    axis_small = normalize(SA[1.0, 2.0, 3.0])
    angle_small = 1e-5
    aa_small = AxisAngle(axis_small, angle_small)
    dcm_small = aa2dcm(aa_small)
    aa_from_dcm_small = dcm2aa(dcm_small)
    @test aa_from_dcm_small.angle ≈ angle_small atol=1e-8
    @test norm(aa_from_dcm_small.axis .- axis_small) < 1e-6

    # Test 8: Round-trip conversion - AA -> DCM -> AA
    aa_original = AxisAngle(normalize(SA[1.0, 2.0, 3.0]), 0.7)
    dcm_intermediate = aa2dcm(aa_original)
    aa_final = dcm2aa(dcm_intermediate)
    @test aa_final.angle ≈ aa_original.angle atol=1e-10
    @test norm(aa_final.axis .- aa_original.axis) < 1e-10

    # Test 9: Type preservation (Float64)
    aa_typed = AxisAngle(normalize(SA[1.0, 1.0, 0.0]), 0.5)
    dcm_typed = aa2dcm(aa_typed)
    aa_from_dcm_typed = dcm2aa(dcm_typed)
    @test typeof(aa_from_dcm_typed) <: AxisAngle{Float64}
    @test typeof(aa_from_dcm_typed.axis) <: SVector{3, Float64}
    @test typeof(aa_from_dcm_typed.angle) <: Float64

end

@testset "dcm2aa numerical stability near pi" begin

    # Baseline (the previous skew-only implementation) can be inaccurate near pi because
    # R23-R32, R31-R13, and R12-R21 collapse toward zero as sin(theta) -> 0. For this test
    # case, the baseline error is around 1e-3 rad.
    function dcm2aa_baseline(dcm::DCM{T}) where {T}
        R = dcm.matrix
        cos_angle = (tr(R) - one(T)) / T(2)
        if cos_angle <= -one(T)
            cos_angle = -one(T)
        elseif cos_angle >= one(T)
            return AA{T}(SA[one(T), zero(T), zero(T)], zero(T))
        end
        angle = acos(cos_angle)
        axis = normalize(
            SVector{3, T}(
                R[2, 3] - R[3, 2],
                R[3, 1] - R[1, 3],
                R[1, 2] - R[2, 1],
            ),
        )
        return AA{T}(axis, angle)
    end

    axis = normalize(SA[0.29375010045565614, 0.33943374564517836, -0.8935858161360756])
    angle = 3.141592653589759
    # pi is 3.141592653589793, so we're testing slightly less.

    aa_true = AA(axis, angle)
    dcm = aa2dcm(aa_true)

    aa_baseline = dcm2aa_baseline(dcm)
    aa_improved = dcm2aa(dcm)

    baseline_err = distance(aa_baseline, aa_true)
    improved_err = distance(aa_improved, aa_true)

    @test baseline_err > 1e-4
    @test improved_err < baseline_err
    @test improved_err < 1e-8

end
