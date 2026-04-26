using GradientOrientations
using Test
using StaticArrays
using LinearAlgebra

@testset "RV" begin

    a = RV(SA[1., 2., 3.])
    b = RV(SA[1., 0., -1.])
    v = SA[5., 2., -4.]
    f = 0.1

    @test rv2erp(one(RV_F64)) ≈ one(ERP_F64)

    @test reframe(a, v) ≈ reframe(rv2erp(a), v)
    @test a ⊗ b ≈ erp2rv(rv2erp(a) ⊗ rv2erp(b))
    @test difference(a, b) ≈ erp2rv(difference(rv2erp(a), rv2erp(b)))
    @test distance(a, b) ≈ distance(rv2erp(a), rv2erp(b))
    @test interpolate(a, b, f) ≈ erp2rv(interpolate(rv2erp(a), rv2erp(b), f))
    @test inv(a) ≈ erp2rv(inv(rv2erp(a)))

    rng = Xoshiro(1)
    rand(rng, RV_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test rv2aa(a) ≈ a atol = tol
    @test rv2dcm(a) ≈ a atol = tol
    @test rv2erp(a) ≈ a atol = tol
    @test rv2rpy(a) ≈ a atol = tol

    vector = SA[0.1, 0.2, 0.3]
    @test distance(RV(vector)) ≈ norm(vector)
    @test distance(RV(-vector)) ≈ norm(vector)
    vector = π/4 * normalize(SA[0.1, 0.2, 0.3])
    @test distance(RV(vector)) ≈ π/4
    @test distance(RV(-vector)) ≈ π/4
    vector = (π + 0.1) * normalize(SA[0.1, 0.2, 0.3])
    @test distance(RV(vector)) ≈ 2π - (π + 0.1)
    @test distance(RV(-vector)) ≈ 2π - (π + 0.1)

end

@testset "rv2dcm" begin

    # Test 1: Zero rotation vector should give identity matrix
    rv_zero = RV(SA[0.0, 0.0, 0.0])
    dcm_identity = rv2dcm(rv_zero)
    identity_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test dcm_identity.matrix ≈ identity_matrix atol=1e-10

    # Test 2: 90 degree rotation about x-axis (π/2 * [1, 0, 0])
    rv_90x = RV(SA[π/2, 0.0, 0.0])
    dcm_90x = rv2dcm(rv_90x)
    expected_90x = @SMatrix [1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 -1.0 0.0]
    @test dcm_90x.matrix ≈ expected_90x atol=1e-10

    # Test 3: 90 degree rotation about y-axis (π/2 * [0, 1, 0])
    rv_90y = RV(SA[0.0, π/2, 0.0])
    dcm_90y = rv2dcm(rv_90y)
    expected_90y = @SMatrix [0.0 0.0 -1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]
    @test dcm_90y.matrix ≈ expected_90y atol=1e-10

    # Test 4: 90 degree rotation about z-axis (π/2 * [0, 0, 1])
    rv_90z = RV(SA[0.0, 0.0, π/2])
    dcm_90z = rv2dcm(rv_90z)
    expected_90z = @SMatrix [0.0 1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0]
    @test dcm_90z.matrix ≈ expected_90z atol=1e-10

    # Test 5: Arbitrary rotation vector
    rv_arbitrary = RV(SA[0.5, 0.3, 0.2])
    dcm_arbitrary = rv2dcm(rv_arbitrary)
    # Check that the result is a proper rotation matrix
    R = dcm_arbitrary.matrix
    identity_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test det(R) ≈ 1.0 atol=1e-10
    @test R' * R ≈ identity_matrix atol=1e-10

    # Test 6: Small angle approximation
    rv_small = RV(SA[1e-5, 2e-5, 3e-5])
    dcm_small = rv2dcm(rv_small)
    R_small = dcm_small.matrix
    identity_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test det(R_small) ≈ 1.0 atol=1e-10
    @test R_small' * R_small ≈ identity_matrix atol=1e-10

    # Test 7: Large angle (180 degrees about [1, 0, 0])
    rv_180x = RV(SA[π, 0.0, 0.0])
    dcm_180x = rv2dcm(rv_180x)
    expected_180x = @SMatrix [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
    @test dcm_180x.matrix ≈ expected_180x atol=1e-10

    # Test 8: Consistency with rv2aa and aa2dcm
    rv_test = RV(SA[0.4, 0.6, 0.2])
    aa_from_rv = rv2aa(rv_test)
    dcm_from_aa = aa2dcm(aa_from_rv)
    dcm_from_rv = rv2dcm(rv_test)
    @test dcm_from_rv.matrix ≈ dcm_from_aa.matrix atol=1e-10

    # Test 9: Type preservation (Float64)
    rv_typed = RV(SA[0.5, 0.5, 0.5])
    dcm_typed = rv2dcm(rv_typed)
    @test typeof(dcm_typed) <: DCM{Float64}
    @test typeof(dcm_typed.matrix) <: SMatrix{3, 3, Float64}

    # Test 10: Orthogonality of result
    identity_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    for angle in [0.1, 0.5, 1.0, π/4, π/2, 2.0]
        for axis_unnorm in [SA[1.0, 0.0, 0.0], SA[0.0, 1.0, 0.0], SA[1.0, 1.0, 0.0], SA[1.0, 1.0, 1.0]]
            axis = normalize(axis_unnorm)
            rv = RV(angle .* axis)
            dcm = rv2dcm(rv)
            R = dcm.matrix
            @test det(R) ≈ 1.0 atol=1e-10 #"Failed for angle=$angle, axis=$axis"
            @test R' * R ≈ identity_matrix atol=1e-10 #"Failed for angle=$angle, axis=$axis"
        end
    end

end
