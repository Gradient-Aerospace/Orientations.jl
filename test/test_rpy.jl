@testset "RPY" begin

    a = rv2rpy(RV(SA[1., 2., 3.]))
    b = rv2rpy(RV(SA[1., 0., -1.]))
    v = SA[5., 2., -4.]
    f = 0.1

    @test rv2erp(zero(RV_F64)) ≈ zero(ERP_F64)

    rpy_tol = 1e-7

    # Note that *all* of these fall back to the ERP implementation anyway, so these aren't
    # tests so much as truisms.
    @test reframe(a, v) ≈ reframe(rpy2erp(a), v)
    @test a ⊗ b ≈ erp2rpy(rpy2erp(a) ⊗ rpy2erp(b)) atol = rpy_tol
    @test difference(a, b) ≈ erp2rpy(difference(rpy2erp(a), rpy2erp(b))) atol = rpy_tol
    @test distance(a, b) ≈ distance(rpy2erp(a), rpy2erp(b))
    @test interpolate(a, b, f) ≈ erp2rpy(interpolate(rpy2erp(a), rpy2erp(b), f)) atol = rpy_tol
    @test inv(a) ≈ erp2rpy(inv(rpy2erp(a))) atol = rpy_tol

    rng = Xoshiro(1)
    rand(rng, RPY_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test rpy2aa(a) ≈ a atol = tol
    @test rpy2dcm(a) ≈ a atol = tol
    @test rpy2erp(a) ≈ a atol = tol
    @test rpy2rv(a) ≈ a atol = tol

    # Conversion to specific types.
    rpy = RPY(0.1, 0.2, 0.3)
    dcm = rpy2dcm(rpy)
    dcm_expected = Rx(rpy.roll) * Ry(rpy.pitch) * Rz(rpy.yaw)
    @test all(dcm.matrix .≈ dcm_expected.matrix)
    erp = rpy2erp(rpy)
    erp_expected = compose(erpx(rpy.roll), erpy(rpy.pitch), erpz(rpy.yaw))
    @test all(smallest(erp) .≈ smallest(erp_expected))

    # Test conversion to/from degrees.
    deg = [10., 20., -30.]
    rpy_deg = convert(RPYDeg, RPY(deg2rad.(deg)...))
    @test rpy_deg.roll ≈ deg[1]
    @test rpy_deg.pitch ≈ deg[2]
    @test rpy_deg.yaw ≈ deg[3]
    rpy = convert(RPY, rpy_deg)
    @test rpy.roll ≈ deg2rad(deg[1])
    @test rpy.pitch ≈ deg2rad(deg[2])
    @test rpy.yaw ≈ deg2rad(deg[3])

end

@testset "RPYDeg zero and rand" begin

    z = zero(RPYDeg)
    @test z == RPYDeg(0.0, 0.0, 0.0)

    z32 = zero(RPYDeg{Float32})
    @test z32 == RPYDeg{Float32}(0f0, 0f0, 0f0)

    rng = Xoshiro(1)

    a = rand(rng, RPYDeg{Float64})
    @test a isa RPYDeg{Float64}
    @test all(isfinite, (a.roll, a.pitch, a.yaw))

    b = rand(rng, RPYDeg{Float32})
    @test b isa RPYDeg{Float32}
    @test all(isfinite, (b.roll, b.pitch, b.yaw))

end
