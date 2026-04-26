@testset "AA" begin

    a = AxisAngle(normalize(SA[1., 2., 3.]), 1.1)
    b = AxisAngle(normalize(SA[1., 0., -1.]), 2.2)
    v = SA[5., 2., -4.]
    f = 0.1

    @test aa2erp(zero(AA_F64)) ≈ zero(ERP_F64)

    @test reframe(a, v) ≈ reframe(aa2erp(a), v)
    @test a ⊗ b ≈ erp2aa(aa2erp(a) ⊗ aa2erp(b))
    @test difference(a, b) ≈ erp2aa(difference(aa2erp(a), aa2erp(b)))
    @test distance(a, b) ≈ distance(aa2erp(a), aa2erp(b))
    @test interpolate(a, b, f) ≈ erp2aa(interpolate(aa2erp(a), aa2erp(b), f))
    @test inv(a) ≈ erp2aa(inv(aa2erp(a)))

    rng = Xoshiro(1)
    rand(rng, AA_F64) # Just test that we can do it.

    # Test conversions to all other types.
    tol = 1e-7
    @test aa2dcm(a) ≈ a atol = tol
    @test aa2erp(a) ≈ a atol = tol
    @test aa2rpy(a) ≈ a atol = tol
    @test aa2rv(a) ≈ a atol = tol

    @test distance(AA(SA[1., 0., 0.], 3.1)) ≈ 3.1
    @test distance(AA(SA[1., 0., 0.], 3.2)) ≈ 2π - 3.2

    c = smallest(AA(SA[1., 0., 0.], 3.1))
    @test c.angle == 3.1
    @test c.axis == SA[1., 0., 0.]
    c = smallest(AA(SA[1., 0., 0.], 3π/2))
    @test c.angle ≈ π/2
    @test c.axis == SA[-1., 0., 0.]

    # Really small rotation
    a = AA(SA[0., 1., 0.], 1e-8)
    erp = aa2erp(a)
    b = erp2aa(erp)
    @test a.axis ≈ b.axis atol = 1e-15
    @test a.angle ≈ b.angle atol = 1e-15

end

@testset "AADeg zero and rand" begin

    z = zero(AADeg)
    @test z == AADeg(SA[1.0, 0.0, 0.0], 0.0)

    z32 = zero(AADeg{Float32})
    @test z32 == AADeg{Float32}(SA[1f0, 0f0, 0f0], 0f0)

    rng = Xoshiro(1)
    a = rand(rng, AADeg{Float64})
    @test a isa AADeg{Float64}
    @test isapprox(norm(a.axis), 1.0; atol = 1e-12)

    b = rand(rng, AADeg{Float32})
    @test b isa AADeg{Float32}
    @test isapprox(norm(b.axis), 1f0; atol = 1f-5)

end

@testset "AADeg AbstractOrientation operations" begin

    a = AADeg(SA[1., 0., 0.], 15.)
    b = AADeg(normalize(SA[1., 2., -1.]), 65.)
    v = SA[5., 2., -4.]
    f = 0.35
    tol = 1e-7

    @test convert(AADeg_F64, convert(ERP_F64, a)) isa AADeg_F64
    @test convert(AADeg_F64, convert(ERP_F64, a)) ≈ a atol = tol

    @test reframe(a, v) ≈ reframe(deg2rad(a), v)

    c = compose(a, b)
    @test c isa AADeg_F64
    @test c ≈ compose(deg2rad(a), deg2rad(b)) atol = tol

    d = difference(a, b)
    @test d isa AADeg_F64
    @test d ≈ difference(deg2rad(a), deg2rad(b)) atol = tol

    @test distance(a) ≈ distance(deg2rad(a))
    @test distance(a, b) ≈ distance(deg2rad(a), deg2rad(b))

    i = interpolate(a, b, f)
    @test i isa AADeg_F64
    @test i ≈ interpolate(deg2rad(a), deg2rad(b), f) atol = tol

    i_long = interpolate(a, b, f; shortest_path = false)
    @test i_long isa AADeg_F64
    @test i_long ≈ interpolate(convert(ERP_F64, a), convert(ERP_F64, b), f; shortest_path = false) atol = tol

    a_inv = inv(a)
    @test a_inv isa AADeg_F64
    @test a_inv ≈ inv(deg2rad(a)) atol = tol

    a32 = AADeg{Float32}(SA[1f0, 0f0, 0f0], 15f0)
    b32 = AADeg{Float32}(SA[0f0, 1f0, 0f0], 45f0)
    @test interpolate(a32, b32, 0.25f0) isa AADeg{Float32}

end
