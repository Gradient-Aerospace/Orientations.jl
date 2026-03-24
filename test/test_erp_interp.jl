@testset "interpolate" begin

    @testset "Basic interpolation" begin

        # Test interpolation between identity and 90-degree rotation about z
        ep_start = ERP(0., 0., 0., 1.)  # Identity
        ep_end = erpz(π/2)  # 90-degree rotation about z

        # At t=0, should get start
        result = interpolate(ep_start, ep_end, 0.0)
        @test distance(result, ep_start) ≈ 0. atol = eps(1.)

        # At t=1, should get end
        result = interpolate(ep_start, ep_end, 1.0)
        @test distance(result, ep_end) ≈ 0. atol = eps(1.)

        # At t=0.5, should get halfway rotation (45 degrees about z)
        result = interpolate(ep_start, ep_end, 0.5)
        expected = erpz(π/4)
        @test distance(result, expected) ≈ 0. atol = eps(1.)

        # At t=0.1, should get a tenth of the rotation (45 degrees about z)
        result = interpolate(ep_start, ep_end, 0.1)
        expected = erpz(π/2 / 10)
        @test distance(result, expected) ≈ 0. atol = eps(1.)

        # Check that result is normalized
        @test norm(result) ≈ 1.0

        # Just test a variety of points that don't start from identity.
        a = erpx(deg2rad(10))
        b = erpx(deg2rad(40))
        @test rad2deg(erp2aa(interpolate(a, b, 0.0)).angle) ≈ 10.
        @test rad2deg(erp2aa(interpolate(a, b, 0.1)).angle) ≈ 10. + 0.1 * 30
        @test rad2deg(erp2aa(interpolate(a, b, 0.5)).angle) ≈ 10. + 0.5 * 30
        @test rad2deg(erp2aa(interpolate(a, b, 0.9)).angle) ≈ 10. + 0.9 * 30
        @test rad2deg(erp2aa(interpolate(a, b, 1.0)).angle) ≈ 40.

        # Also, test that zero rotation doesn't break anything.
        @test distance(interpolate(a, a, 0.4), a) ≈ 0. atol = eps(2π)

    end

    @testset "Invalid t values" begin

        ep_start = ERP(0., 0., 0., 1.)
        ep_end = ERP(0., 0., 0., 1.)

        # t < 0 should throw error
        @test_throws DomainError interpolate(ep_start, ep_end, -0.1)

        # t > 1 should throw error
        @test_throws DomainError interpolate(ep_start, ep_end, 1.1)

    end

    @testset "Shortest path behavior" begin

        # Create two rotations that are more than 180 degrees apart
        ep_start = ERP(0., 0., 0., 1.)  # Identity
        ep_end_close = erpz(π/4)        # 45 degrees (close)
        ep_end_far = erpz(7π/4)         # 315 degrees (far, but same orientation)

        # With shortest_path = true (default), both should give same result
        result_close = interpolate(ep_start, ep_end_close, 0.5, shortest_path = true)
        result_far = interpolate(ep_start, ep_end_far, 0.5, shortest_path = true)

        # The results should be approximately the same (accounting for ERP ambiguity)
        # They represent the same rotation
        angle_close = distance(result_close, ep_start)
        angle_far = distance(result_far, ep_start)
        @test angle_close ≈ angle_far

    end

    @testset "Long way around behavior" begin

        ep_start = ERP(0., 0., 0., 1.)
        ep_end = other(erpz(π/2))

        # Take the long way around
        result = interpolate(ep_start, ep_end, 0.5, shortest_path = false)
        @test distance(result, ep_start) ≈ 3π/4

    end

    @testset "Interpolation with opposite quaternions" begin

        # Test that opposite quaternions (representing same rotation) are handled correctly
        ep1 = normalize(ERP(1., 2., 3., 4.))
        ep2 = other(ep1)

        result = interpolate(ep1, ep2, 0.5)
        @test distance(result, ep1) ≈ 0. atol = eps(1.)

    end

end
