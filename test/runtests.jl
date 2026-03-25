using Test
using Random
using GradientOrientations
using GradientOrientations: AbstractOrientation
using LinearAlgebra
using StaticArrays

import GradientOrientations: crs3

# By default, this is pretty tight.
Base.isapprox(a::ERP, b::ERP; atol = eps(4π)) = distance(normalize(a), normalize(b)) <= atol

# Obviously, by relying on ERP conversion to test all orientations, we should do a good job
# of making sure the ERP stuff is right on its own (separate from these tests).
function Base.isapprox(a::AbstractOrientation, b::AbstractOrientation; kwargs...)
    return isapprox(
        convert(ERP, a),
        convert(ERP, b);
        kwargs...
    )
end

include("test_aa.jl")
include("test_dcm.jl")
include("test_erp.jl")
include("test_rv.jl")
include("test_rpy.jl")

# Now that we've tested ERPs, let's test the other types by converting to/from ERPs.

@testset "all operations on all types" begin

    # For binary operations, we'll need some reference to use.
    ref_erp = normalize(ERP(1., 2., 3., 4.))
    ref_aa = erp2aa(ref_erp)
    ref_dcm = erp2dcm(ref_erp)
    ref_rpy = erp2rpy(ref_erp)
    ref_rv = erp2rv(ref_erp)

    v = SA[2., -3., 4.]

    v_tol = 1e-6
    erp_tol = 1e-6

    # Let's construct a bunch of different orientations and see that operations using them
    # Doing this with roll, pitch, and yaw let's us test a lot of annoying conditions
    # directly.
    for yaw in LinRange(-2π, 2π, 21)
        for pitch in LinRange(-π/2, π/2, 11)
            for roll in LinRange(-π, π, 21)

                rpy = RPY(roll, pitch, yaw)
                erp = rpy2erp(rpy)
                aa = erp2aa(erp)
                dcm = erp2dcm(erp)
                rv = erp2rv(erp)

                v_expected = reframe(erp, v)
                @test reframe(aa, v) ≈ v_expected atol = v_tol
                @test reframe(dcm, v) ≈ v_expected atol = v_tol
                @test reframe(rpy, v) ≈ v_expected atol = v_tol
                @test reframe(rv, v) ≈ v_expected atol = v_tol

                erp_expected = compose(erp, ref_erp)
                @test compose(aa, ref_aa) ≈ erp_expected atol = erp_tol
                @test compose(dcm, ref_dcm) ≈ erp_expected atol = erp_tol
                @test compose(rpy, ref_rpy) ≈ erp_expected atol = erp_tol
                @test compose(rv, ref_rv) ≈ erp_expected atol = erp_tol

                erp_expected = difference(erp, ref_erp)
                @test difference(aa, ref_aa) ≈ erp_expected atol = erp_tol
                @test difference(dcm, ref_dcm) ≈ erp_expected atol = erp_tol
                @test difference(rpy, ref_rpy) ≈ erp_expected atol = erp_tol
                @test difference(rv, ref_rv) ≈ erp_expected atol = erp_tol

                angle_expected = distance(erp, ref_erp)
                @test distance(aa, ref_aa) ≈ angle_expected atol = erp_tol
                @test distance(dcm, ref_dcm) ≈ angle_expected atol = erp_tol
                @test distance(rpy, ref_rpy) ≈ angle_expected atol = erp_tol
                @test distance(rv, ref_rv) ≈ angle_expected atol = erp_tol

                for f in LinRange(0., 1., 5)
                    erp_expected = interpolate(erp, ref_erp, f)
                    @test interpolate(aa, ref_aa, f) ≈ erp_expected atol = erp_tol
                    @test interpolate(dcm, ref_dcm, f) ≈ erp_expected atol = erp_tol
                    @test interpolate(rpy, ref_rpy, f) ≈ erp_expected atol = erp_tol
                    if !isapprox(interpolate(rpy, ref_rpy, f), erp_expected; atol = erp_tol)
                        @show rpy
                        @show ref_rpy
                        @show erp
                        @show ref_erp
                        @show f
                    end
                    @test interpolate(rv, ref_rv, f) ≈ erp_expected atol = erp_tol
                end

                erp_expected = inv(erp)
                @test inv(aa) ≈ erp_expected atol = erp_tol
                @test inv(dcm) ≈ erp_expected atol = erp_tol
                @test inv(rpy) ≈ erp_expected atol = erp_tol
                @test inv(rv) ≈ erp_expected atol = erp_tol

                eltype(aa) == Float64
                eltype(dcm) == Float64
                eltype(erp) == Float64
                eltype(rpy) == Float64
                eltype(rv) == Float64

            end
        end
    end

end

@testset "specialty conversions" begin

    rpy = RPYDeg_F64(0., 0., 90.)
    erp = convert(ERP, rpy)
    @test erp ≈ erpz(π/2)

    aa = AADeg_F64(SA[1., 0., 0.], 90)
    erp = convert(ERP, aa)
    @test erp ≈ erpx(π/2)

end
