using Test
using Random
using GradientOrientations
using GradientOrientations: AbstractOrientation
using LinearAlgebra
using StaticArrays

import GradientOrientations: crs3

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
                    @test interpolate(rv, ref_rv, f) ≈ erp_expected atol = erp_tol
                end

                erp_expected = inv(erp)
                @test inv(aa) ≈ erp_expected atol = erp_tol
                @test inv(dcm) ≈ erp_expected atol = erp_tol
                @test inv(rpy) ≈ erp_expected atol = erp_tol
                @test inv(rv) ≈ erp_expected atol = erp_tol

                @test aa isa AA_F64
                @test dcm isa DCM_F64
                @test erp isa ERP_F64
                @test rpy isa RPY_F64
                @test rv isa RV_F64

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

@testset "orientation isapprox" begin

    erp = normalize(ERP(1., -2., 3., -4.))
    @test erp ≈ other(erp)
    @test erp2aa(erp) ≈ erp
    @test erp2dcm(erp) ≈ erp
    @test erp2rpy(erp) ≈ erp
    @test erp2rv(erp) ≈ erp
    @test erpx(1e-8) ≈ one(ERP_F64) atol = 1.1e-8
    @test !isapprox(erpx(1e-8), one(ERP_F64); atol = 0.9e-8)

end

@testset "identity orientation" begin

    @test one(ERP_F64) == ERP(0., 0., 0., 1.)
    @test one(AA_F64) == AA(SA[1., 0., 0.], 0.)
    @test one(DCM_F64) == DCM(SMatrix{3, 3, Float64, 9}(I))
    @test one(RPY_F64) == RPY(0., 0., 0.)
    @test one(RV_F64) == RV(SA[0., 0., 0.])
    @test one(ERP) == one(ERP_F64)
    @test one(AA) == one(AA_F64)
    @test one(AADeg) == one(AADeg_F64)
    @test one(DCM) == one(DCM_F64)
    @test one(RPY) == one(RPY_F64)
    @test one(RPYDeg) == one(RPYDeg_F64)
    @test one(RV) == one(RV_F64)
    @test identity_orientation(ERP_F64) == one(ERP_F64)
    @test identity_orientation(erpx(0.2)) == one(ERP_F64)
    @test_throws MethodError zero(ERP_F64)
    @test_throws MethodError zero(AA_F64)
    @test_throws MethodError zero(DCM_F64)
    @test_throws MethodError zero(RPY_F64)
    @test_throws MethodError zero(RV_F64)

end

@testset "conversion holes" begin

    aa = AxisAngle(2. .* [1., 0., 0.], 0.5)
    @test aa isa AA_F64
    @test aa.axis == SA[2., 0., 0.]
    @test AA{Float32}([1., 0., 0.], 0.5) isa AA{Float32}

    rv = RV([1., 2., 3.])
    @test rv == RV(SA[1., 2., 3.])
    @test RV{Float32}([1., 2., 3.]) isa RV{Float32}

    m32 = Matrix{Float32}(I, 3, 3)
    @test DCM(m32) isa DCM{Float32}
    @test DCM{Float64}(m32) isa DCM_F64
    @test convert(DCM, m32) ≈ one(DCM{Float32})
    @test convert(DCM_F64, m32) ≈ one(DCM_F64)

    aa_deg = AADeg([1., 0., 0.], 90.)
    @test convert(DCM, aa_deg) ≈ erpx(π/2)
    @test convert(RV, aa_deg) ≈ erpx(π/2)
    @test convert(RPY, aa_deg) ≈ erpx(π/2)
    @test convert(DCM{Float32}, aa_deg) isa DCM{Float32}

    rpy_deg = RPYDeg(1., 2., 3.)
    @test convert(AA, rpy_deg) ≈ convert(ERP, rpy_deg)
    @test convert(DCM, rpy_deg) ≈ convert(ERP, rpy_deg)
    @test convert(RV, rpy_deg) ≈ convert(ERP, rpy_deg)
    @test convert(RV{Float32}, rpy_deg) isa RV{Float32}

    rpy32 = RPY{Float32}(0.1f0, 0.2f0, 0.3f0)
    @test convert(AA_F64, rpy32) isa AA_F64
    @test convert(DCM_F64, rpy32) isa DCM_F64
    @test convert(ERP_F64, rpy32) isa ERP_F64
    @test convert(RPY_F64, rpy32) isa RPY_F64
    @test convert(RV_F64, rpy32) isa RV_F64
    @test convert(AADeg_F64, rpy32) isa AADeg_F64
    @test convert(RPYDeg_F64, rpy32) isa RPYDeg_F64

end
