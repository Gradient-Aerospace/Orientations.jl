# GradientOrientations.jl

*NOTE*: This is a sandbox for now. If we like this package, we should release it as open source under Orientations.jl.

This package is useful for specifying the orientation of a thing relative to another thing and for operating on those orientations.

The available orientation types are:

* `AxisAngle` (aka `AA`)
* `DirectionCosineMatrix` (aka `DCM`)
* `EulerRodriguesParameters` (aka `ERP`)
* `RotationVector` (aka `RV`)
* `RollPitchYaw` (aka `RPY`)
* `AxisAngleDeg` (aka `AADeg`)
* `RollPitchYawDeg` (aka `RPYDeg`)

## Basic Example

Suppose that frame B is rotated 45 degrees from frame A about the (common) `x` axis. Here's how we specify that, using the AxisAngle orientation type as an example:

```julia
using GradientOrientations
using StaticArrays # For the static vector type
aa_b_wrt_a = AxisAngle(SA[1., 0., 0.], deg2rad(45))
```

Suppose we know how some vector is expressed in frame A. The way that vector would be expressed in frame B is given by the `reframe` function:

```julia
v_in_b = reframe(aa_b_wrt_a, v_in_a)
```

In this sense, orientations imply frame rotation (passive rotations).

If frame C is rotated 50 degrees from frame B about the axis `1/√3 * [1., 1., 1.]` (in either B or C, the rotation axis being necessarily common to both), we could specify frame C with respect to (wrt) frame A as using the `compose` function.

```julia
aa_c_wrt_b = AxisAngle(1/√3 * [1., 1., 1.], deg2rad(50))
aa_c_wrt_a = compose(aa_c_wrt_b, aa_b_wrt_a)
```

and of course we now have:

```julia
v_in_c = reframe(aa_c_wrt_a, v_in_a)
```

All orientation types behave in exactly this way, with the same order of operations, etc. In general, you don't have to even know which type of orientation you're working with to use it correctly.

## Conversions

Converting from one type to another can be accomplished with the `Base.convert` function, as in:

```julia
dcm_b_wrt_a = convert(DCM, aa_b_wrt_a)
erp_b_wrt_a = convert(ERP, aa_b_wrt_a)
rv_b_wrt_a  = convert(RV,  aa_b_wrt_a)
rpy_b_wrt_a = convert(RPY, aa_b_wrt_a)
```

Those `convert` methods actually call the underlying functions for conversion.

```julia
dcm_b_wrt_a = aa2dcm(aa_b_wrt_a)
erp_b_wrt_a = aa2erp(aa_b_wrt_a)
rv_b_wrt_a = aa2rv(aa_b_wrt_a)
rpy_b_wrt_a = aa2rpy(aa_b_wrt_a)
```

Similar functions with obvious names exist to convert between all available types.

Where a direct, numerically stable conversion from one type directly to another exists, that will be the underlying method. Where a direct conversion does not exist, the right-hand type will be converted to Euler-Rodrigues Parameters and then from there to the desired type. (All orientation types are expected to have a path to `EulerRodriguesParameters`.)

Note that the RollPitchYaw type implies that the intended frame is oriented `roll` around the x axis of a frame that's oriented `pitch` around the y axis of a frame that's oriented `yaw` around the z axis of the reference frame. That is, the following are equivalent:

```julia
rpy_b_wrt_a = RPY(0.1, 0.2, 0.3)
erp_b_wrt_a = compose(erpx(0.1), erpy(0.2), erpz(0.3))
rpy_b_wrt_a = erp2rpy(erp_b_wrt_a)
```

## Operations

* `reframe(b_wrt_a, v_a)`: If `v_a` is a vector expressed in frame A and `b_wrt_a` is the orientation of frame B wrt frame A (any type of orientation will do), then this returns the same vector expressed in B.
* `compose(c_wrt_b, b_wrt_a)`: Returns the orientation of C wrt A with the same type as the inputs.
* `difference(c_wrt_a, b_wrt_a)`: Returns the orientation of C wrt B with the same type as the inputs.
* `distance`: Returns the rotation angle of the orientation, the "smallest way around", in radians.
* `interpolate(o1, o2, f)`: Interpolates (spherically) from orientation 1 to orientation 2 (both wrt the same reference) using `f` in the inclusive range [0, 1].
* `Base.inv(b_wrt_a)`: Inverts the orientation, returning A wrt B.
* `Base.zero(type)`: Returns an orientation of the given type with zero rotation from the reference.
* `Random.rand(rng, type)`: Returns a random orientation of the given type (any of the available orientation types) drawn uniformly from SO(3).

The `⊗` operator (`\otimes`) can also be used for composition. That is, `a ⊗ b == compose(a, b)`.

For `EulerRodriguesParameters`:

* `rate(erp_b_wrt_a, omega_b_wrt_a_in_b)`: Returns the derivative over time of the given ERP using the given rotation rate. The rotation rate must be expressed in frame B; if you have `omega_b_wrt_a_in_a`, convert it first with `reframe(erp_b_wrt_a, omega_b_wrt_a_in_a)`. The return type will be an ERP type, though obviously it will not have a unit norm and will not, itself, represent an orientation. This is useful for numerical integration. While all types could conceivably implement this function, this package only implements it for `EulerRodriguesParameters` since it has the best numerical properties.
* `other`: Returns an ERP that rotates "the other way around".
* `smallest`: Returns an ERP that is equivalent to the input but is the "shortest way around".
* `LinearAlgebra.normalize`: Returns a normalized version of an ERP -- useful after operations such as numerical integration.

The `DirectionCosineMatrix` also functions like a matrix in that it supports the `*` operator. So you can compose DCMs using *, and you can "reframe" a vector using `*` as well.

## Construction

Here are examples of constructing each type.

For `AxisAngle`, provide the axis first and then the angle:

```julia
AA(SA[1., 0., 0.], 0.)
```

The constructor stores the axis exactly as provided; it does not normalize the axis.
Conversions and orientation operations that require a unit axis normalize internally, and
`LinearAlgebra.normalize(aa)` returns an equivalent `AxisAngle` with a unit axis.

For convenience, `aax`, `aay`, and `aaz` exist to create orientations rotated about a given primary axis.

For `DirectionCosineMatrix`, give the matrix itself:

```julia
DCM(
    @SMatrix [
        1. 0. 0.;
        0. 1. 0.;
        0. 0. 1.;
    ]
)
```

For convenience, `Rx`, `Ry`, and `Rz` exist to create DCMs rotated about a given primary axis.

For `EulerRodriguesParameters`, give each element. The scalar is last. The following corresponds to "no rotation":

```julia
ERP(0., 0., 0., 1.)
ERP(; x = 0., y = 0., z = 0., s = 1.)
```

Similarly, `erpx`, `erpy`, and `erpz` all exist.

For `RotationVector`, give the vector:

```julia
RV(SA[0., 0., 0.])
```

For `RollPitchYaw`, provide roll, pitch, and yaw:

```julia
RPY(0., 0., 0.)
RPY(; roll = 0., pitch = 0., yaw = 0.)
```

## Human Interface Types

There are a couple of types that store their angles in degrees for human input and output. For example, `RPYDeg` lets the user specify roll, pitch, and yaw in degrees, and `AADeg` does the same for axis-angle orientations.

These types support the same `AbstractOrientation` operations as their radian counterparts, and operations that return orientations preserve the degree-valued type. They can also be converted to their radian forms via `convert(RPY, my_rpy_deg)`, `convert(AA, my_aa_deg)`, or `deg2rad`. The `distance` functions still return radians, consistent with the rest of the package.

## Where Are The Quaternions?

The reader may notice that there is no "quaternion" type. This is just nomenclature. The word "quaternion" is used with so many different conventions that any time the word "quaternion" appears on an interface between two systems, a healthy discussion (and usually several examples) are necessary to describe what is actually meant. This package chooses to bypass that nomenclatural minefield. Instead, this package uses "Euler-Rodrigues Parameters", which is unambiguously described by Shuster in "A Survey of Attitude Representations", freely available [here](https://www.malcolmdshuster.com/Doorway_Pubs-1970-1998.htm). This type, which he also calls "the quaternion of rotation" (JPL conventions), has consistent rules and conversions to the other types implemented here. More to the point, this type has all of the advantages of (in fact, *is*) a quaternion of rotation (a minimal representation of SO(3), numerical stability, good behavior under numerical integration, good conversion to and from other types) without the confusion surrounding the conventions (scalar first or last? is vector or frame rotation implied? what are the rules of composition, and is that the same thing as quaternion multiplication? if you right-multiply a vector, is that the same as right-multiplying the same vector with a DCM obtained from the conversion of the quaternion to a DCM?).

### Exporting to Eigen

Eigen is a C++ library for linear algebra, and it features a Quaternion type. The following two snippets give the same results:

```julia
using GradientOrientations
using StaticArrays
erp_b_wrt_a = ERP(x = -0.39015349272073274, y = 0.44885619975114965, z = -0.5331968363460537, s = 0.6016722511245876)
v_a = SA[ -0.09342155668355605, -0.5691321723949833, 0.18864021876469472]
v_b = reframe(erp_b_wrt_a, v_a)
```

This gives `SA[0.5384388474664465, -0.27833775489103457, -0.02891108412057084]`.

For eigen, the scalar moves first, and the x, y, and z components are negated (expressing the inverted rotation we'd expect for a library that implements vector/active rotations).

```cpp
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

int main() {
  Eigen::Quaterniond q_a_to_b(0.6016722511245876, 0.39015349272073274, -0.44885619975114965, 0.5331968363460537);
  Eigen::Vector3d v_a(-0.09342155668355605, -0.5691321723949833, 0.18864021876469472);
  Eigen::Vector3d v_b = q_a_to_b * v_a;
  std::cout << v_b << std::endl;
}
```

This gives the same output as the above.

The composition rules then follow the same pattern:

```julia
erp_c_wrt_b = ERP(x = -0.26328146006533937, y = 0.3194557972432347, z = -0.547884663730034, s = 0.7269479084796793)
erp_c_wrt_a = compose(erp_c_wrt_b, erp_b_wrt_a)
v_c = reframe(erp_c_wrt_a, v_a)
```

which gives `SA[0.37890885798068275, 0.28686794099866375, 0.37730493081162575]`.

For Eigen, we add the following to the end of `main`:

```cpp
  Eigen::Quaterniond q_b_to_c(0.7269479084796793, 0.26328146006533937, -0.3194557972432347, 0.547884663730034);
  Eigen::Quaterniond q_a_to_c = q_b_to_c * q_a_to_b;
  Eigen::Vector3d v_c = q_a_to_c * v_a;
  std::cout << v_c << std::endl;
```

which gives the same result.

Similarly, converting `erp_b_wrt_a` to a DCM gives the same results as converting `q_a_to_b` to a Matrix3d in Eigen.

This is the export process:

```julia
function export_erp_to_eigen_quat(erp::ERP)
    return [erp.s, -erp.x, -erp.y, -erp.z]
end
```

### Exporting to Rotations.jl

Rotations.jl is an excellent package that implements *active rotations*. When constructing its rotations, you say, "This *does something* to a vector." That is, `RotX(pi/4) * v` rotates `v` by `pi/4` radians about the x axis. This is exactly the opposite of GradientOrientations, where `erp_BA = erpx(pi/4)` says, "frame B is rotated pi/4 radians from frame A", and if `v_A` is a vector expressed in frame A, then `v_B = reframe(erp_BA, v_A)` is that same vector expressed in frame B.

The following give the same results:

```julia
using GradientOrientations
using StaticArrays
erp_b_wrt_a = ERP(x = -0.39015349272073274, y = 0.44885619975114965, z = -0.5331968363460537, s = 0.6016722511245876)
v_a = SA[ -0.09342155668355605, -0.5691321723949833, 0.18864021876469472]
v_b = reframe(erp_b_wrt_a, v_a)
```

This gives `SA[0.5384388474664465, -0.27833775489103457, -0.02891108412057084]`.

```julia
using Rotations
using StaticArrays

# Move the scalar first, and swap the sign of the remaining elements (invert the rotation).
q_a_to_b = QuatRotation(0.6016722511245876, 0.39015349272073274, -0.44885619975114965, 0.5331968363460537)
v_a = SA[ -0.09342155668355605, -0.5691321723949833, 0.18864021876469472] # Same as above.
v_b = q_a_to_b * v_a
```

This gives `SA[0.5384388474664465, -0.27833775489103457, -0.02891108412057084]`, the same as the above.

Extending the examples for composition (quaternion multiplication in Rotations.jl), here's GradientOrientations.jl:

```julia
erp_c_wrt_b = ERP(x = -0.26328146006533937, y = 0.3194557972432347, z = -0.547884663730034, s = 0.7269479084796793)
erp_c_wrt_a = compose(erp_c_wrt_b, erp_b_wrt_a)
v_c = reframe(erp_c_wrt_a, v_a)
```

which gives `SA[0.37890885798068275, 0.28686794099866375, 0.37730493081162575]`.

For Rotations.jl:

```julia
# Move the scalar first, and swap the sign of the remaining elements (invert the rotation).
q_b_to_c = QuatRotation(0.7269479084796793, 0.26328146006533937, -0.3194557972432347, 0.547884663730034)
q_a_to_c = q_b_to_c * q_a_to_b
v_c = q_a_to_c * v_a
```

which also gives `SA[0.37890885798068275, 0.28686794099866375, 0.37730493081162575]`.

Note that `erp2dcm(erp_b_wrt_a)` and `RotMatrix(q_a_to_b)` also give the same matrices.

This is the export process:

```julia
function export_erp_to_rotations_quatrotation(erp::ERP)
    return QuatRotation(erp.s, -erp.x, -erp.y, -erp.z)
end
```
