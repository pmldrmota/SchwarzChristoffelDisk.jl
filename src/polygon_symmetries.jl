export PolygonSymmetries

using StaticArrays

"Reflect `point ∈ C` about an axis in `C`"
function reflect(axis, point)
    normalised_axis = axis / abs(axis)
    normalised_axis * conj(normalised_axis' * point)
end

"""Constructive polygon symmetry information

In order to construct a symmetric polygon, it is sufficient to know
the order of rotational symmetry and/or one mirror symmetry axis.
Any two non-identical mirror symmetries generate a rotational symmetry,
therefore one mirror symmetry and one rotational symmetry are enough.
In order to relate the total number of vertices to independent free
parameters, we need to know how many vertices are fixed to be on symmetry
axes.

:rotational_symmetry: Number of discrete rotation steps until a full
    2π rotation is achieved.
:has_mirror_symmetry: If the polygon has a mirror symmetry. If `true`, there are
    `rotational_symmetry` mirror symmetries in total.
:points_on_axes: Number of vertices ∈ {0,1,2} lying on the principal mirror
    symmetry axes.
"""
struct PolygonSymmetries
    rotational_symmetry::Int64
    has_mirror_symmetry::Bool
    points_on_axes::Int64
end

function PolygonSymmetries(poly::Polygon{N,W,B,M}) where {N,W,B,M}
    # go through all possible rotational symmetries
    rot_order = 1
    for k = N:-1:2
        # rotation by 2π/k
        poly_rot = cispi(2 / k) * poly.w
        # check permutations for congruence
        if any(all(poly_rot .≈ circshift(poly.w, i)) for i = 0:N-1)
            rot_order = k
            break
        end
    end

    function is_symmetry_axis(axis)
        poly_reflect = reverse(reflect.(axis, poly.w))
        # check permutations for congruence
        any(all(poly_reflect .≈ circshift(poly.w, i)) for i = 0:N-1)
    end

    midpoints = (poly.w .+ circshift(poly.w, 1)) / 2
    potential_symmetry_points = sort([poly.w; midpoints], by = angle)
    idx_1 = findfirst(is_symmetry_axis, potential_symmetry_points)
    has_mirror_symmetry = !isnothing(idx_1)
    axis_1 = potential_symmetry_points[idx_1]

    points_on_axes = 0
    if has_mirror_symmetry
        points_on_axes += count(≈(axis_1), poly.w)
    end
    if rot_order > 1
        # there must be a second mirror symmetry
        idx_2 = findnext(is_symmetry_axis, potential_symmetry_points, idx_1+1)
        @assert !isnothing(idx_2) "Expected a second mirror symmetry from rotation"
        axis_2 = potential_symmetry_points[idx_2]
        points_on_axes += count(≈(axis_2), poly.w)
    else
        # count all points on the single mirror symmetry axis
        points_on_axes += count(≈(-axis_1), poly.w)
    end

    PolygonSymmetries(rot_order, has_mirror_symmetry, points_on_axes)
end
