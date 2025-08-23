export PolygonSymmetries

using StaticArrays

"Reflect `point ∈ C` about an axis in `C`"
reflect(axis, point) = axis * conj(axis' * point) / abs2(axis)

"""Constructive polygon symmetry information

In order to construct a symmetric polygon, it is sufficient to know
the order of rotational symmetry and/or one mirror symmetry axis.
Any two non-identical mirror symmetries generate a rotational symmetry,
therefore one mirror symmetry and one rotational symmetry are enough.
In order to relate the total number of vertices to independent free
parameters, we need to know how many vertices are fixed to be on symmetry
axes.

:rotational_order: Number of discrete rotation steps until a full
    2π rotation is achieved.
:has_mirror_symmetry: If the polygon has a mirror symmetry. If `true`, there are
    `rotational_order` mirror symmetries in total.
:points_on_axes: Number of vertices ∈ {0,1,2} lying on the principal mirror
    symmetry axes.
"""
struct PolygonSymmetries
    rotational_order::Int64
    has_mirror_symmetry::Bool
    points_on_axes::Int64
end

function PolygonSymmetries(poly::Polygon{N}) where {N}
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

    # go through all possible mirror symmetries
    midpoints = (poly.w .+ circshift(poly.w, 1)) / 2
    potential_symmetry_points = sort([poly.w; midpoints], by = angle)

    function is_symmetry_axis(axis)
        poly_reflect = reverse(reflect.(axis, poly.w))
        # check permutations for congruence
        any(all(poly_reflect .≈ circshift(poly.w, i)) for i = 0:N-1)
    end

    idx_1 = findfirst(is_symmetry_axis, potential_symmetry_points)
    has_mirror_symmetry = !isnothing(idx_1)
    points_on_axes = 0
    if has_mirror_symmetry
        axis = potential_symmetry_points[idx_1]
        points_on_axes += count(==(axis), poly.w)

        # search for second mirror symmetry
        idx_2 = findnext(is_symmetry_axis, potential_symmetry_points, idx_1+1)
        if isnothing(idx_2)
            # count all points on the single mirror symmetry axis
            points_on_axes += count(==(-axis), poly.w)
        else
            points_on_axes += count(==(potential_symmetry_points[idx_2]), poly.w)
        end
    end

    PolygonSymmetries(rot_order, has_mirror_symmetry, points_on_axes)
end
