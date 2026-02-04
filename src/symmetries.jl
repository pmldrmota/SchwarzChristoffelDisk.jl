export NoSymmetry,
    CyclicSymmetry,
    BilateralSymmetry,
    DihedralSymmetry,
    PolygonSymmetry,
    classify_symmetry,
    get_axis

using StaticArrays
import Base

abstract type AbstractSymmetry end

"""Rotational symetry

Type parameter `{R}` encodes the number of steps making up a full 2π rotation.
"""
struct CyclicSymmetry{R} <: AbstractSymmetry end

"""No symmetry (alias for CyclicSymmetry with `R=1`)
"""
const NoSymmetry = CyclicSymmetry{1}

"""Regular polygon symmetry

Type parameter `{R}` encodes the number of steps making up a full 2π rotation.
Type parameter `{P}` encodes the number of polygon nodes on the symmetry axis.
The field `axis` stores one axis of symmetry. The other ones can be constructed
by rotation (see `symmetry_axes`).
If `P > 0`, then `axis` is assumed to coincide with a vertex.
"""
struct DihedralSymmetry{R,P,T<:Complex} <: AbstractSymmetry
    axis::T

    function DihedralSymmetry{R,P}(axis) where {R,P}
        axis = complex(axis)
        abs(axis) < √eps() && throw(ArgumentError("axis must be nonzero"))
        new{R,P,typeof(axis)}(axis)
    end
end

"""Mirror symetry (alias for DihedralSymmetry with `R=1`)

Type parameter `{P}` encodes the number of polygon nodes on the symmetry axis.
The field `axis` stores the axis of symmetry.
If `P > 0`, then `axis` is assumed to coincide with a vertex.
"""
const BilateralSymmetry{P,T} = DihedralSymmetry{1,P,T}
BilateralSymmetry{P}(axis) where {P} = DihedralSymmetry{1,P}(axis)

"""Wrapper for symmetry of a Polygon

The extra element `first_independent_vertex` stores the index
"""
struct PolygonSymmetry{S<:AbstractSymmetry}
    symmetry::S
    first_independent_vertex::Int
end

"Getter functions"
get_axis(s::DihedralSymmetry) = s.axis
get_axis(ps::PolygonSymmetry{<:DihedralSymmetry}) = ps.symmetry.axis

"Helper function to determine on which side of the axis a point is"
function which_side(axis, point)
    u = imag(axis' * point)
    abs(u) < abs(axis) * √eps() ? 0 : sign(u)
end
is_on(axis, point) = which_side(axis, point) |> iszero

"Equality check taking into account symmetric transformations"
Base.:(==)(a::DihedralSymmetry{R,P}, b::DihedralSymmetry{R,P}) where {R,P} =
    any(i -> is_on(a.axis, cispi(i // R) * b.axis), 1:R)

"List of all symmetry axes for a DihedralSymmetry"
symmetry_axes(sym::DihedralSymmetry{R}) where {R} =
    SVector{R}(cispi(k // R) * sym.axis for k ∈ 0:(R-1))
symmetry_axes(psym::PolygonSymmetry) = symmetry_axes(psym.symmetry)

"Number of independent vertices of N-vertex polygon of typed symmetry"
num_independent_vertices(N, ::CyclicSymmetry{R}) where {R} = N ÷ R
num_independent_vertices(N, ::DihedralSymmetry{R,P}) where {R,P} = P + (N ÷ R - P) ÷ 2

"""Classify the symmetry group of a polygon

:param w: vertices of the polygon
:param β: left-turn angles at the nodes of the polygon
:param ℓ: edge lengths of the polygon

:returns: instance of PolygonSymmetry
"""
function classify_symmetry(
    w::StaticVector{N},
    β::StaticVector{N},
    ℓ::StaticVector{N},
) where {N}
    # want to preserve angle information on infinite vertices
    # replace all infinite edges with a unique finite length
    my_inf = 1 + sum(filter(isfinite, ℓ))  # ∉ ℓ
    ℓ = [
        if isinf(ℓ[i])
            my_inf + abs2(isinf(w[i]) ? w[mod1(i+1,N)] : w[i])
        else
            ℓ[i]
        end for i ∈ eachindex(ℓ)
    ]
    L = SVector(ℓ .* cispi.(β))
    M = SVector(ℓ .* cispi.(circshift(β, -1)) |> reverse)

    congruent_circshift(A, B, k) = all(B[i] ≈ A[mod1(i - k, N)] for i ∈ 1:N)

    rot_order = 1 + count(k -> congruent_circshift(L, L, k), 1:(N-1))
    refl₁ = findfirst(m -> congruent_circshift(L, M, m), 1:N)

    has_rotation = rot_order > 1
    has_mirror = !isnothing(refl₁)
    if has_mirror
        wc(i) = w[mod1(i, N)]
        circular_average(a, b) = sqrt(a * b / (abs(a) * abs(b)))
        first_indep(refl) = mod1(1 - cld(refl - 1, 2), N)

        axis(refl) =
            if iseven(refl)
                Δ = refl ÷ 2
                v = wc(1 - Δ)
                isfinite(v) ? v : circular_average(wc(0 - Δ), wc(2 - Δ))
            else
                Δ₊ = (refl + 1) ÷ 2
                Δ₋ = (refl - 1) ÷ 2
                circular_average(wc(1 - Δ₊), wc(1 - Δ₋))
            end

        if has_rotation
            # With 2 adjacent mirror axes, we can just count how many of them
            # include a vertex.
            refl₂ = findnext(m -> congruent_circshift(L, M, m), 1:N, refl₁ + 1)
            points_on_axes = iseven(refl₁) + iseven(refl₂)
            # report the axis which contains the point
            refl = iseven(refl₁) ? refl₁ : refl₂
            s = DihedralSymmetry{rot_order,points_on_axes}(axis(refl))
            PolygonSymmetry(s, first_indep(refl))
        else
            # With 1 mirror axis only, a polygon with an odd number of vertices
            # must have exactly 1 vertex on the axis. Otherwise, there are either
            # 0 or 2 vertices on the axis.
            points_on_axes = isodd(N) ? 1 : (iseven(refl₁) ? 2 : 0)
            s = BilateralSymmetry{points_on_axes}(axis(refl₁))
            # report the index which is on the axis
            refl = isodd(N) && isodd(refl₁) ? refl₁ + N : refl₁
            PolygonSymmetry(s, first_indep(refl))
        end
    else
        s = CyclicSymmetry{rot_order}()
        PolygonSymmetry(s, 1)
    end
end
