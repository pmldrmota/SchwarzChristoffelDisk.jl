export NoSymmetry, CyclicSymmetry, BilateralSymmetry, DihedralSymmetry, classify_symmetry

using StaticArrays
import Base

abstract type AbstractSymmetry end

"""No symmetry"""
struct NoSymmetry <: AbstractSymmetry end

"""Rotational symetry

Type parameter `{R}` encodes the number of steps making up a full 2π rotation.
"""
struct CyclicSymmetry{R} <: AbstractSymmetry end

"""Mirror symetry

Type parameter `{P}` encodes the number of polygon nodes on the symmetry axis.
The field `axis` stores the axis of symmetry.
If `P > 0`, then `axis` is assumed to coincide with a vertex.
"""
struct BilateralSymmetry{P,T} <: AbstractSymmetry
    axis::Complex{T}

    function BilateralSymmetry{P}(axis::Complex{T}) where {P,T}
        abs(axis) < √eps() && throw(ArgumentError("axis must be nonzero"))
        new{P,T}(axis)
    end
end

"""Regular polygon symmetry

Type parameter `{R}` encodes the number of steps making up a full 2π rotation.
Type parameter `{P}` encodes the number of polygon nodes on the symmetry axis.
The field `axis` stores one axis of symmetry. The other ones can be constructed
by rotation.
If `P > 0`, then `axis` is assumed to coincide with a vertex.
"""
struct DihedralSymmetry{R,P,T} <: AbstractSymmetry
    axis::Complex{T}

    function DihedralSymmetry{R,P}(axis::Complex{T}) where {R,P,T}
        abs(axis) < √eps() && throw(ArgumentError("axis must be nonzero"))
        new{R,P,T}(axis)
    end
end

is_on(axis, point) = abs(imag(axis' * point)) < abs(axis) * √eps()

Base.:(==)(a::BilateralSymmetry{P}, b::BilateralSymmetry{P}) where {P} =
    is_on(a.axis, b.axis)
Base.:(==)(a::DihedralSymmetry{R,P}, b::DihedralSymmetry{R,P}) where {R,P} =
    any(i -> is_on(a.axis, cispi(2i // R) * b.axis), 1:R)

"""Classify the symmetry group of a polygon

:param w: vertices of the polygon
:param β: left-turn angles at the nodes of the polygon
:returns: instance of a subtype of AbstractSymmetry
"""
function classify_symmetry(w::StaticVector{N}, β::StaticVector{N}, ℓ::StaticVector{N}) where {N}
    # want to preserve angle information on infinite vertices
    # replace all infinite edges with a unique finite length
    my_inf = 1 + sum(filter(isfinite, ℓ))  # ∉ ℓ
    ℓ = SVector(replace(x -> isinf(x) ? my_inf : x, ℓ))
    L = SVector(ℓ .* cispi.(β))
    M = SVector(ℓ .* cispi.(circshift(β, -1)) |> reverse)

    congruent_circshift(A, B, k) = all(B[i] ≈ A[mod1(i - k, N)] for i ∈ 1:N)

    rot_order = 1 + count(k -> congruent_circshift(L, L, k), 1:(N-1))
    refl₁ = findfirst(m -> congruent_circshift(L, M, m), 1:N)

    has_rotation = rot_order > 1
    has_mirror = !isnothing(refl₁)
    if has_mirror
        wc(i) = w[mod1(i, N)]
        circular_average(a, b) = cis((angle(a) + angle(b)) / 2)

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
            DihedralSymmetry{rot_order,points_on_axes}(axis(iseven(refl₁) ? refl₁ : refl₂))
        else
            # With 1 mirror axis only, a polygon with an odd number of vertices
            # must have exactly 1 vertex on the axis. Otherwise, there are either
            # 0 or 2 vertices on the axis.
            points_on_axes = isodd(N) ? 1 : (iseven(refl₁) ? 2 : 0)
            BilateralSymmetry{points_on_axes}(axis(refl₁))
        end
    else
        if has_rotation
            CyclicSymmetry{rot_order}()
        else
            NoSymmetry()
        end
    end
end
