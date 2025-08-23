export Polygon

using StaticArrays

abstract type AbstractSymmetry end

"""No symmetry"""
struct NoSymmetry <: AbstractSymmetry end
"""Rotational symetry

Type parameter `{R}` encodes the number of steps making up a full 2π rotation.
"""
struct CyclicSymmetry{R} <: AbstractSymmetry end

"""Mirror symetry

Type parameter `{P}` encodes the number of polygon nodes on the symmetry axis.
The field `α` stores the angle of the axis with respect to the real axis.
"""
struct BilateralSymmetry{P,A} <: AbstractSymmetry
    α::A
end

"""Regular polygon symmetry

Type parameter `{R}` encodes the number of steps making up a full 2π rotation.
Type parameter `{P}` encodes the number of polygon nodes on the symmetry axis.
The field `α` stores the smallest out of all angles that any of the mirror symmetry axes
subtend with the real axis.
"""
struct DihedralSymmetry{R,P,A} <: AbstractSymmetry
    α::A
end

"""Classify the symmetry group of a polygon

:param w: vertices of the polygon
:returns: instance of a subtype of AbstractSymmetry
"""
function classify_symmetry(w::SVector{N}) where {N}
    # todo: is_congruent can be optimised to avoid allocations in circshift.
    is_congruent(nodes) = any(all(nodes .≈ circshift(w, i)) for i ∈ 0:N-1)
    is_rotational_symmetry(k) = is_congruent(cispi(2 / k) * w)
    is_mirror_symmetry(α) = is_congruent(reverse(@. α * conj(α' * w) / abs2(α)))

    # classify rotational symmetry
    rotational_order = N + 1 - something(findfirst(is_rotational_symmetry, N:-1:2), N)
    has_rotation = rotational_order > 1

    # candidate mirror symmetry axes = vertices ∪ edge midpoints
    midpoints = @. (w + w[mod1(2:N+1, N)]) / 2
    axes = sort([w; midpoints], by = angle)

    # find first symmetry axis, if any
    idx₁ = findfirst(is_mirror_symmetry, axes)
    has_mirror = !isnothing(idx₁)
    points_on_axes = 0

    if has_mirror
        axis₁ = axes[idx₁]
        points_on_axes += count(==(axis₁), w)

        # look for another axis
        idx₂ = findnext(is_mirror_symmetry, axes, idx₁ + 1)
        axis₂ = isnothing(idx₂) ? -axis₁ : axes[idx₂]
        points_on_axes += count(==(axis₂), w)
    end

    if has_rotation && has_mirror
        α = angle(axes[idx₁])
        DihedralSymmetry{rotational_order,points_on_axes,typeof(α)}(α)
    elseif !has_rotation && has_mirror
        α = angle(axes[idx₁])
        BilateralSymmetry{points_on_axes,typeof(α)}(α)
    elseif has_rotation && !has_mirror
        CyclicSymmetry{rotational_order}()
    else
        NoSymmetry()
    end
end

struct Polygon{N,S<:AbstractSymmetry,W,F}
    w::SVector{N,W}  # vertices
    s::S
    β::SVector{N,F}  # left-turn angles
    ℓ::SVector{N,F}  # length of edges [i,i+1]

    function Polygon(w::SVector{N,W}, β::SVector{N,F}, ℓ::SVector{N,F}) where {N,W,F}
        sym = classify_symmetry(w)
        new{N,typeof(sym),W,F}(w, sym, β, ℓ)
    end
end

function Polygon(w::SVector{N,W}) where {N,W}
    @assert count(isinf, w) == 0 "must specify angles if there are infinities"
    # preallocate output
    β = zeros(N)
    ℓ = zeros(N)
    for (i, wi) ∈ enumerate(w)
        post = w[mod1(i + 1, N)] - wi
        pre = wi - w[mod1(i - 1, N)]
        β[i] = angle(pre' * post) / π
        ℓ[i] = abs(post)
    end
    Polygon(w, SVector{N}(β), SVector{N}(ℓ))
end

function Polygon(w::SVector{N,W}, β::SVector{N,F}) where {N,W,F}
    ∑β = sum(β)
    @assert ∑β ≈ 2 "wrong angles (∑β=$∑β)"
    k = 0
    # preallocate output
    ℓ = zeros(N)
    for (i, wi) ∈ enumerate(w)
        post = mod1(i + 1, N)
        @assert !(isinf(wi) && isinf(w[post])) "remove consecutive infinity at k=$post"
        if isinf(w[mod1(i - 2, N)]) && !isinf(w[mod1(i - 1, N)]) && !isinf(wi)
            # find circshift such that w[N-1] is an infinity and w[1] and w[N] are finite
            k = i - 1
        end
        ℓ[i] = abs(w[post] - wi)
    end
    if k != 0
        w = circshift(w, -k)
        ℓ = circshift(ℓ, -k)
        β = circshift(β, -k)
    end
    Polygon(SVector{N,W}(w), SVector{N,F}(β), SVector{N,F}(ℓ))
end
