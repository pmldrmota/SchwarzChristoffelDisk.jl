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
:param β: left-turn angles at the nodes of the polygon
:returns: instance of a subtype of AbstractSymmetry
"""
function classify_symmetry(w::SVector{N}, β::SVector{N}) where {N}
    # precompute circshifts
    circshifts_w = SMatrix{N}((w[mod1(i - k, N)] for k ∈ 1:N for i ∈ 1:N)...)
    circshifts_β = SMatrix{N}((β[mod1(i - k, N)] for k ∈ 1:N for i ∈ 1:N)...)

    isinf_or_nan(x) = isinf(x) || isnan(x)
    is_approx_equal(w₁, β₁, w₂, β₂) =
        isinf_or_nan(w₁) && isinf_or_nan(w₂) && mod(β₁, 2) == mod(β₂, 2) || w₁ ≈ w₂

    is_congruent(nodes, turn_angles) = any(
        all(is_approx_equal.(nodes, turn_angles, sw, sβ)) for
        (sw, sβ) ∈ zip(eachcol(circshifts_w), eachcol(circshifts_β))
    )

    is_rotational_symmetry(k) = is_congruent(cispi(2 // k) * w, β)
    # still not correct because equal β at opposite infinities does not mean
    # the lines coming from ∞ are mirror images of each other (there can be an offset).
    # this wouldn't be captured by the finite endpoints if they themselves are on the axis.
    is_mirror_symmetry(α) = is_congruent(reverse(@. α * conj(α' * w) / abs2(α)), reverse(β))

    # classify rotational symmetry
    rotational_order = N + 1 - something(findfirst(is_rotational_symmetry, N:-1:2), N)
    has_rotation = rotational_order > 1

    # candidate mirror symmetry axes: finite points and edge normals
    # if w[i] is its own mirror image, then the axis would be along w[i]
    finite_w = filter(!isinf, w)
    # if w[1] has a mirror image, w[idx], then the axis would be normal to w[1] - w[idx]
    edge_normals = (finite_w[2:end] .- finite_w[1]) * cispi(1 // 2)
    # filter out those that are identical under 180° rotation and eliminate duplicates
    # we can't use `unique()` because we have numerical precision.
    axes = []
    for ax ∈ [finite_w; edge_normals]
        ϕ = angle(ax)
        ax_upper = (ϕ ≥ 0 ? ax : -ax) / abs(ax)
        if isnothing(findfirst(≈(ax_upper), axes))
            push!(axes, ax_upper)
        end
    end
    axes = sort(axes, by = angle)

    # find first symmetry axis, if any
    idx₁ = findfirst(is_mirror_symmetry, axes)
    has_mirror = !isnothing(idx₁)

    on_positive_axis(ax) = p -> ax * conj(p) ≈ abs(p)
    points_on_axes = 0

    if has_mirror
        axis₁ = axes[idx₁]
        points_on_axes += count(on_positive_axis(axis₁), finite_w)
        # look for another axis
        idx₂ = findnext(is_mirror_symmetry, axes, idx₁ + 1)
        # if no other found, then we need to use the negative part of axis₁
        axis₂ = isnothing(idx₂) ? -axis₁ : axes[idx₂]
        points_on_axes += count(on_positive_axis(axis₂), finite_w)
    end

    # are polygons whose finite vertices are symmetric also symmetric in the
    # prevertices even if the overall polygon doesn't have that symmetry?
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
        sym = classify_symmetry(w, β)
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
