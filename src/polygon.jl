export Polygon

using StaticArrays

abstract type AbstractSymmetry end
struct NoSymmetry <: AbstractSymmetry end
struct CyclicSymmetry{R} <: AbstractSymmetry end
struct BilateralSymmetry{P} <: AbstractSymmetry end
struct DihedralSymmetry{R,P} <: AbstractSymmetry end

struct Polygon{N,S<:AbstractSymmetry,W,B,L}
    w::SVector{N,W}  # vertices
    sym::S  # symmetry
    β::SVector{N,B}  # left-turn angles
    ℓ::SVector{N,L}  # length of edges [i,i+1]
end

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
        DihedralSymmetry{rotational_order,points_on_axes}()
    elseif !has_rotation && has_mirror
        BilateralSymmetry{points_on_axes}()
    elseif has_rotation && !has_mirror
        CyclicSymmetry{rotational_order}()
    else
        NoSymmetry()
    end
end

function Polygon(w::SVector{N,W}, β::SVector{N,B}, ℓ::SVector{N,L}) where {N,W,B,L}
    sym = classify_symmetry(w)
    Polygon{N,typeof(sym),W,B,L}(w, sym, β, ℓ)
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

function Polygon(w::SVector{N,W}, β::SVector{N,B}) where {N,W,B}
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
    Polygon(SVector{N,W}(w), SVector{N,B}(β), SVector{N}(ℓ))
end
