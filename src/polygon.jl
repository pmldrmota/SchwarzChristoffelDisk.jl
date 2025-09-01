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
The field `axis` stores the axis of symmetry.
"""
struct BilateralSymmetry{P,T} <: AbstractSymmetry
    axis::Complex{T}
    BilateralSymmetry{P}(axis::Complex{T}) where {P,T} = new{P,T}(axis)
end

"""Regular polygon symmetry

Type parameter `{R}` encodes the number of steps making up a full 2π rotation.
Type parameter `{P}` encodes the number of polygon nodes on the symmetry axis.
The field `axis` stores one axis of symmetry. The other ones can be constructed
by rotation.
"""
struct DihedralSymmetry{R,P,T} <: AbstractSymmetry
    axis::Complex{T}
    DihedralSymmetry{R,P}(axis::Complex{T}) where {R,P,T} = new{R,P,T}(axis)
end

"""Classify the symmetry group of a polygon

:param w: vertices of the polygon
:param β: left-turn angles at the nodes of the polygon
:returns: instance of a subtype of AbstractSymmetry
"""
function classify_symmetry(w::SVector{N}, β::SVector{N}, ℓ::SVector{N}) where {N}
    # want to preserve angle information on infinite vertices
    # replace all infinite edges with a unique finite length
    my_inf = sum(filter(isfinite, ℓ))  # ∉ ℓ
    ℓ = SVector(replace(x -> isinf(x) ? my_inf : x, ℓ))
    L = SVector(ℓ .* cispi.(β))
    M = SVector(ℓ .* cispi.(circshift(β, -1)) |> reverse)

    congruent_circshift(A, B, k) = all(B[i] ≈ A[mod1(i - k, N)] for i ∈ 1:N)

    rot_order = 1 + count(k -> congruent_circshift(L, L, k), 1:N-1)
    refl₁ = findfirst(m -> congruent_circshift(L, M, m), 1:N)

    has_rotation = rot_order > 1
    has_mirror = !isnothing(refl₁)
    if has_mirror
        wc(i) = w[mod1(i, N)]
        axis = if iseven(refl₁)
            Δ = refl₁ ÷ 2
            v = wc(1 - Δ)
            isfinite(v) ? v : wc(0 - Δ) + wc(2 - Δ)
        else
            Δ₊ = (refl₁ + 1) ÷ 2
            Δ₋ = (refl₁ - 1) ÷ 2
            wc(1 - Δ₊) + wc(1 - Δ₋)
        end

        if has_rotation
            # With 2 adjacent mirror axes, we can just count how many of them
            # include a vertex.
            refl₂ = findnext(m -> congruent_circshift(L, M, m), 1:N, refl₁ + 1)
            points_on_axes = iseven(refl₁) + iseven(refl₂)
            DihedralSymmetry{rot_order,points_on_axes}(axis)
        else
            # With 1 mirror axis only, a polygon with an odd number of vertices
            # must have exactly 1 vertex on the axis. Otherwise, there are either
            # 0 or 2 vertices on the axis.
            points_on_axes = isodd(N) ? 1 : (iseven(refl₁) ? 2 : 0)
            BilateralSymmetry{points_on_axes}(axis)
        end
    else
        if has_rotation
            CyclicSymmetry{rot_order}()
        else
            NoSymmetry()
        end
    end
end

struct Polygon{N,S<:AbstractSymmetry,W,F}
    w::SVector{N,W}  # vertices
    s::S
    β::SVector{N,F}  # left-turn angles
    ℓ::SVector{N,F}  # length of edges [i,i+1]

    function Polygon(w::SVector{N,W}, β::SVector{N,F}, ℓ::SVector{N,F}) where {N,W,F}
        sym = classify_symmetry(w, β, ℓ)
        new{N,typeof(sym),W,F}(w, sym, β, ℓ)
    end
end

function Polygon(w::SVector{N,W}) where {N,W}
    @assert length(w) > 2 "Polygon must have at least 3 nodes"
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
    @assert length(w) > 2 "Polygon must have at least 3 nodes"
    ∑β = sum(β)
    @assert ∑β ≈ 2 "wrong angles (∑β=$∑β)"
    @assert all(>(1), diff(findall(isinf, w))) "remove consecutive infinities"

    wc(i) = w[mod1(i, N)]

    ℓ = [abs(wc(i + 1) - w[i]) for i ∈ eachindex(w)]

    # find circshift such that w[N-1] is an infinity and w[N] and w[1] are finite
    k = findfirst(i -> isinf(wc(i - 2)) && isfinite(wc(i - 1)) && isfinite(wc(i)), 1:N)
    if isnothing(k)
        # If there is no sequence [Inf, Finite, Finite] to be found, then there
        # is either no infinity or not one connected section consisting of at least
        # 2 vertices - which is problematic as that section would have a rotational
        # degree of freedom.
        @assert count(isinf, w) == 0 "provide at least 2 connected finite vertices"
    else
        w = circshift(w, k)
        ℓ = circshift(ℓ, k)
        β = circshift(β, k)
    end

    # check supplied angles β for inconsistencies with supplied vertices w
    # these are the inclines of all segments from β
    γ = mod.(angle(w[1] - w[end]) .+ π .* cumsum([0; β[1:end-1]]), 2π)
    # these are the indices of finite segments
    idxs = findall(i -> isfinite(w[i]) && isfinite(wc(i - 1)), 1:N)
    # these are the inclines of the finite segments
    α = [angle(w[i] - wc(i - 1)) for i ∈ idxs]
    @assert all(γ[idxs] .≈ α) "inconsistent w and β"

    Polygon(SVector{N,W}(w), SVector{N,F}(β), SVector{N,F}(ℓ))
end
