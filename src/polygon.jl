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
function classify_symmetry(w::SVector{N}, β::SVector{N}, ℓ::SVector{N}) where {N}
    # want to preserve angle information on infinite vertices
    # replace all infinite edges with a unique finite length
    my_inf = sum(filter(!isinf, ℓ))  # ∉ ℓ
    ℓ = SVector(replace(x -> isinf(x) ? my_inf : x, ℓ))
    L = MVector(ℓ .* cispi.(β))
    M = MVector(ℓ .* cispi.(circshift(β, -1)) |> reverse)

    rot_equivalent(A, B, rot_a) = all(B .≈ circshift(A, rot_a))

    rot_order = 1 + count(k -> rot_equivalent(L, L, k), 1:N-1)
    reflections = [m for m ∈ 0:N-1 if rot_equivalent(L, M, m)]

    has_rotation = rot_order > 1
    has_mirror = length(reflections) > 0
    if has_rotation && has_mirror
        # With 2 adjacent mirror axes, we can just count how many of them
        # include a vertex.
        points_on_axes = count(r -> iseven(r + N), reflections[1:2])
        α = NaN
        DihedralSymmetry{rot_order,points_on_axes,typeof(α)}(α)
    elseif !has_rotation && has_mirror
        # With 1 mirror axis only, a polygon with an odd number of vertices
        # must have exactly 1 vertex on the axis. Otherwise, there are either
        # 0 or 2 vertices on the axis.
        points_on_axes = isodd(N) ? 1 : (iseven(reflections[1] + N) ? 2 : 0)
        α = NaN
        BilateralSymmetry{points_on_axes,typeof(α)}(α)
    elseif has_rotation && !has_mirror
        CyclicSymmetry{rot_order}()
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
        sym = classify_symmetry(w, β, ℓ)
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
    @assert all(>(1), diff(findall(isinf, w))) "remove consecutive infinities"
    has_2_finite_connected_nodes = false
    k = 0
    # preallocate output
    ℓ = zeros(N)
    for (i, wi) ∈ enumerate(w)
        post₁ = w[mod1(i + 1, N)]
        post₂ = w[mod1(i + 2, N)]
        if !isinf(post₁) && !isinf(post₂)
            has_2_finite_connected_nodes = true
            if isinf(wi)
                # found circshift such that w[N-1] is an infinity and w[1] and w[N] are finite
                k = i + 1
            end
        end
        ℓ[i] = abs(post₁ - wi)
    end
    @assert has_2_finite_connected_nodes "provide at least 2 connected finite vertices"
    if k ≠ 0
        w = circshift(w, -k)
        ℓ = circshift(ℓ, -k)
        β = circshift(β, -k)
    end
    # check supplied angles β for inconsistencies with supplied vertices w
    # these are the inclines of all segments from β
    γ = mod.(angle(w[1] - w[end]) .+ π .* cumsum([0; β[1:end-1]]), 2π)
    # these are the indices of finite segments
    idxs = findall(i -> !isinf(w[i]) && !isinf(w[mod1(i-1,N)]), 1:N)
    # these are the inclines of the finite segments
    α = [angle(w[i] - w[mod1(i-1,N)]) for i ∈ idxs]
    @assert all(γ[idxs] .≈ α) "inconsistent w and β"
    Polygon(SVector{N,W}(w), SVector{N,F}(β), SVector{N,F}(ℓ))
end
