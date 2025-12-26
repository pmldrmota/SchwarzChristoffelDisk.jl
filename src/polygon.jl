export Polygon, first_independent_vertex, num_independent_vertices

using StaticArrays

struct Polygon{N,S<:AbstractSymmetry,W<:Complex,F}
    w::SVector{N,W}  # vertices w.r.t. centre of mass
    s::S
    β::SVector{N,F}  # left-turn angles
    ℓ::SVector{N,F}  # length of edges [i,i+1]

    function Polygon(
        w::SVector{N,W},
        s::S,
        β::SVector{N,F},
        ℓ::SVector{N,F},
    ) where {N,S,W,F}
        N < 3 && throw(ArgumentError("Polygon must have at least 3 nodes"))
        double_∞ = findfirst(i -> isinf(w[i]) && isinf(w[mod1(i + 1, N)]), 1:N)
        isnothing(double_∞) || throw(ArgumentError("remove consecutive infinities"))
        ∑β = sum(β)
        ∑β ≈ 2 || throw(ArgumentError("wrong angles (∑β=$∑β)"))
        new{N,S,W,F}(complex.(w), s, β, ℓ)
    end
end

Polygon(w::SVector{N}, β::SVector{N,F}, ℓ::SVector{N,F}) where {N,F} =
    Polygon(complex.(w), classify_symmetry(w, β, ℓ), β, ℓ)

function calc_β_ℓ(w::SVector{N,<:Complex}, β_lu::Dict{Int,<:Number}) where {N}
    # preallocate output
    β = Vector{Float64}(undef, N)
    ℓ = Vector{Float64}(undef, N)
    for (i, wi) ∈ enumerate(w)
        post = w[mod1(i + 1, N)] - wi
        pre = wi - w[mod1(i - 1, N)]
        β[i] = if isfinite(pre) && isfinite(post)
            angle(pre' * post) / π
        elseif i ∈ keys(β_lu)
            β_lu[i]
        elseif isinf(wi)
            NaN
        else
            throw("require β[$i] next to infinity")
        end
        ℓ[i] = abs(post)
    end
    # Infer β at infinity
    count_nan = count(isnan, β)
    if count_nan != 0
        if count_nan == 1
            idx = findfirst(isnan, β)
            β[idx] = 2 - sum(filter(!isnan, β))
        else
            throw("Cannot infer β at multiple infinities")
        end
    end

    (SVector{N}(β), SVector{N}(ℓ))
end

function Polygon(w::SVector{N}, β_lu::Dict{Int,<:Number} = Dict{Int,Float64}()) where {N}
    w = complex.(w)
    Polygon(SVector{N}(w), calc_β_ℓ(w, β_lu)...)
end

function Polygon(
    w_base::SVector{B},
    s::CyclicSymmetry{R},
    β_lu_base::Dict{Int,<:Number} = Dict{Int,Float64}(),
) where {B,R}
    N = B * R
    w = Vector{ComplexF64}(undef, N)
    for i ∈ 0:(R-1)
        r = cispi(2i // R)
        w[(1+B*i):(B*(i+1))] .= r .* w_base
    end
    w = SVector{N}(w)
    β_lu = Dict{Int,Float64}()
    for (k, β) ∈ pairs(β_lu_base)
        for i ∈ 0:(R-1)
            β_lu[B*i+k] = β
        end
    end
    Polygon(w, s, calc_β_ℓ(w, β_lu)...)
end

function reflect(axis, point)
    normalised_axis = axis / abs(axis)
    normalised_axis * conj(normalised_axis' * point)
end

"""
Assumes that the mirror image continues in the same order as the base,
i.e., the last vertex in the base is connected to the first vertex in the
image.
"""
make_mirror(w::SVector{B}, s::BilateralSymmetry{0}) where {B} =
    SVector{2B}(w..., reflect.(s.axis, w)[end:-1:1]...)
make_mirror(w::SVector{B}, s::BilateralSymmetry{1}) where {B} =
    SVector{2B-1}(w..., reflect.(s.axis, w)[(end-1):-1:1]...)
make_mirror(w::SVector{B}, s::BilateralSymmetry{2}) where {B} =
    SVector{2B-2}(w..., reflect.(s.axis, w)[(end-1):-1:2]...)

function mirror_β_lu(K, β_lu_base)
    β_lu = Dict{Int,Float64}()
    for (k, β) ∈ pairs(β_lu_base)
        β_lu[k] = β
        β_lu[K-k] = β
    end
    β_lu
end

function Polygon(
    w_base,
    symmetry::BilateralSymmetry{P},
    β_lu_base::Dict{Int,<:Number} = Dict{Int,Float64}(),
) where {P}
    w = make_mirror(w_base, symmetry)
    β_lu = mirror_β_lu(length(w) + min(1, P), β_lu_base)
    Polygon(w, symmetry, calc_β_ℓ(w, β_lu)...)
end

function Polygon(
    w_base::SVector{B,W},
    symmetry::DihedralSymmetry{R,P},
    β_lu_base::Dict{Int,<:Number} = Dict{Int,Float64}(),
) where {B,W,R,P}
    # relies on assumption that axis includes a vertex if P = 1.
    w_rotbase = make_mirror(w_base, BilateralSymmetry{P}(symmetry.axis))
    if P == 2
        w_rotbase = SVector{2B-P-1,W}(w_rotbase[1:(end-1)]...)
    end
    β_lu_rotbase = mirror_β_lu(length(w_rotbase) + 1, β_lu_base)
    temp = Polygon(w_rotbase, CyclicSymmetry{R}(), β_lu_rotbase)
    Polygon(temp.w, symmetry, temp.β, temp.ℓ)
end

num_independent_vertices(::Polygon{N,NoSymmetry}) where {N} = N
num_independent_vertices(::Polygon{N,CyclicSymmetry{R}}) where {N,R} = N ÷ R
num_independent_vertices(::Polygon{N,<:BilateralSymmetry{P}}) where {N,P} = P + (N - P) ÷ 2
num_independent_vertices(::Polygon{N,<:DihedralSymmetry{R,P}}) where {N,R,P} =
    P + (N ÷ R - P) ÷ 2

"Symmetry start index"
first_independent_vertex(::Polygon{<:Any,NoSymmetry}) = 1
first_independent_vertex(::Polygon{<:Any,<:CyclicSymmetry}) = 1

function is_left(axis, point)
    u = imag(conj(axis) * point)
    u < 0 && abs(u) > abs(axis) * √eps()
end

""" Bilateral symmetry start index

For P=1 and P=2 we search for patterns
     l.r
     l∞r
    l∞.∞r
"""
first_independent_vertex(poly::Polygon{N,<:BilateralSymmetry}) where {N} =
    findfirst(i -> begin
        w = poly.w[i]
        if isinf(w)
            w₋ = poly.w[mod1(i - 1, N)]
            w₊ = poly.w[mod1(i + 1, N)]
            is_left(poly.s.axis, w₋) != is_left(poly.s.axis, w₊)
        else
            v = poly.s.axis / abs(poly.s.axis) * abs(w)
            v ≈ w || v ≈ -w
        end
    end, 1:N)

""" Bilateral symmetry start index

For P=0 we search for the finite vertices whose edge is divided by the axis.
pattern: lr
"""
first_independent_vertex(poly::Polygon{N,<:BilateralSymmetry{0}}) where {N} = findfirst(
    i -> begin
        w₋ = poly.w[mod1(i - 1, N)]
        w = poly.w[i]
        isfinite(w₋) &&
            isfinite(w) &&
            (is_left(poly.s.axis, w₋) != is_left(poly.s.axis, w))
    end,
    1:N,
)

""" Dihedral symmetry start index

see BilateralSymmetry
"""
function first_independent_vertex(
    poly::Polygon{<:Any,<:DihedralSymmetry{<:Any,P}},
) where {P}
    p = Polygon(poly.w, BilateralSymmetry{P}(poly.s.axis), poly.β, poly.ℓ)
    first_independent_vertex(p)
end
