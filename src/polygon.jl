export Polygon, first_independent_vertex, num_independent_vertices

using StaticArrays

struct Polygon{N,S<:AbstractSymmetry,W,F}
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
        new{N,S,W,F}(w, s, β, ℓ)
    end
end

Polygon(w, β, ℓ) = Polygon(w, classify_symmetry(w, β, ℓ), β, ℓ)

function calc_β_ℓ(w::SVector{N,W}) where {N,W}
    @assert count(isinf, w) == 0 "must specify angles if there are infinities"
    # preallocate output
    β = Vector{Float64}(undef, N)
    ℓ = Vector{Float64}(undef, N)
    for (i, wi) ∈ enumerate(w)
        post = w[mod1(i + 1, N)] - wi
        pre = wi - w[mod1(i - 1, N)]
        β[i] = angle(pre' * post) / π
        ℓ[i] = abs(post)
    end
    (SVector{N}(β), SVector{N}(ℓ))
end

Polygon(w::SVector{N,W}) where {N,W} = Polygon(w, calc_β_ℓ(w)...)

function Polygon(w::SVector{N,W}, β::SVector{N,F}) where {N,W,F}
    ℓ = [abs(w[mod1(i + 1, N)] - w[i]) for i ∈ 1:N]
    # todo: check supplied angles β for inconsistencies with supplied vertices w
    Polygon(SVector{N,W}(w), SVector{N,F}(β), SVector{N,F}(ℓ))
end

function Polygon(w_base::SVector{B,W}, s::CyclicSymmetry{R}) where {B,W,R}
    N = B * R
    w = SVector{N,W}(ntuple(i -> w_base[mod1(i, B)], N))
    Polygon(w, s, calc_β_ℓ(w)...)
end

function reflect(axis, point)
    normalised_axis = axis / abs(axis)
    normalised_axis * conj(normalised_axis' * point)
end

make_mirror(w::SVector{B,W}, s::BilateralSymmetry{0}) where {B,W} =
    SVector{2B,W}(w..., reflect.(s.axis, w)[end:-1:1])

make_mirror(w::SVector{B,W}, s::BilateralSymmetry{2}) where {B,W} =
    SVector{2B-2,W}(w..., reflect.(s.axis, w)[end-1:-1:2])

function make_mirror(w::SVector{B,W}, s::BilateralSymmetry{1}) where {B,W}
    w_r = reflect.(s.axis, w)
    if w_r[1] ≈ w[1]
        SVector{2B-1,W}(w..., w_r[end:-1:2]...)
    else
        SVector{2B-1,W}(w..., w_r[end-1:-1:1]...)
    end
end

function Polygon(w_base, symmetry::BilateralSymmetry)
    w = make_mirror(w_base, symmetry)
    Polygon(w, symmetry, calc_β_ℓ(w)...)
end

function Polygon(w_base::SVector{B,W}, symmetry::DihedralSymmetry{R,P}) where {B,W,R,P}
    # relies on assumption that axis includes a vertex if P = 1.
    w_rotbase = make_mirror(w_base, BilateralSymmetry{P}(symmetry.axis))
    if P == 2
        w_rotbase = SVector{2B-3,W}(w_rotbase[1:end-1]...)
    end
    Polygon(w_rotbase, CyclicSymmetry{R}())
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
function first_independent_vertex(poly::Polygon{<:Any,<:DihedralSymmetry{<:Any,P}}) where {P}
    p = Polygon(poly.w, BilateralSymmetry{P}(poly.s.axis), poly.β, poly.ℓ)
    first_independent_vertex(p)
end
