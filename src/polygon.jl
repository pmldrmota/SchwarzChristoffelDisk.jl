export Polygon, first_independent_vertex, num_independent_vertices

using StaticArrays

struct Polygon{N,S<:AbstractSymmetry,W<:Complex,F,G}
    w::SVector{N,W}  # vertices w.r.t. centre of mass
    s::S
    β::SVector{N,F}  # left-turn angles
    ℓ::SVector{N,G}  # length of edges [i,i+1]

    function Polygon(
        w::StaticVector{N},
        s::S,
        β::StaticVector{N},
        ℓ::StaticVector{N},
    ) where {N,S}
        N < 3 && throw(ArgumentError("Polygon must have at least 3 nodes"))
        double_∞ = findfirst(i -> isinf(w[i]) && isinf(w[mod1(i + 1, N)]), 1:N)
        isnothing(double_∞) || throw(ArgumentError("found consecutive infinities"))
        for (i, (wi, βi)) ∈ enumerate(zip(w, β))
            if isfinite(wi)
                (-1 ≤ βi < 1) || throw(ArgumentError("β[$i] ∉ [-1, 1) at node $i ($wi)"))
            else
                (1 ≤ βi ≤ 3) || throw(ArgumentError("β[$i] ∉ [1, 3] at node $i (∞)"))
            end
        end
        ∑β = sum(β)
        isapprox(∑β, 2; rtol = 1e-5) || throw(ArgumentError("wrong angles (∑β=$∑β)"))
        new{N,S,eltype(w),eltype(β),eltype(ℓ)}(SVector(w), s, SVector(β), SVector(ℓ))
    end
end

"""Calculate left-turn angles and edge lengths from vertices `w`

The returned β may contain `NaN`s.
"""
function calc_β_ℓ(w::StaticVector{N}, β_lu::Dict{Int,<:Number}) where {N}
    # preallocate output
    β = MVector{N,Float64}(undef)
    ℓ = MVector{N,Float64}(undef)
    for (i, wi) ∈ enumerate(w)
        post = w[mod1(i + 1, N)] - wi
        pre = wi - w[mod1(i - 1, N)]
        β[i] = if isfinite(pre) && isfinite(post)
            α = angle(pre' * post) / π
            if haskey(β_lu, i) && !isapprox(β_lu[i], α; rtol = 1e-4, atol = 1e-5)
                @warn "β_lu[$i]=$(β_lu[i]) inconsistent with α=$α; using α"
            end
            α
        elseif haskey(β_lu, i)
            β_lu[i]
        else
            NaN
        end
        ℓ[i] = abs(post)
    end
    (β, ℓ)
end

"Construct Polygon and infer symmetry from nodes and angles"
function Polygon(w::SVector{N}, β_lu::Dict{Int,<:Number} = Dict{Int,Float64}()) where {N}
    w = complex.(w)
    (β, ℓ) = calc_β_ℓ(w, β_lu)
    num_missing = count(isnan, β)
    if num_missing > 0
        if num_missing == 1
            idx = findfirst(isnan, β)
            β[idx] = 2 - sum(filter(!isnan, β))
        else
            throw(ArgumentError("missing $(num_missing-1) left-turn angles"))
        end
    end
    s = classify_symmetry(w, β, ℓ)
    Polygon(w, s, β, ℓ)
end

function make_rotation!(w_base::SVector{B}, β_lu, ::CyclicSymmetry{R}) where {B,R}
    N = B * R
    # apply rotational symmetry to vertices
    w = MVector{N,ComplexF64}(undef)
    for i ∈ 0:(R-1)
        r = cispi(2i // R)
        w[(1+B*i):(B*(i+1))] .= r .* w_base
    end
    # extend left-turn angle lookup to full polygon
    for (k, β) ∈ pairs(β_lu)
        for i ∈ 1:(R-1)
            β_lu[B*i+k] = β
        end
    end
    SVector(w)
end
make_rotation!(w::SVector, β_lu, ::NoSymmetry) = w

"Construct Polygon with cyclic symmetry"
function Polygon(
    w_base::SVector{B},
    s::CyclicSymmetry{R},
    β_lu::Dict{Int,<:Number} = Dict{Int,Float64}(),
) where {B,R}
    w = make_rotation!(w_base, β_lu, s)
    # calculate left-turn angles and edge lengths for full polygon
    (β, ℓ) = calc_β_ℓ(w, β_lu)
    # try to infer the value of missing left-turn angles (indicated as NaNs)
    num_missing = count(isnan, β)
    if num_missing > 0
        if num_missing == R
            missing_value = (2 - sum(filter(!isnan, β))) / R
            β[findall(isnan, β)] .= missing_value
        else
            nd = num_missing ÷ R - 1
            throw(ArgumentError("missing $nd left-turn angles"))
        end
    end
    Polygon(w, s, β, ℓ)
end

"""
Assumes that the mirror image continues in the same order as the base,
i.e., the last vertex in the base is connected to the first vertex in the
image.
"""
function make_mirror!(w::SVector{B}, β_lu, s::DihedralSymmetry{<:Any,P}) where {B,P}
    normalised_axis = s.axis / abs(s.axis)
    mirror = map(point -> if isinf(point)
        point
    else
        normalised_axis * conj(normalised_axis' * point)
    end, w)
    # Is [begin] or [end] on the symmetry axis?
    # There cannot be 2 Infs next to each other left and right of the axis,
    # therefore if there is an Inf at the boundary, it has to be on the axis.
    rng = if P == 0
        B:-1:1
    elseif P == 1
        if is_on(s.axis, w[begin]) || isinf(w[begin])
            B:-1:2
        else
            (B-1):-1:1
        end
    elseif P == 2
        (B-1):-1:2
    end
    for (k, β) ∈ pairs(β_lu)
        β_lu[1+B+rng.start-k] = β
    end
    SVector{2B-P}(w..., mirror[rng]...)
end

function Polygon(
    w_base::SVector{B,W},
    symmetry::DihedralSymmetry{R,P},
    β_lu::Dict{Int,<:Number} = Dict{Int,Float64}(),
) where {B,W,R,P}
    w_rotbase = make_mirror!(w_base, β_lu, symmetry)
    w = make_rotation!(w_rotbase, β_lu, CyclicSymmetry{R}())

    (β, ℓ) = calc_β_ℓ(w, β_lu)
    num_missing = count(isnan, β)
    if num_missing > 0
        missing_value = (2 - sum(filter(!isnan, β)))
        if num_missing == R
            β[findall(isnan, β)] .= missing_value / R
        elseif num_missing == 2R && count(isinf, w_base) == 1
            β[findall(isnan, β)] .= missing_value / 2R
        else
            throw(ArgumentError("missing left-turn angles"))
        end
    end
    Polygon(w, symmetry, β, ℓ)
end

num_independent_vertices(N, ::CyclicSymmetry{R}) where {R} = N ÷ R
num_independent_vertices(N, ::DihedralSymmetry{R,P}) where {R,P} = P + (N ÷ R - P) ÷ 2
num_independent_vertices(poly::Polygon{N}) where {N} = num_independent_vertices(N, poly.s)

"""First independent vertex

For P=0 we search for the finite vertices whose edge is divided by the axis.
pattern: lr

For P=1 and P=2 we search for patterns
     l.r
     l∞r
    l∞.∞r
"""
first_independent_vertex(::Polygon) = 1
function first_independent_vertex(poly::Polygon{N,<:DihedralSymmetry{R}}) where {N,R}
    axes = [cispi(k // R) * poly.s.axis for k ∈ 0:(R-1)]
    findfirst(
        i -> begin
            w = poly.w[i]
            if isinf(w)
                w₋ = poly.w[mod1(i - 1, N)]
                w₊ = poly.w[mod1(i + 1, N)]
                any(ax -> which_side(ax, w₋) * which_side(ax, w₊) < 0, axes)
            else
                any(ax -> is_on(ax, w), axes)
            end
        end,
        1:N,
    )
end
function first_independent_vertex(poly::Polygon{N,<:DihedralSymmetry{R,0}}) where {N,R}
    axes = [cispi(k // R) * poly.s.axis for k ∈ 0:(R-1)]
    findfirst(i -> begin
        w₋ = poly.w[mod1(i - 1, N)]
        w = poly.w[i]
        if isinf(w₋) || isinf(w)
            false
        else
            any(ax -> which_side(ax, w₋) * which_side(ax, w) < 0, axes)
        end
    end, 1:N)
end
