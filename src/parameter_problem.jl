export sc_parameter_problem

using NLsolve, StaticArrays

"""Convert `y` to `θ` such that `θ` is sorted

:param y: (n-1)-element array of unconstrained parameters
"""
function y_to_θ(y)
    y2 = cumsum([1; exp.(-cumsum(y))])
    t = 2 / y2[end]
    t .* y2
end

"""Fix the scaling constant of the Schwarz-Christoffel transformation and set the singularities
"""
function sc_fix!(f, θ, kN, wN)
    f.z .= cispi.(θ)
    any(isnan, f.z) && error("sc_fix! got NaN output from θ=$θ.")

    # evaluate image vertices with unit constant
    f.c = 1
    # fix scaling constant such that wN is correct
    f.c = wN / sc_trafo(f, f.z[kN])
end

# Dispatch free parameters based on the typed symmetry.
free_params(::Polygon{N,NoSymmetry}) where {N} = @MVector zeros(N - 1)
function free_params(poly::Polygon{N,CyclicSymmetry{R}}) where {N,R}
    V = N ÷ R
    if V == 1
        # the location of this vertex is taken care of by the integration constant
        @MVector zeros(0)
    elseif V == 2
        z = (isinf(poly.w[1]) || isinf(poly.w[2])) ? 0 : 1
        @MVector zeros(z)
    else
        @MVector zeros(V-1)
    end
end
free_params(::Polygon{N,<:DihedralSymmetry{R,P}}) where {N,R,P} =
    @MVector zeros((N ÷ R - P) ÷ 2)

prevertex_repeat(v::MVector{V}, ::Val{R}) where {V,R} =
    MVector{R * V}(ntuple(i -> v[mod1(i, V)], R * V))

# Dispatch the prevertices based on the typed symmetry.
prevertex_params(v::MVector{V}, ::CyclicSymmetry{R}) where {V,R} =
    prevertex_repeat(MVector{V+1}(v..., -sum(v)), Val(R))
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,0}) where {R,V} =
    prevertex_repeat(MVector{2V}(v..., -v[end:-1:1]...), Val(R))
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,1}) where {R,V} =
    prevertex_repeat(MVector{2V+1}(0, v..., -v[end:-1:1]...), Val(R))
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,2}) where {R,V} =
    prevertex_repeat(MVector{2V+2}(0, v..., 0, -v[end:-1:1]...), Val(R))

circshift_noalloc_poplast(v::MVector{V}, Δ::Int) where {V} =
    MVector{V-1}(ntuple(i -> v[mod1(i - Δ, V)], V - 1))

prevertices(v, ::NoSymmetry, ::Int) = v |> y_to_θ
prevertices(v, s::AbstractSymmetry, Δ::Int) =
    circshift_noalloc_poplast(prevertex_params(v, s), Δ) |> y_to_θ

""" ProblemIndices

  - Δ: The prevertex shift Δ is (idx₁-1), where idx₁ is the first independent
    vertex of the polygon.
  - kN: The index of the vertex which will be fixed by the integration constant
    of the Schwarz-Christoffel transformation.
  - k_fix: List of indices of vertices which will be fixed in both coordinates.
  - k_len: List of indices of edges (starting at the vertex with that index) which
    will be fixed to the desired length.
"""
struct ProblemIndices{X,L}
    Δ::Int  # prevertex shift
    kN::Int  # vertex fixed by integration constant
    k_fix::SVector{X,Int}  # vertices fixed by cost function
    k_len::SVector{L,Int}  # edges fixed by cost function
end

"""Find the first occurrence for which predicate evaluates to true

The start and steps limit the elements considered to [start:(start+steps)].
"""
function findnext_circ(
    predicate::Function,
    A::StaticVector{N},
    start::Integer,
    steps::Integer = N,
) where {N}
    ff = findfirst(i -> predicate(A[mod1(start - 1 + i, N)]), 1:steps)
    isnothing(ff) && return nothing
    mod1(start - 1 + ff, N)
end

"""Find all occurrences for which predicate evaluates to true

The start and steps limit the elements considered to [start:(start+steps)].
The search results are naturally in circular order.
"""
function findall_circ(
    predicate::Function,
    A::StaticVector{N},
    start::Integer,
    steps::Integer = N,
) where {N}
    ff = findall(i -> predicate(A[mod1(start - 1 + i, N)]), 1:steps)
    mod1.(start - 1 .+ ff, N)
end

function ProblemIndices(poly::Polygon{N,CyclicSymmetry{R}}) where {N,R}
    num_free = length(free_params(poly))

    if num_free == 1
        # shortcut for trivial case with 1 free parameter -> fix any finite edge
        kN = findfirst(isfinite, poly.ℓ)
        return ProblemIndices{0,1}(0, kN, SVector{0}(), SA[kN])
    end

    num_independent = num_independent_vertices(poly)
    k∞_base = findall_circ(isinf, poly.w, 1, num_independent)
    num_∞_base = length(k∞_base)
    num_fix = max(1, num_∞_base - (R > 1))
    num_len = num_free - 2 * num_fix

    if num_∞_base == 0
        # standard case: fix one edge with kN and k_fix, and rest with k_len
        return ProblemIndices{1,num_len}(0, N, SA[1], 1:num_len)
    end

    kN = findfirst(isfinite, poly.ℓ)
    k₁ = mod1(kN+1, N)
    k∞ = findall_circ(isinf, poly.w, k₁)
    k_fix = [k₁; mod1.(k∞[1:(num_fix-1)] .+ 1, N)]
    k_len = findall_circ(isfinite, poly.ℓ, k₁)[1:num_len]

    ProblemIndices{num_fix,num_len}(0, kN, k_fix, k_len)
end

function ProblemIndices(poly::Polygon{N,<:DihedralSymmetry{R,P}}) where {N,R,P}
    # For P=0, there is definitely no infinity on any symmetry axis.
    # Therefore, the rotational degree of freedom is constrained by kN and the symmetric
    # vertex, which means we don't need k_fix in the same segment.
    # For P=1, idx₁ is always on axis and on the other axis, there is definitely no infinity.
    # Therefore, the rotational DOF is again constrained from that segment.
    idx₁ = first_independent_vertex(poly)  # always on axis
    # the total number of parameters and constraints needs to be `2 * length(k_fix) + length(k_len)`
    num_free = length(free_params(poly))

    if num_free == 1
        # shortcut for trivial case with 1 free parameter -> fix one finite edge
        start = mod1(idx₁ + (P > 0), N)
        kN = findnext_circ(isfinite, poly.ℓ, start)
        return ProblemIndices{0,1}(idx₁-1, kN, SVector{0}(), SA[kN])
    end

    # however, because each fixed vertex comes with 2 constraints, we need to round down, in case
    # we hit the upper limit from the number of segments below
    max_num_fix = fld(num_free, 2)
    ∞_on_axes = num_infs_on_axes(poly)
    num_∞ = count(isinf, poly.w)
    num_segments = (num_∞ - R * ∞_on_axes) ÷ 2R + 1
    min_num_fix = min(max_num_fix, num_segments - 1)
    num_len = num_free - 2 * min_num_fix
    num_independent = num_independent_vertices(poly)

    # start out with the minimal number of k_fix (one per segment)
    k_fix = (findall_circ(isinf, poly.w, idx₁+1) .+ 1)[1:min_num_fix]
    # start out with the maximum number of k_len in the symmetry base
    k_len = findall_circ(isfinite, poly.ℓ, idx₁, num_independent)
    # if there are not enough k_len in the symmetry base, add k_fix
    if length(k_len) < num_len
        push!(k_fix, mod1(popfirst!(k_len) + 1, N))
        num_len -= 2
    end
    keepat!(k_len, 1:num_len)

    if num_∞ == 0
        kN = k_len[1]
    else
        # choose kN such that it fixes the first segment but doesn't coincide with
        # the symmetry axis unless absolutely necessary
        ∞₁ = findnext_circ(isinf, poly.w, idx₁+1)
        kN = mod1(∞₁ - 1, N)
        while kN ∈ k_fix || mod1(kN - 1, N) ∈ k_len || isinf(poly.w[kN])
            kN = mod1(kN - 1, N)
        end
    end
    ProblemIndices{length(k_fix),length(k_len)}(idx₁-1, kN, k_fix, k_len)
end

function cost_function!(F, x, f, poly, idxs::ProblemIndices{0,1})
    sc_fix!(f, prevertices(x, poly.s, idxs.Δ), idxs.kN, poly.w[idxs.kN])
    k = idxs.k_len[1]
    F[1] = abs(sc_segment(f, k)) - poly.ℓ[k]
end

function cost_function!(F, x, f, poly, idxs::ProblemIndices{X,L}) where {X,L}
    sc_fix!(f, prevertices(x, poly.s, idxs.Δ), idxs.kN, poly.w[idxs.kN])

    # fix one vertex per component - 2 real conditions each
    @inbounds for (i, k) ∈ enumerate(idxs.k_fix)
        dk = poly.w[k] - sc_trafo(f, f.z[k])
        F[2i-1] = real(dk)
        F[2i] = imag(dk)
    end
    # fix n-2m-1 side lengths - 1 real condition each
    @inbounds for (j, k) ∈ enumerate(idxs.k_len)
        F[2X+j] = abs(sc_segment(f, k)) - poly.ℓ[k]
    end
end

function solve_parameter_problem(::StaticVector{0}, poly::Polygon{N}) where {N}
    # trivial case does not require NLsolve
    k = findfirst(isfinite, poly.w)
    θ₀ = y_to_θ(zeros(N - 1))
    f = SchwarzChristoffel(θ₀, poly.β)
    sc_fix!(f, θ₀, k, poly.w[k])
    sc_test_ok(f, poly.w) || @warn "parameter_problem failed"
    (nothing, f)
end

function solve_parameter_problem(x₀::StaticVector{F}, poly) where {F}
    idxs = ProblemIndices(poly)
    @show idxs
    f = SchwarzChristoffel(prevertices(x₀, poly.s, idxs.Δ), poly.β)
    # solve
    sol = nlsolve((F, x) -> cost_function!(F, x, f, poly, idxs), x₀; ftol = 1e-10)
    # apply solution
    sc_fix!(f, prevertices(sol.zero, poly.s, idxs.Δ), idxs.kN, poly.w[idxs.kN])
    sc_test_ok(f, poly.w) || @warn "parameter_problem failed"
    (sol, f)
end

# Dispatch on number of free parameters
sc_parameter_problem(poly::Polygon{N}) where {N} =
    solve_parameter_problem(free_params(poly), poly)
