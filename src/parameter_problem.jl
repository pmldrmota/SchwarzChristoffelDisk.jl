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

# dispatch free parameters based on typed symmetry
free_params(::Polygon{N,NoSymmetry}) where {N} = @MVector zeros(N - 1)
cyclic_free_params(::Val{1}) = @MVector zeros(0)
cyclic_free_params(::Val{2}) = @MVector zeros(1)
cyclic_free_params(::Val{V}) where {V} = @MVector zeros(V)
free_params(::Polygon{N,CyclicSymmetry{R}}) where {N,R} = cyclic_free_params(Val(N ÷ R))
free_params(::Polygon{N,<:DihedralSymmetry{R,P}}) where {N,R,P} =
    @MVector zeros((N ÷ R - P) ÷ 2)

# For cyclic symmetry, it is important that the cost function connects vertices
# across the base, where the base is implicitly defined in this prevertex_params
# as vertices 1:R (to naturally work with DihedralSymmetry too).
prevertex_params(v::MVector{V}, ::CyclicSymmetry{R}) where {V,R} =
    MVector{R * V}(ntuple(i -> v[mod1(i, V)], R * V))
# Special case: cyclic symmetry with 2 vertices in base always gives [A,-A,A,-A...]
prevertex_params(v::MVector{1}, ::CyclicSymmetry{R}) where {R} =
    MVector{2R}(ntuple(i -> iseven(i) ? v[1] : -v[1], 2R))
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,0}) where {R,V} =
    prevertex_params(MVector{2V}(v..., -v[end:-1:1]...), CyclicSymmetry{R}())
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,1}) where {R,V} =
    prevertex_params(MVector{2V+1}(0, v..., -v[end:-1:1]...), CyclicSymmetry{R}())
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,2}) where {R,V} =
    prevertex_params(MVector{2V+2}(0, v..., 0, -v[end:-1:1]...), CyclicSymmetry{R}())

circshift_noalloc_poplast(v::MVector{V}, Δ::Int) where {V} =
    MVector{V-1}(ntuple(i -> v[mod1(i - Δ, V)], V - 1))

prevertices(v, ::NoSymmetry, ::Int) = v |> y_to_θ
prevertices(v, s::CyclicSymmetry, ::Int) =
    circshift_noalloc_poplast(prevertex_params(v, s), 0) |> y_to_θ
prevertices(v, s::DihedralSymmetry{<:Any,0}, kN) =
    circshift_noalloc_poplast(prevertex_params(v, s), kN-1) |> y_to_θ
prevertices(v, s::DihedralSymmetry, kN) =
    circshift_noalloc_poplast(prevertex_params(v, s), kN-2) |> y_to_θ

struct ProblemIndices{X,L}
    kN::Int
    k_fix::SVector{X,Int}
    k_len::SVector{L,Int}

    function ProblemIndices(poly::Polygon{N}) where {N}
        idx₁ = first_independent_vertex(poly)
        num_free = length(free_params(poly))
        num_infs = count(isinf, poly.w)

        if num_infs == 0
            kN = idx₁
            if num_free == 1
                return new{0,1}(kN, SVector{0}(), SVector{1}(kN))
            else
                k_fix = SVector{1}(mod1(kN + 1, N))
                num_ℓ = num_free - 2
                k_len = SVector{num_ℓ}(mod1(kN + i, N) for i ∈ 1:num_ℓ)
                return new{1,num_ℓ}(kN, k_fix, k_len)
            end
        else
            throw(ErrorException("ProblemIndices not implemented for infinities"))
        end
    end
end

function cost_function!(F, x, f, poly, idxs::ProblemIndices{0,1})
    sc_fix!(f, prevertices(x, poly.s, idxs.kN), idxs.kN, poly.w[idxs.kN])
    k = idxs.k_len[1]
    F[1] = abs(sc_segment(f, k)) - poly.ℓ[k]
end

function cost_function!(F, x, f, poly, idxs::ProblemIndices{X,L}) where {X,L}
    sc_fix!(f, prevertices(x, poly.s, idxs.kN), idxs.kN, poly.w[idxs.kN])

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
    f = SchwarzChristoffel(prevertices(x₀, poly.s, idxs.kN), poly.β)
    # solve
    sol = nlsolve((F, x) -> cost_function!(F, x, f, poly, idxs), x₀; ftol = 1e-10)
    # apply solution
    sc_fix!(f, prevertices(sol.zero, poly.s, idxs.kN), idxs.kN, poly.w[idxs.kN])
    sc_test_ok(f, poly.w) || @warn "parameter_problem failed"
    (sol, f)
end

# Dispatch on number of free parameters
sc_parameter_problem(poly::Polygon{N}) where {N} =
    solve_parameter_problem(free_params(poly), poly)
