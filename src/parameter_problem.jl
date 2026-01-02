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
free_params(::Polygon{N,CyclicSymmetry{R}}) where {N,R} = @MVector zeros(N ÷ R)
free_params(::Polygon{N,<:DihedralSymmetry{R,P}}) where {N,R,P} =
    @MVector zeros((N ÷ R - P) ÷ 2)

prevertex_params(v::MVector{V}, ::CyclicSymmetry{R}) where {V,R} =
    MVector{R * V}(ntuple(i -> v[mod1(i, V)], R * V))
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,0}) where {R,V} =
    prevertex_params(MVector{2V}(v..., -v[end:-1:1]...), CyclicSymmetry{R}())
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,1}) where {R,V} =
    prevertex_params(MVector{2V+1}(v..., 0, -v[end:-1:1]...), CyclicSymmetry{R}())
prevertex_params(v::MVector{V}, ::DihedralSymmetry{R,2}) where {R,V} =
    prevertex_params(MVector{2V+2}(v..., 0, -v[end:-1:1]..., 0), CyclicSymmetry{R}())

circshift_noalloc_poplast(v::MVector{V}, Δ::Int) where {V} =
    MVector{V - 1}(ntuple(i -> v[mod1(i - Δ, V)], V - 1))

prevertices(v, ::NoSymmetry, ::Int) = v |> y_to_θ
prevertices(v, s::DihedralSymmetry{<:Any,2}, kN) =
    circshift_noalloc_poplast(prevertex_params(v, s), kN) |> y_to_θ
prevertices(v, s, kN) = circshift_noalloc_poplast(prevertex_params(v, s), kN+1) |> y_to_θ

struct ProblemIndices{X,L}
    kN::Int
    k_fix::SVector{X,Int}
    k_len::SVector{L,Int}

    function ProblemIndices(poly::Polygon{N}) where {N}
        idx₁ = first_independent_vertex(poly)
        n = num_independent_vertices(poly)
        range = SVector{n}(mod1(idx₁ + i, N) for i ∈ 0:(n-1))

        # last one reaches outside the symmetry-generating set
        k_len = filter(k -> isfinite(poly.ℓ[k]), range)
        # use the first one as kN
        kN = popfirst!(k_len)
        # the connected vertex to fix the rotation and translation
        k1 = mod1(kN + 1, N)

        # we need one fixed vertex per segment between infinities
        k_fix = Int64[k1]
        # add the one before each infinity / last in range
        k_d = range[1]
        skip_current = (k_d == k1)
        for k ∈ range[2:end]
            if k == k1
                skip_current = true
            end
            if isinf(poly.w[k]) || k == range[end]
                if !skip_current
                    push!(k_fix, k_d)
                end
                skip_current = false
            end
            k_d = k
        end

        # need indices of finite ℓ
        num_x = length(k_fix)
        k_fix = SVector{num_x}(k_fix)
        num_ℓ = length(free_params(poly)) - 2 * num_x
        k_len = SVector{num_ℓ}(k_len[1:num_ℓ])
        new{num_x,num_ℓ,Int64}(kN, k_fix, k_len)
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

function solve_parameter_problem(x₀::StaticVector{0}, poly)
    # trivial case does not require NLsolve
    k = findfirst(isfinite, poly.w)
    θ₀ = prevertices(x₀, poly.s, k)
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
