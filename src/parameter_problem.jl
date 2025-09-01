export sc_parameter_problem

using SchwarzChristoffelDisk: Polygon
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
function sc_fix!(f, θ, wN)
    f.z .= cispi.(θ)
    @assert !any(isnan.(f.z)) "sc_fix! got NaN output $(f.z) from $θ"

    # evaluate image vertices with unit constant
    f.c = 1
    # fix scaling constant such that wN is correct
    f.c = wN / sc_trafo(f, f.z[end])
end

# dispatch free parameters based on typed symmetry
free_params(::Polygon{N,NoSymmetry}) where {N} = @MVector zeros(N - 1)
free_params(::Polygon{N,CyclicSymmetry{R}}) where {N,R} = @MVector zeros(N ÷ R)
free_params(::Polygon{N,<:BilateralSymmetry{P}}) where {N,P} = @MVector zeros((N - P) ÷ 2)
free_params(::Polygon{N,<:DihedralSymmetry{R,P}}) where {N,R,P} =
    @MVector zeros((N ÷ R - P) ÷ 2)

# dispatch prevertices using typed symmetry
prevertices(v::MVector, ::NoSymmetry) = v |> y_to_θ
prevertices(v::MVector{V,T}, ::CyclicSymmetry{R}) where {V,T,R} =
    MVector{R * V - 1,T}(ntuple(i -> v[mod1(i, V)], R * V - 1)) |> y_to_θ
prevertices(v::MVector{V,T}, ::BilateralSymmetry{0}) where {V,T} =
    MVector{2 * V - 1,T}(v..., -v[end:-1:2]...) |> y_to_θ
prevertices(v::MVector{V,T}, ::BilateralSymmetry{1}) where {V,T} =
    MVector{2 * V,T}(v..., 0, -v[end:-1:2]...) |> y_to_θ
prevertices(v::MVector{V,T}, ::BilateralSymmetry{2}) where {V,T} =
    MVector{2 * V + 1,T}(v..., 0, -v[end:-1:1]...) |> y_to_θ
prevertices(v::MVector{V,T}, ::DihedralSymmetry{R,0}) where {V,T,R} =
    prevertices(MVector{2 * V,T}(v..., -v[end:-1:1]...), CyclicSymmetry{R}())
prevertices(v::MVector{V,T}, ::DihedralSymmetry{R,1}) where {V,T,R} =
    prevertices(MVector{2 * V + 1,T}(v..., 0, -v[end:-1:1]...), CyclicSymmetry{R}())
prevertices(v::MVector{V,T}, ::DihedralSymmetry{R,2}) where {V,T,R} =
    prevertices(MVector{2 * V + 2,T}(v..., 0, -v[end:-1:1]..., 0), CyclicSymmetry{R}())

function sc_parameter_problem(poly::Polygon{N}) where {N}
    x₀ = free_params(poly)
    f = SchwarzChristoffel(prevertices(x₀, poly.s), poly.β)

    k_inf = findall(isinf, poly.w)
    # vertex positions
    k_fix = tuple(1, (k_inf[2:end] .- 1)...)
    # n-2m-1 side lengths
    k_fin = tuple(findall(isfinite, poly.ℓ)[1:(N-2*length(k_fix)-1)]...)

    function cost_function!(F, x)
        sc_fix!(f, prevertices(x, poly.s), poly.w[end])

        if !isnothing(F)
            # fix one vertex per component - 2 real conditions each
            @inbounds for (i, k) ∈ enumerate(k_fix)
                dk = poly.w[k] - sc_trafo(f, f.z[k])
                F[2i-1] = real(dk)
                F[2i] = imag(dk)
            end
            # fix n-2m-1 side lengths - 1 real condition each
            m = length(k_fix)
            @inbounds for (j, k) ∈ enumerate(k_fin)
                F[2m+j] = abs(sc_segment(f, k)) - poly.ℓ[k]
            end
        end
    end

    # solve
    sol = nlsolve(cost_function!, x₀)
    # apply solution
    sc_fix!(f, prevertices(sol.zero, poly.s), poly.w[end])

    # test
    !sc_test_ok(f, poly.w) && @warn "parameter_problem failed"

    (sol, f)
end
