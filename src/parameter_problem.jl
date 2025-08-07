export sc_parameter_problem

using SchwarzChristoffelDisk: Polygon
using NLsolve, StaticArrays

"""Convert `y` to `θ` such that `θ` is sorted

:param y: (n-1)-element array of unconstrained parameters
"""
function y_to_θ(y)
    y2 = 1 .+ cumsum(exp.(-cumsum(y)))
    t = 2π / y2[end]
    [t; t .* y2]
end

"""Fix the scaling constant of the Schwarz-Christoffel transformation and set the singularities from `y`
"""
function sc_fix!(f, y, wN)
    y2 = cumsum([1; exp.(-cumsum(y))])
    t = 2 / y2[end]
    @. f.z = cispi(t * y2)
    @assert !any(isnan(x) for x ∈ f.z) "sc_fix! got NaN output $(f.z) from $y"

    # evaluate image vertices with unit constant
    f.c = 1
    # fix scaling constant such that wN is correct
    f.c = wN / sc_trafo(f, f.z[end])
end

function sc_parameter_problem(polygon::Polygon{N,W,B,L}) where {N,W,B,L}
    y₀ = @MVector zeros(N - 1)
    f = SchwarzChristoffel(y_to_θ(y₀), polygon.β)

    k_inf = findall(isinf, polygon.w)
    # vertex positions
    k_fix = tuple([1; k_inf[2:end] .- 1]...)
    # n-2m-1 side lengths
    k_fin = tuple(findall(!isinf, polygon.ℓ)[1:(N-2*length(k_fix)-1)]...)

    function cost_function!(F, y)
        sc_fix!(f, y, polygon.w[end])

        if !isnothing(F)
            # fix one vertex per component - 2 real conditions each
            @inbounds for (i, k) ∈ enumerate(k_fix)
                dk = polygon.w[k] - sc_trafo(f, f.z[k])
                F[2i-1] = real(dk)
                F[2i] = imag(dk)
            end
            # fix n-2m-1 side lengths - 1 real condition each
            m = length(k_fix)
            @inbounds for (j, k) ∈ enumerate(k_fin)
                F[2m+j] = abs(sc_segment(f, k)) - polygon.ℓ[k]
            end
        end
    end

    # solve
    sol = nlsolve(cost_function!, y₀)
    # apply solution
    sc_fix!(f, sol.zero, polygon.w[end])

    # test
    !sc_test_ok(f, polygon.w) && @warn "parameter_problem failed"

    (sol, f)
end
