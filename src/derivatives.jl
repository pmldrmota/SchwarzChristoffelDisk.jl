export sc_first_derivative, sc_taylor_series

using TaylorSeries

sc_term(zk, α, zv) = (1 - zv / zk)^α

"""Derivative of the Schwarz-Christoffel transformation

:param z: source vertices
:param β: left-turn angles in units of π
:param zv: source point to evaluate at
:param k: index to skip (used for Gauss-Jacobi quadrature)
"""
function sc_derivative(z, β, zv, k=nothing)
    val = one(ComplexF64)
    for (i, (zk, βk)) in enumerate(zip(z, β))
        i == k && continue
        val *= sc_term(zk, -βk, zv)
    end
    val
end

"n-th derivative of `sc_term`"
function sc_term_derivative(zk, α, zv, n)
    p = prod(α + 1 - i for i in 1:n)
    (-1 / zk)^n * sc_term(zk, α - n, zv) * p
end

"""Evaluate the derivative

:param zv: the function argument
:param k: index of the singularity to ignore
"""
function sc_first_derivative(f, zv, k=nothing)
    f.c * sc_derivative(f.z, f.β, zv, k)
end

"""Derivative of the logarithm of the Schwarz-Christoffel derivative

The logarithmic derivative turns the recursive product rule into a single sum.
"""
function log_second_derivative(f, zv)
    mapreduce(((zk, βk),) -> sc_term_derivative(zk, -βk, zv, 1) * sc_term(zk, +βk, zv),
        +,
        zip(f.z, f.β);
        init=zero(ComplexF64))
end

function quotient_rule(f, zv)
    mapreduce(((zk, βk),) -> begin
            g = sc_term(zk, +βk, zv)  # inverse denominator
            g * sc_term_derivative(zk, -βk, zv, 2) -
            (g * sc_term_derivative(zk, -βk, zv, 1))^2
        end,
        +, zip(f.z, f.β); init=zero(ComplexF64))
end

"""Evaluate the second derivative"""
function sc_second_derivative(f, zv, d1=sc_first_derivative(f, zv))
    d1 * log_second_derivative(f, zv)
end

"""Evaluate the third derivative"""
function sc_third_derivative(f, zv, d1=sc_first_derivative(f, zv), d2=sc_second_derivative(f, zv, d1))
    d1 * quotient_rule(f, zv) + d2^2 / d1
end

"""nth-order Taylor series around `zv`

Use `inverse()` to get the Taylor series of the inverse transformation.
"""
function sc_taylor_series(f, zv, n = 3)
    if n == 0
        return Taylor1([0])
    end

    ∂₁ = sc_first_derivative(f, zv)
    if n == 1
        return Taylor1([0, ∂₁])
    end

    ∂₂ = sc_second_derivative(f, zv, ∂₁)
    if n == 2
        return Taylor1([0, ∂₁, ∂₂ / 2])
    end

    ∂₃ = sc_third_derivative(f, zv, ∂₁, ∂₂)
    if n == 3
        return Taylor1([0, ∂₁, ∂₂ / 2, ∂₃ / 6])
    end

    throw("$n-th order Taylor series not implemented")
end
