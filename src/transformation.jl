export SchwarzChristoffel, sc_trafo, sc_test_ok

using StaticArrays, FastGaussQuadrature

#: Number of points for Gauss quadrature
const NQUAD = 7

struct StalledException <: Exception
    message::String
end

mutable struct SchwarzChristoffel{Z,B,C,QC}
    z::Z
    β::B
    c::C
    quadcache::QC
end

"""Compute the Gauss quadrature rules needed to compute the Schwarz-Christoffel map
"""
function assemble_quadcache(β)
    jacobi_weights(β) = β < 1 ? SVector{NQUAD}.(gaussjacobi(NQUAD, -β, 0)) : nothing
    qjac = tuple(map(jacobi_weights, β)...)
    qleg = SVector{NQUAD}.(gausslegendre(NQUAD))

    (jacobi_quadrules = qjac, legendre_quadrule = qleg)
end

"""Constructor of the Schwarz-Christoffel derivative

:param θ: the angles of the singularities, in ascending order
:param β: the turning angles of the polygon
:param c: scaling and rotation constant
"""
function SchwarzChristoffel(θ::AbstractVector, β::AbstractVector, c = one(ComplexF64))
    SchwarzChristoffel(cis.(θ), β, c, assemble_quadcache(β))
end

"""Find the distance to the closest singularity

:param zv: the starting point
:param k: index of the singularity to ignore
"""
distance_to_singularity(f, zv, ::Nothing) = sqrt(minimum(abs2.(zv .- f.z)))
function distance_to_singularity(f, zv, k)
    dists = abs2.(zv .- f.z)
    @views sqrt(minimum([dists[begin:(k-1)]; dists[(k+1):end]]))
end

"""Integrate using Gauss quadrature

Takes into account the change of integration limits.

:param integrand: integrand
:param zm: midpoint (b+a)/2
:param zd: radius (b-a)/2
:param (xq,wq): quadrature points and weights
the points are in (-1, 1)
"""
gauss_quadrature(integrand, zm, zd, (xq, wq)) =
    zd * mapreduce(((xqk, wqk),) -> integrand(zm + zd * xqk) * wqk, +, zip(xq, wq))

function sc_integrate(f, za, zb, ::Nothing)
    gauss_quadrature(
        zv -> sc_first_derivative(f, zv),
        (zb + za) / 2,
        (zb - za) / 2,
        f.quadcache.legendre_quadrule,
    )
end

function sc_integrate(f, za, zb, kb)
    jacobi_quadrule = f.quadcache.jacobi_quadrules[kb]
    if isnothing(jacobi_quadrule)
        return complex(Inf)
    else
        zd = (zb - za) / 2
        zm = (zb + za) / 2
        gq = gauss_quadrature(zv -> sc_first_derivative(f, zv, kb), zm, zd, jacobi_quadrule)
        abs(zd)^(-f.β[kb]) * gq
    end
end

"""Compound Gauss-(Jacobi) integral

`zb` may be a singular point (index `kb`). If it is not, then pass `kb=nothing`.
`z0` is not allowed to be a singular point and must be within the section of
the disk enclosed by the connections from the nearest singularities to the origin.
"""
function sc_compound_gauss_jacobi(f, z0, zb, kb)
    @assert !isnan(z0) && !isnan(zb) "sc_compound_gauss_jacobi got NaN input"
    # integrate from z0 to zb
    integral = zero(zb)
    if zb ≈ z0
        return integral
    end

    # remaining distance
    drem = abs(z0 - zb)

    for _ ∈ 1:100
        # maximum integration interval
        dmax = 2 * distance_to_singularity(f, zb, kb)
        if iszero(dmax)
            # This error is thrown if due to finite numerical precision, the
            # distance to the closest singularity is smaller than the floating
            # point resolution of the prevertex coordinates.
            # The solution would be to parameterise everything in terms of the
            # difference to a fixed "anchor", taken to be the closest singularity.
            throw(StalledException("Compound Gauss-Jacobi integration stalled"))
        end

        # integration start point
        za = zb + min(drem, dmax) * (z0 - zb) / drem
        # integrate
        integral += sc_integrate(f, za, zb, kb)

        if dmax ≥ drem
            return integral
        end

        # use legendre quadrule from now on (if not already)
        kb = nothing
        zb = za
        drem = abs(z0 - zb)
    end

    throw(ErrorException("sc_compound_gauss_jacobi failed at ($z0,$zb,$kb); drem = $drem"))
end

"""Integrate along segment connecting adjacent singularities `(ka, ka+1)`
"""
function sc_segment(f, ka)
    kb = mod1(ka + 1, length(f.z))
    # Taking the midpoint doesn't seem to be accurate close to the boundary.
    # zhalf = (f.z[kb]+f.z[ka])/2
    zhalf = zero(f.z[kb])
    sc_compound_gauss_jacobi(f, zhalf, f.z[kb], kb) -
    sc_compound_gauss_jacobi(f, zhalf, f.z[ka], ka)
end

"""Schwarz-Christoffel transformation
"""
function sc_trafo(f::SchwarzChristoffel, zb::Number)
    k = findfirst(zb == zz for zz ∈ f.z)
    sc_compound_gauss_jacobi(f, zb, k)
end

"""Test whether the Schwarz-Christoffel transformation maps to the vertices `w`

:returns: true if test passed
"""
function sc_test_ok(f::SchwarzChristoffel, w::AbstractVector)
    sum(abs2(sc_trafo(f, f.z[k]) - wk) for (k, wk) ∈ enumerate(w) if !isinf(wk)) <
    length(w)^2 * 1e-6
end
