module SCTest

using Test
using SchwarzChristoffelDisk
using TaylorSeries
using StaticArrays

@testset "Derivatives" begin
    z1 = complex(0.2, 0.4)
    z2 = complex(0.2, 0.1)
    θ = [0, π / 2, π]
    β = fill(2 // 3, 3)
    sc = SchwarzChristoffel(θ, β)

    @test SchwarzChristoffelDisk.sc_term(1im, -1 // 2, z1) ≈
          complex(1.2411967672541269, -0.20141850719855625)

    @test SchwarzChristoffelDisk.sc_term_derivative(1im, -1 // 2, z1, 1) ≈
          complex(-0.4613630722124488, -0.880542948640956)

    @test SchwarzChristoffelDisk.sc_term_derivative(1im, -1 // 2, z1, 2) ≈
          complex(-1.635199330282814, 1.6984741239587269)

    @test SchwarzChristoffelDisk.sc_term_derivative(1im, -1 // 2, z1, 3) ≈
          complex(8.413277127698743, 4.008904833612143)

    @test SchwarzChristoffelDisk.sc_derivative(cis.(θ), β, z2) ≈
          complex(1.0691869766674431, -0.1270819350025273)

    @test SchwarzChristoffelDisk.sc_first_derivative(sc, z2) ≈
          complex(1.0691869766674431, -0.1270819350025273)

    @test SchwarzChristoffelDisk.sc_second_derivative(sc, z2) ≈
          complex(0.04884175161498165, -0.6101257588856169)

    @test SchwarzChristoffelDisk.sc_third_derivative(sc, z2) ≈
          complex(0.5180525610080138, 0.3545585060173819)

    for n = 0:3
        taylor = @test_nowarn sc_taylor_series(sc, z2, n)
        @test getcoeff(taylor, 0) ≈ 0
        if n ≥ 1
            @test getcoeff(taylor, 1) ≈ complex(1.0691869766674431, -0.12708193500252732)
        end
        if n ≥ 2
            @test getcoeff(taylor, 2) ≈ complex(0.024420875807490826, -0.30506287944280847)
        end
        if n ≥ 3
            @test getcoeff(taylor, 3) ≈ complex(0.08634209350133563, 0.05909308433623032)
        end
    end
end

@testset "Polygon" begin
    # Test whether error is thrown if N < 3
    @test_throws ArgumentError Polygon(SA[-1, 1])

    # NoSymmetry (finite)
    poly = Polygon(SA[0, 1, 1+1im, 2im])
    @test all(poly.β .== [0.5, 0.5, 0.25, 0.75])
    @test all(poly.ℓ .== [1, 1, sqrt(2), 2])
    @test poly.s isa NoSymmetry
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) == 1

    # NoSymmetry (infinite)
    @test_throws ArgumentError Polygon(SA[0, 1, 1+1im, Inf, 2im])
    poly = Polygon(SA[0, 1, 1+1im, Inf, 2im], Dict(3 => -0.25, 4 => 1.25, 5 => 0))
    @test all(poly.β .== [0.5, 0.5, -0.25, 1.25, 0])
    @test all(poly.ℓ .== [1, 1, Inf, Inf, 2])
    @test poly.s isa NoSymmetry
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) == 1

    # CyclicSymmetry (finite)
    poly = Polygon(SA[1im, -1+1im], CyclicSymmetry{2}())
    ϕ = acos(1/sqrt(5))/π
    κ = sqrt(5)
    @test all(poly.β .== [ϕ, 1-ϕ, ϕ, 1-ϕ])
    @test all(poly.ℓ .== [1, κ, 1, κ])
    @test poly.s isa CyclicSymmetry{2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) == 1

    # CyclicSymmetry (infinite)
    poly = Polygon(SA[1im, Inf, -1+1im], CyclicSymmetry{2}(), Dict(1 => -ϕ, 3 => ϕ))
    @test all(poly.β .== [-ϕ, 1, ϕ, -ϕ, 1, ϕ])
    @test all(poly.ℓ .== [Inf, Inf, sqrt(5), Inf, Inf, sqrt(5)])
    @test poly.s isa CyclicSymmetry{2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) == 1

    # BilateralSymmetry{0} (finite)
    poly = Polygon(SA[-1+1im, -2], BilateralSymmetry{0}(1.0im))
    @test all(poly.β .== [0.25, 0.75, 0.75, 0.25])
    @test all(poly.ℓ .== [sqrt(2), 4, sqrt(2), 2])
    @test poly.s isa BilateralSymmetry{0}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 3]

    # BilateralSymmetry{0} (infinite)
    poly =
        Polygon(SA[-1+1im, Inf, -2], BilateralSymmetry{0}(1im), Dict(1 => -0.25, 3 => 0.25))
    @test all(poly.β .== [-0.25, 1, 0.25, 0.25, 1, -0.25])
    @test all(poly.ℓ .== [Inf, Inf, 4, Inf, Inf, 2])
    @test poly.s isa BilateralSymmetry{0}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 4]

    # BilateralSymmetry{1} (finite)
    poly = Polygon(SA[1im, -1], BilateralSymmetry{1}(1im))
    @test all(poly.β .== [0.5, 0.75, 0.75])
    @test all(poly.ℓ .== [sqrt(2), 2, sqrt(2)])
    @test poly.s isa BilateralSymmetry{1}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 3]

    # BilateralSymmetry{1} (infinity not on axis)
    poly = Polygon(SA[1im, Inf, -1], BilateralSymmetry{1}(1im), Dict(1 => -0.5, 3 => -0.25))
    @test all(poly.β .== [-0.5, 1.5, -0.25, -0.25, 1.5])
    @test all(poly.ℓ .== [Inf, Inf, 2, Inf, Inf])
    @test poly.s isa BilateralSymmetry{1}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) == 1

    # BilateralSymmetry{1} (infinity on axis)
    poly = Polygon(SA[Inf, -1], BilateralSymmetry{1}(1im), Dict(2 => 0.25))
    @test all(poly.β .== [1.5, 0.25, 0.25])
    @test all(poly.ℓ .== [Inf, 2, Inf])
    @test poly.s isa BilateralSymmetry{1}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) == 1

    # BilateralSymmetry{2} (finite)
    poly = Polygon(SA[1im, -1, -2im], BilateralSymmetry{2}(1im))
    @test all(poly.β .== [0.5, 0.75-ϕ, 2ϕ, 0.75-ϕ])
    @test all(poly.ℓ .== [sqrt(2), sqrt(5), sqrt(5), sqrt(2)])
    @test poly.s isa BilateralSymmetry{2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 3]

    # BilateralSymmetry{2} (one infinity on axis)
    poly = Polygon(SA[1im, -1, Inf], BilateralSymmetry{2}(1im), Dict(2 => 0.25))
    @test all(poly.β .== [0.5, 0.25, 1, 0.25])
    @test all(poly.ℓ .== [sqrt(2), Inf, Inf, sqrt(2)])
    @test poly.s isa BilateralSymmetry{2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 3]

    # BilateralSymmetry{2} (both infinities on axis)
    poly = Polygon(
        SA[Inf, -1+1im, -1-1im, Inf],
        BilateralSymmetry{2}(1im),
        Dict(1 => 2, 2 => -0.5, 3 => -0.25),
    )
    @test all(poly.β .== [2, -0.5, -0.25, 1.5, -0.25, -0.5])
    @test all(poly.ℓ .== [Inf, 2, Inf, Inf, 2, Inf])
    @test poly.s isa BilateralSymmetry{2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 4]

    # DihedralSymmetry{2,0} (finite)
    poly = Polygon(SA[0.5+1im], DihedralSymmetry{2,0}(1im))
    @test all(poly.β .== [0.5, 0.5, 0.5, 0.5])
    @test all(poly.ℓ .== [1, 2, 1, 2])
    @test poly.s isa DihedralSymmetry{2,0}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)

    # DihedralSymmetry{2,0} (infinite)
    poly = Polygon(
        SA[0.5+1im, Inf, 1+2im],
        DihedralSymmetry{2,0}(-1im),
        Dict(1 => -0.5, 2 => 1.25, 3 => -0.25),
    )
    @test all(
        poly.β .==
        [-0.5, 1.25, -0.25, -0.25, 1.25, -0.5, -0.5, 1.25, -0.25, -0.25, 1.25, -0.5],
    )
    @test all(poly.ℓ .== [Inf, Inf, 2, Inf, Inf, 2, Inf, Inf, 2, Inf, Inf, 2])
    @test poly.s isa DihedralSymmetry{2,0}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 4, 7, 10]

    # DihedralSymmetry{2,1} (finite)
    poly = Polygon(SA[1+1im, 2im], DihedralSymmetry{2,1}(1im))
    @test all(poly.β .== [0.25, 0.5, 0.25, 0.25, 0.5, 0.25])
    @test all(poly.ℓ .== [sqrt(2), sqrt(2), 2, sqrt(2), sqrt(2), 2])
    @test poly.s isa DihedralSymmetry{2,1}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 2, 4, 5]

    # DihedralSymmetry{2,1} (infinite)
    poly = Polygon(SA[1+1im, Inf], DihedralSymmetry{2,1}(1im), Dict(1 => -0.25))
    @test all(poly.β .== [-0.25, 1.5, -0.25, -0.25, 1.5, -0.25])
    @test all(poly.ℓ .== [Inf, Inf, 2, Inf, Inf, 2])
    @test poly.s isa DihedralSymmetry{2,1}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 2, 4, 5]

    # DihedralSymmetry{4,2} (finite)
    poly = Polygon(SA[1+1im, 0.5+1im, 2im], DihedralSymmetry{4,2}(2im))
    @test all(
        poly.β .== [0.5, -ϕ, 2ϕ, -ϕ, 0.5, -ϕ, 2ϕ, -ϕ, 0.5, -ϕ, 2ϕ, -ϕ, 0.5, -ϕ, 2ϕ, -ϕ],
    )
    @test all(poly.ℓ .== [1, κ, κ, 1, 1, κ, κ, 1, 1, κ, κ, 1, 1, κ, κ, 1] / 2)
    @test poly.s isa DihedralSymmetry{4,2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test isodd(first_independent_vertex(poly))

    # DihedralSymmetry{2,2} (finite)
    poly = Polygon(SA[2, 1+0.5im, 2.5im], DihedralSymmetry{2,2}(2im))
    @test all(poly.β .≈ [2ϕ, 0.5-2ϕ, 2ϕ, 0.5-2ϕ, 2ϕ, 0.5-2ϕ, 2ϕ, 0.5-2ϕ])
    @test all(poly.ℓ .== [κ/2, κ, κ, κ/2, κ/2, κ, κ, κ/2])
    @test poly.s isa DihedralSymmetry{2,2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test isodd(first_independent_vertex(poly))

    # DihedralSymmetry{2,2} (infinite, infinities not on axes)
    poly = Polygon(
        SA[2, 1+1im, Inf, 2.5im],
        DihedralSymmetry{2,2}(1im),
        Dict(2 => -0.5, 4 => -0.9),
    )
    @test all(poly.β .== [0.5, -0.5, 1.2, -0.9, 1.2, -0.5, 0.5, -0.5, 1.2, -0.9, 1.2, -0.5])
    @test all(
        poly.ℓ .==
        [sqrt(2), Inf, Inf, Inf, Inf, sqrt(2), sqrt(2), Inf, Inf, Inf, Inf, sqrt(2)],
    )
    @test poly.s isa DihedralSymmetry{2,2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 4, 7, 10]

    # DihedralSymmetry{2,2} (infinite, infinity on one axis)
    poly = Polygon(SA[2, 1+0.5im, Inf], DihedralSymmetry{2,2}(2im), Dict(2 => -0.75))
    @test all(poly.β .≈ [2ϕ, -0.75, 2.5-2ϕ, -0.75, 2ϕ, -0.75, 2.5-2ϕ, -0.75])
    @test all(poly.ℓ .== [κ/2, Inf, Inf, κ/2, κ/2, Inf, Inf, κ/2])
    @test poly.s isa DihedralSymmetry{2,2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test isodd(first_independent_vertex(poly))

    # DihedralSymmetry{2,2} (infinite, infinities on both axes)
    poly = Polygon(
        SA[Inf, 2+0.5im, 1+1.5im, Inf],
        DihedralSymmetry{2,2}(2im),
        Dict(2 => -0.5, 3 => -0.25, 4 => 1.2),
    )
    @test all(
        poly.β .== [1.3, -0.5, -0.25, 1.2, -0.25, -0.5, 1.3, -0.5, -0.25, 1.2, -0.25, -0.5],
    )
    @test all(
        poly.ℓ .==
        [Inf, sqrt(2), Inf, Inf, sqrt(2), Inf, Inf, sqrt(2), Inf, Inf, sqrt(2), Inf],
    )
    @test poly.s isa DihedralSymmetry{2,2}
    @test poly.s == classify_symmetry(poly.w, poly.β, poly.ℓ)
    @test first_independent_vertex(poly) ∈ [1, 4, 7, 10]
end

@testset "Transformation" begin
    # Square
    f = SchwarzChristoffel(SA[0, π/2, π, 3π/2], SA[1/2, 1/2, 1/2, 1/2])
    a = 1.3110287771438591
    @test sc_test_ok(f, [a, 1im * a, -a, -1im * a])
end

@testset "ParameterProblem" begin
    # Finite case
    poly = Polygon(SA[1.0+0.5im, 0.1+1.0im, -(1.0+0.5im), -(0.1+1.0im)])
    @test_nowarn sc_parameter_problem(poly)

    # Infinite case
    poly = Polygon(SA[1.0, 1.0im, Inf], Dict([1, 2] .=> [1/3, 1/3]))
    @test_nowarn sc_parameter_problem(poly)
end

@testset "InverseTransformation" begin
    f = SchwarzChristoffel(SA[0, π/2, π], SA[1, 1/2, 1/2])
    z = -0.876 + 0.229
    w = sc_trafo(f, z)
    @test sc_inv(f, w) ≈ z
end

end  # module
