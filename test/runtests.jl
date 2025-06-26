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
    # Finite case
    poly = Polygon(SA[0, 1, 1im])
    @test poly.β == [0.5, 0.75, 0.75]
    @test poly.ℓ == [1.0, sqrt(2), 1.0]

    # Infinite case
    poly = Polygon(SA[1.0, 1.0im, Inf], SA[1/4, 1/4, 2-2/4])
    @test poly.β == [0.25, 1.5, 0.25]
    @test poly.ℓ == [Inf, Inf, sqrt(2)]
end

@testset "Transformation" begin
    # Square
    f = SchwarzChristoffel(SA[0, π/2, π, 3π/2], SA[1/2, 1/2, 1/2, 1/2])
    a = 1.3110287771438591
    @test sc_test_ok(f, [a, 1im * a, -a, -1im * a])
end

@testset "ParameterProblem" begin
    # Finite case
    poly = Polygon(SA[1.0+0.5im, 0.1+1.0im, -(1.0 + 0.5im), -(0.1 + 1.0im)])
    @test_nowarn sc_parameter_problem(poly)

    # Infinite case
    poly = Polygon(SA[1.0, 1.0im, Inf], SA[1/3, 1/3, 2-2/3])
    @test_nowarn sc_parameter_problem(poly)
end

end  # module
