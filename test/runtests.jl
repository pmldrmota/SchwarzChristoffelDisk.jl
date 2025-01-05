module SCTest

using Test
using SchwarzChristoffelDisk
using TaylorSeries

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

    taylor = sc_taylor_series(sc, z2)
    @test getcoeff(taylor, 1) ≈ complex(1.0691869766674431, -0.12708193500252732)
    @test getcoeff(taylor, 2) ≈ complex(0.024420875807490826, -0.30506287944280847)
    @test getcoeff(taylor, 3) ≈ complex(0.08634209350133563, 0.05909308433623032)
end

end  # module
