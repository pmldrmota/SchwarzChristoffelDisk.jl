module SCTest

using Test
using SchwarzChristoffelDisk
using TaylorSeries
using StaticArrays

# hypotenuse and angle for right-angled triangle with cathetes 1 and 2
const κ = sqrt(5)
const ϕ = acos(1/κ)/π

# @testset "Derivatives" begin
#     z1 = complex(0.2, 0.4)
#     z2 = complex(0.2, 0.1)
#     θ = [0, π / 2, π]
#     β = fill(2 // 3, 3)
#     sc = SchwarzChristoffel(θ, β)

#     @test SchwarzChristoffelDisk.sc_term(1im, -1 // 2, z1) ≈
#           complex(1.2411967672541269, -0.20141850719855625)

#     @test SchwarzChristoffelDisk.sc_term_derivative(1im, -1 // 2, z1, 1) ≈
#           complex(-0.4613630722124488, -0.880542948640956)

#     @test SchwarzChristoffelDisk.sc_term_derivative(1im, -1 // 2, z1, 2) ≈
#           complex(-1.635199330282814, 1.6984741239587269)

#     @test SchwarzChristoffelDisk.sc_term_derivative(1im, -1 // 2, z1, 3) ≈
#           complex(8.413277127698743, 4.008904833612143)

#     @test SchwarzChristoffelDisk.sc_derivative(cis.(θ), β, z2) ≈
#           complex(1.0691869766674431, -0.1270819350025273)

#     @test SchwarzChristoffelDisk.sc_first_derivative(sc, z2) ≈
#           complex(1.0691869766674431, -0.1270819350025273)

#     @test SchwarzChristoffelDisk.sc_second_derivative(sc, z2) ≈
#           complex(0.04884175161498165, -0.6101257588856169)

#     @test SchwarzChristoffelDisk.sc_third_derivative(sc, z2) ≈
#           complex(0.5180525610080138, 0.3545585060173819)

#     for n = 0:3
#         taylor = @test_nowarn sc_taylor_series(sc, z2, n)
#         @test getcoeff(taylor, 0) ≈ 0
#         if n ≥ 1
#             @test getcoeff(taylor, 1) ≈ complex(1.0691869766674431, -0.12708193500252732)
#         end
#         if n ≥ 2
#             @test getcoeff(taylor, 2) ≈ complex(0.024420875807490826, -0.30506287944280847)
#         end
#         if n ≥ 3
#             @test getcoeff(taylor, 3) ≈ complex(0.08634209350133563, 0.05909308433623032)
#         end
#     end
# end

# @testset "Polygon" begin

#     classify_symmetry(p::Polygon) = SchwarzChristoffelDisk.classify_symmetry(p.w, p.β, p.ℓ)

#     # Test whether error is thrown if N < 3
#     @test_throws ArgumentError Polygon(SA[-1, 1])

#     @testset "NoSymmetry" begin
#         @testset "finite" begin
#             @testset let poly = Polygon(SA[0, 1, 1+1im, 2im])
#                 @test all(poly.β .== [0.5, 0.5, 0.25, 0.75])
#                 @test all(poly.ℓ .== [1, 1, sqrt(2), 2])
#                 @test poly.s isa NoSymmetry
#                 @test poly.s == classify_symmetry(poly)
#                 @test first_independent_vertex(poly) == 1
#             end
#         end
#         @testset "infinite" begin
#             @test_throws ArgumentError Polygon(SA[0, 1, 1+1im, Inf, 2im])
#             @testset let poly =
#                     Polygon(SA[0, 1, 1+1im, Inf, 2im], Dict(3 => -0.25, 4 => 1.25, 5 => 0))
#                 @test all(poly.β .== [0.5, 0.5, -0.25, 1.25, 0])
#                 @test all(poly.ℓ .== [1, 1, Inf, Inf, 2])
#                 @test poly.s isa NoSymmetry
#                 @test poly.s == classify_symmetry(poly)
#                 @test first_independent_vertex(poly) == 1
#             end
#         end
#     end

#     @testset "CyclicSymmetry" begin
#         @testset "finite" begin
#             @testset let poly = Polygon(SA[1im, -1+1im], CyclicSymmetry{2}())
#                 @test all(poly.β .== [ϕ, 1-ϕ, ϕ, 1-ϕ])
#                 @test all(poly.ℓ .== [1, κ, 1, κ])
#                 @test poly.s isa CyclicSymmetry{2}
#                 @test poly.s == classify_symmetry(poly)
#                 @test first_independent_vertex(poly) == 1
#             end
#         end
#         @testset "infinite" begin
#             @testset let poly = Polygon(
#                     SA[1im, Inf, -1+1im],
#                     CyclicSymmetry{2}(),
#                     Dict(1 => -ϕ, 3 => ϕ),
#                 )
#                 @test all(poly.β .== [-ϕ, 1, ϕ, -ϕ, 1, ϕ])
#                 @test all(poly.ℓ .== [Inf, Inf, κ, Inf, Inf, κ])
#                 @test poly.s isa CyclicSymmetry{2}
#                 @test poly.s == classify_symmetry(poly)
#                 @test first_independent_vertex(poly) == 1
#             end
#         end
#     end

#     @testset "BilateralSymmetry" begin
#         @testset "finite" begin
#             @testset "P=0" begin
#                 @testset let poly = Polygon(SA[-1+1im, -2], BilateralSymmetry{0}(1im))
#                     @test all(poly.β .== [0.25, 0.75, 0.75, 0.25])
#                     @test all(poly.ℓ .== [sqrt(2), 4, sqrt(2), 2])
#                     @test poly.s isa BilateralSymmetry{0}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 3]
#                 end
#             end
#             @testset "P=1" begin
#                 @testset let poly = Polygon(SA[1im, -1], BilateralSymmetry{1}(1im))
#                     @test all(poly.β .== [0.5, 0.75, 0.75])
#                     @test all(poly.ℓ .== [sqrt(2), 2, sqrt(2)])
#                     @test poly.s isa BilateralSymmetry{1}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 3]
#                 end
#             end
#             @testset "P=2" begin
#                 @testset let poly = Polygon(SA[1im, -1, -2im], BilateralSymmetry{2}(1im))
#                     @test all(poly.β .== [0.5, 0.75-ϕ, 2ϕ, 0.75-ϕ])
#                     @test all(poly.ℓ .== [sqrt(2), κ, κ, sqrt(2)])
#                     @test poly.s isa BilateralSymmetry{2}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 3]
#                 end
#             end
#         end

#         @testset "infinite" begin
#             @testset "P=0" begin
#                 @testset let poly = Polygon(
#                         SA[-1+1im, Inf, -2],
#                         BilateralSymmetry{0}(1im),
#                         Dict(1 => -0.25, 3 => 0.25),
#                     )
#                     @test all(poly.β .== [-0.25, 1, 0.25, 0.25, 1, -0.25])
#                     @test all(poly.ℓ .== [Inf, Inf, 4, Inf, Inf, 2])
#                     @test poly.s isa BilateralSymmetry{0}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 4]
#                 end
#             end
#             @testset "P=1 (infinity not on axis)" begin
#                 @testset let poly = Polygon(
#                         SA[1im, Inf, -1],
#                         BilateralSymmetry{1}(1im),
#                         Dict(1 => -0.5, 3 => -0.25),
#                     )
#                     @test all(poly.β .== [-0.5, 1.5, -0.25, -0.25, 1.5])
#                     @test all(poly.ℓ .== [Inf, Inf, 2, Inf, Inf])
#                     @test poly.s isa BilateralSymmetry{1}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) == 1
#                 end
#             end
#             @testset "P=1 (infinity on axis)" begin
#                 @testset let poly =
#                         Polygon(SA[Inf, -1], BilateralSymmetry{1}(1im), Dict(2 => 0.25))
#                     @test all(poly.β .== [1.5, 0.25, 0.25])
#                     @test all(poly.ℓ .== [Inf, 2, Inf])
#                     @test poly.s isa BilateralSymmetry{1}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) == 1
#                 end
#             end
#             @testset "P=2 (one infinity on axis)" begin
#                 @testset let poly = Polygon(
#                         SA[1im, -1, Inf],
#                         BilateralSymmetry{2}(1im),
#                         Dict(2 => 0.25),
#                     )
#                     @test all(poly.β .== [0.5, 0.25, 1, 0.25])
#                     @test all(poly.ℓ .== [sqrt(2), Inf, Inf, sqrt(2)])
#                     @test poly.s isa BilateralSymmetry{2}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 3]
#                 end
#             end
#             @testset "P=2 (both infinities on axis)" begin
#                 @testset let poly = Polygon(
#                         SA[Inf, -1+1im, -1-1im, Inf],
#                         BilateralSymmetry{2}(1im),
#                         Dict(1 => 2, 2 => -0.5, 3 => -0.25),
#                     )
#                     @test all(poly.β .== [2, -0.5, -0.25, 1.5, -0.25, -0.5])
#                     @test all(poly.ℓ .== [Inf, 2, Inf, Inf, 2, Inf])
#                     @test poly.s isa BilateralSymmetry{2}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 4]
#                 end
#             end
#         end
#     end

#     @testset "DihedralSymmetry" begin
#         @testset "finite" begin
#             @testset "P=0" begin
#                 @testset let poly = Polygon(SA[0.5+1im], DihedralSymmetry{2,0}(1im))
#                     @test all(poly.β .== [0.5, 0.5, 0.5, 0.5])
#                     @test all(poly.ℓ .== [1, 2, 1, 2])
#                     @test poly.s isa DihedralSymmetry{2,0}
#                     @test poly.s == classify_symmetry(poly)
#                 end
#             end
#             @testset "P=1" begin
#                 @testset let poly = Polygon(SA[1+1im, 2im], DihedralSymmetry{2,1}(1im))
#                     b = [0.25, 0.5, 0.25]
#                     @test all(poly.β .== [b; b])
#                     l = [sqrt(2), sqrt(2), 2]
#                     @test all(poly.ℓ .== [l; l])
#                     @test poly.s isa DihedralSymmetry{2,1}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 2, 4, 5]
#                 end
#             end
#             @testset "P=2" begin
#                 @testset let poly =
#                         Polygon(SA[1+1im, 0.5+1im, 2im], DihedralSymmetry{4,2}(2im))
#                     b = [0.5, -ϕ, 2ϕ, -ϕ]
#                     @test all(poly.β .== [b; b; b; b])
#                     l = [1, κ, κ, 1]
#                     @test all(poly.ℓ .== [l; l; l; l] / 2)
#                     @test poly.s isa DihedralSymmetry{4,2}
#                     @test poly.s == classify_symmetry(poly)
#                     @test isodd(first_independent_vertex(poly))
#                 end
#                 @testset let poly =
#                         Polygon(SA[2, 1+0.5im, 2.5im], DihedralSymmetry{2,2}(2im))
#                     b = [2ϕ, 0.5-2ϕ, 2ϕ, 0.5-2ϕ]
#                     @test all(poly.β .≈ [b; b])
#                     l = [κ/2, κ, κ, κ/2]
#                     @test all(poly.ℓ .== [l; l])
#                     @test poly.s isa DihedralSymmetry{2,2}
#                     @test poly.s == classify_symmetry(poly)
#                     @test isodd(first_independent_vertex(poly))
#                 end
#             end
#         end
#         @testset "infinite" begin
#             @testset "P=0" begin
#                 @testset let poly = Polygon(
#                         SA[0.5+1im, Inf, 1+2im],
#                         DihedralSymmetry{2,0}(-1im),
#                         Dict(1 => -0.5, 2 => 1.25, 3 => -0.25),
#                     )
#                     b = [-0.5, 1.25, -0.25, -0.25, 1.25, -0.5]
#                     @test all(poly.β .== [b; b])
#                     l = [Inf, Inf, 2, Inf, Inf, 2]
#                     @test all(poly.ℓ .== [l; l])
#                     @test poly.s isa DihedralSymmetry{2,0}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 4, 7, 10]
#                 end
#             end
#             @testset "P=1 (infinity on axis)" begin
#                 @testset let poly = Polygon(
#                         SA[1+1im, Inf],
#                         DihedralSymmetry{2,1}(1im),
#                         Dict(1 => -0.25),
#                     )
#                     @test all(poly.β .== [-0.25, 1.5, -0.25, -0.25, 1.5, -0.25])
#                     @test all(poly.ℓ .== [Inf, Inf, 2, Inf, Inf, 2])
#                     @test poly.s isa DihedralSymmetry{2,1}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 2, 4, 5]
#                 end
#             end
#             @testset "P=2 (infinities not on axes)" begin
#                 @testset let poly = Polygon(
#                         SA[2, 1+1im, Inf, 2.5im],
#                         DihedralSymmetry{2,2}(1im),
#                         Dict(2 => -0.5, 4 => -0.9),
#                     )
#                     b = [0.5, -0.5, 1.2, -0.9, 1.2, -0.5]
#                     @test all(poly.β .== [b; b])
#                     l = [sqrt(2), Inf, Inf, Inf, Inf, sqrt(2)]
#                     @test all(poly.ℓ .== [l; l])
#                     @test poly.s isa DihedralSymmetry{2,2}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 4, 7, 10]
#                 end
#             end
#             @testset "P=2 (one infinity on axis)" begin
#                 @testset let poly = Polygon(
#                         SA[2, 1+0.5im, Inf],
#                         DihedralSymmetry{2,2}(2im),
#                         Dict(2 => -0.75),
#                     )
#                     b = [2ϕ, -0.75, 2.5-2ϕ, -0.75]
#                     @test all(poly.β .≈ [b; b])
#                     l = [κ/2, Inf, Inf, κ/2]
#                     @test all(poly.ℓ .== [l; l])
#                     @test poly.s isa DihedralSymmetry{2,2}
#                     @test poly.s == classify_symmetry(poly)
#                     @test isodd(first_independent_vertex(poly))
#                 end
#             end
#             @testset "P=2 (both infinities on axes)" begin
#                 @testset let poly = Polygon(
#                         SA[Inf, 2+0.5im, 1+1.5im, Inf],
#                         DihedralSymmetry{2,2}(2im),
#                         Dict(2 => -0.5, 3 => -0.25, 4 => 1.2),
#                     )
#                     b = [1.3, -0.5, -0.25, 1.2, -0.25, -0.5]
#                     @test all(poly.β .== [b; b])
#                     l = [Inf, sqrt(2), Inf, Inf, sqrt(2), Inf]
#                     @test all(poly.ℓ .== [l; l])
#                     @test poly.s isa DihedralSymmetry{2,2}
#                     @test poly.s == classify_symmetry(poly)
#                     @test first_independent_vertex(poly) ∈ [1, 4, 7, 10]
#                 end
#             end
#         end
#     end
# end

# @testset "Transformation" begin
#     # Square
#     f = SchwarzChristoffel(SA[0, π/2, π, 3π/2], SA[1/2, 1/2, 1/2, 1/2])
#     a = 1.3110287771438591
#     @test sc_test_ok(f, [a, 1im * a, -a, -1im * a])
# end

@testset "ParameterProblem" begin

    function test_circshift_poly(poly::Polygon{N}) where {N}
        for k ∈ 0:(N-1)
            @testset let poly_shifted = Polygon(
                    SVector{N}(ntuple(i -> poly.w[mod1(i - k, N)], N)),
                    poly.s,
                    SVector{N}(ntuple(i -> poly.β[mod1(i - k, N)], N)),
                    SVector{N}(ntuple(i -> poly.ℓ[mod1(i - k, N)], N)),
                )
                (_, f) = sc_parameter_problem(poly_shifted)
                @test sc_test_ok(f, poly_shifted.w)
            end
        end
    end

    @testset "NoSymmetry" begin
        test_circshift_poly.([
            # finite
            Polygon(SA[-1-1im, 1-1im, 1+1im, -1+2im]),
            Polygon(SA[1.0, 1.0im, -1.4-0.4im, -0.23-1im, 0.8-2im]),
            Polygon(SA[0.2-0.5im, 1-0.5im, 2im, -1-0.5im, 0.2-0.5im, 0.2-0.1im]),
            # infinite
            Polygon(
                SA[1.0im, -1+1im, -Inf, -1im, 1-1im, Inf],
                Dict(1 => -0.2, 2 => -0.2, 4 => -0.41, 5 => -0.4, 6 => 1.6),
            ),
            Polygon(
                SA[1.0im, -1+1im, -Inf, -1im, 0.25-1im, 0.5-1.5im, 0.75-1im, 1-1im, Inf],
                Dict(1 => -0.2, 2 => -0.2, 4 => -0.41, 8 => -0.3, 9 => 1.5),
            ),
        ])
    end

    @testset "CyclicSymmetry" begin
        test_circshift_poly.([
            # finite
            Polygon(SA[2+0.1im], CyclicSymmetry{3}()),
            Polygon(SA[1], CyclicSymmetry{4}()),
            Polygon(SA[1], CyclicSymmetry{8}()),
            Polygon(SA[1im, -1+1im], CyclicSymmetry{2}()),
            Polygon(SA[1im, -1+1im], CyclicSymmetry{3}()),
            Polygon(SA[1.0+0.5im, 0.1+1.0im], CyclicSymmetry{2}()),
            Polygon(SA[2, 2+0.3im], CyclicSymmetry{4}()),
            Polygon(SA[2+0.1im, 2+0.2im], CyclicSymmetry{3}()),
            Polygon(SA[1, 0.3*cispi(0.1), 0.8*cispi(0.2)], CyclicSymmetry{3}()),
            Polygon(SA[2, 2+0.2im, 1+1im, 1+2im], CyclicSymmetry{3}()),
            # infinite
            Polygon(SA[1+0.23im, Inf], CyclicSymmetry{4}(), Dict(1 => -0.8)),
            Polygon(
                SA[1-0.5im, 1+2im, Inf],
                CyclicSymmetry{4}(),
                Dict(1 => -0.5, 2 => -0.123),
            ),
            Polygon(SA[1im, Inf, -1+1im], CyclicSymmetry{2}(), Dict(1 => -ϕ, 3 => ϕ)),
            Polygon(
                SA[1-0.5im, 0.3, 1.2+0.5im, Inf],
                CyclicSymmetry{4}(),
                Dict(1 => -0.25, 3 => 0.2),
            ),
            Polygon(
                SA[1-0.5im, 1-0.25im, 0.7+0.3im, 1+1im, Inf],
                CyclicSymmetry{4}(),
                Dict(1 => -0.2, 4 => -0.25),
            ),
            Polygon(
                SA[1-0.5im, 1+0.5im, Inf],
                CyclicSymmetry{4}(),
                Dict(1 => -0.5, 2 => -0.5),
            ),
            Polygon(
                SA[1-0.5im, 0.8, 1+0.5im, Inf],
                CyclicSymmetry{4}(),
                Dict(1 => -0.25, 3 => -0.25),
            ),
            Polygon(
                SA[1-0.5im, 0.75-0.25im, 0.75+0.25im, 1+0.5im, Inf],
                CyclicSymmetry{4}(),
                Dict(1 => -0.25, 4 => -0.25),
            ),
        ])
    end

    @testset "BilateralSymmetry" begin
        # merge with DihedralSymmetry?
        @testset "P=0" begin
            test_circshift_poly.([
                Polygon(
                    SA[-1+0.5im, -1.2-0.2im, -1.5+0.8im, -2-0.5im],
                    BilateralSymmetry{0}(1im),
                ),
            ])
        end
        @testset "P=1" begin
            test_circshift_poly.([
                Polygon(SA[2im, -1+1im, -1.5+1im, -2.1-0.5im], BilateralSymmetry{1}(1im)),
                Polygon(SA[1.0im, -1-2im, 1-2im]),
                Polygon(SA[Inf, -1], BilateralSymmetry{1}(1im), Dict(2 => 0.25)),
                Polygon(SA[Inf, -1-1im], BilateralSymmetry{1}(1im), Dict(2 => 0.25)),
                Polygon(SA[1-1im, 2, 1, 2, 2im], BilateralSymmetry{1}(1im)),
            ])
        end
        @testset "P=2" begin
            test_circshift_poly.([
                Polygon(
                    SA[2im, -1+1im, -1.5+1im, -2.1-0.5im, -1im],
                    BilateralSymmetry{2}(1im),
                ),
                Polygon(SA[2, 1+1im, -2, 1-1im]),
                Polygon(SA[1im, -1, Inf], BilateralSymmetry{2}(1im), Dict(2 => 0.25)),
                Polygon(SA[1im, -2, Inf], BilateralSymmetry{2}(1im), Dict(2 => 0.25)),
                Polygon(
                    SA[1im, -1.5+2im, -2, Inf],
                    BilateralSymmetry{2}(1im),
                    Dict(3 => -0.25),
                ),
            ])
        end
    end

    @testset "DihedralSymmetry" begin
        @testset "P=0" begin
            test_circshift_poly.([
                # finite
                Polygon(SA[0.5+1im], DihedralSymmetry{2,0}(1im)),
                Polygon(SA[1], DihedralSymmetry{2,0}(1+2im)),
                Polygon(
                    SA[1+0.2im, 0.5+1im, 0.2+1.5im, 0.12+1im],
                    DihedralSymmetry{2,0}(1im),
                ),
                Polygon(SA[2-0.4im], DihedralSymmetry{3,0}(2.0+0im)),
                # infinite
                Polygon(
                    SA[2+1im, Inf, 0.5+2.5im],
                    DihedralSymmetry{2,0}(0.1+1im),
                    Dict(1 => -0.5, 3 => -0.5),
                ),
                Polygon(
                    SA[1+0.3im, 2+0.5im, 2+1im, Inf, 2+2im, Inf, 2+3im],
                    DihedralSymmetry{2,0}(1im),
                    Dict(3 => -0.5, 4 => 1, 5 => -0.8, 6=>1.3, 7 => -0.5),
                ),
                Polygon(
                    SA[2+0.2im, 1.5+1im, Inf, 0.5+2.5im],
                    DihedralSymmetry{2,0}(1im),
                    Dict(2 => -0.5, 4 => -0.3),
                ),
                Polygon(
                    SA[2+0.2im, 1.5+1im, Inf, 0.5+2.5im, 0.2+2.1im],
                    DihedralSymmetry{2,0}(1im),
                    Dict(2 => -0.5, 4 => 0.1),
                ),
            ])
        end
        @testset "P=1" begin
            test_circshift_poly.([
                # finite
                Polygon(SA[1+1im, 2im], DihedralSymmetry{2,1}(1im)),
                Polygon(SA[1+1im, 1.2+0.9im, 1.7+1im, 2im], DihedralSymmetry{2,1}(1im)),
                # infinite
                Polygon(SA[2+1im, Inf], DihedralSymmetry{2,1}(0.1+1im), Dict(1 => -0.2)),
                Polygon(
                    SA[2+1im, 1+2im, Inf],
                    DihedralSymmetry{2,1}(1im),
                    Dict(2 => -0.25),
                ),
                Polygon(
                    SA[1+1im, Inf, 2im],
                    DihedralSymmetry{2,1}(1im),
                    Dict(1 => -0.5, 3=>-0.5),
                ),
                Polygon(
                    SA[1+1im, 2+2im, Inf, 2im],
                    DihedralSymmetry{2,1}(1im),
                    Dict(2 => -0.25, 4=>-0.5),
                ),
                Polygon(
                    SA[1+1im, 2+2im, Inf, 1+3im, 2im],
                    DihedralSymmetry{2,1}(1im),
                    Dict(2 => -0.25, 4=>-0.25),
                ),
                Polygon(
                    SA[1+1im, 2+2im, Inf, 1+3im, Inf],
                    DihedralSymmetry{2,1}(1im),
                    Dict(2 => -0.25, 3=>1, 4=>-0.5, 5=>1),
                ),
                Polygon(
                    SA[1+1im, 2+2im, Inf, 2+3im, 1+4im, Inf],
                    DihedralSymmetry{2,1}(1im),
                    Dict(2 => -0.25, 3=>1, 4=>-0.25, 5=>-0.25, 6=>1),
                ),
                Polygon(
                    SA[2+1im, Inf, 1+2im, Inf],
                    DihedralSymmetry{2,1}(1im),
                    Dict(1 => -0.5, 2 => 1, 3 => -0.8),
                ),
                Polygon(
                    SA[2+1im, Inf, 1+2im, 1+2.4im, Inf, 2im],
                    DihedralSymmetry{2,1}(1im),
                    Dict(1 => -0.5, 2 => 1, 3 => -0.5, 4 => -0.2, 5=>1.1, 6 => -0.8),
                ),
                Polygon(
                    SA[2+1im, Inf, 1+2im, 1+2.4im, Inf],
                    DihedralSymmetry{2,1}(1im),
                    Dict(1 => -0.5, 2 => 1, 3 => -0.5, 4 => -0.2),
                ),
                Polygon(
                    SA[2+1im, Inf, 1+2im, 1+2.4im, 1im],
                    DihedralSymmetry{2,1}(1im),
                    Dict(1 => -0.5, 2 => 1),
                ),
                Polygon(
                    SA[2+0.2im, 1.5+1im, Inf, 0.5+2.5im, 0.2+2.1im, Inf],
                    DihedralSymmetry{2,1}(1im),
                    Dict(
                        2 => -0.5,
                        3 => 1.0173595508115136,
                        4 => 0.1,
                        5 => -0.9,
                        6 => 1.209664,
                        7 => -0.9,
                    ),
                ),
            ])
        end
        @testset "P=2" begin
            test_circshift_poly.([
                # finite
                Polygon(SA[1+1im, 0.5+1im, 2im], DihedralSymmetry{4,2}(2im)),
                Polygon(SA[2, 1+0.5im, 2.5im], DihedralSymmetry{2,2}(2im)),
                # infinite
                Polygon(SA[2, 1+1im, Inf], DihedralSymmetry{2,2}(1im), Dict(2 => -0.25)),
                Polygon(
                    SA[2, 1+1im, 0.5+1im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(3 => -0.5),
                ),
                Polygon(SA[Inf, 1+1im, 2im], DihedralSymmetry{2,2}(1im), Dict(2 => -0.25)),
                Polygon(
                    SA[Inf, 2+1im, 1+1im, 2im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(2 => -0.25),
                ),
                Polygon(
                    SA[Inf, 2+1im, 1+2im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2 => -0.25, 3 => -0.25, 4 => 1),
                ),
                Polygon(
                    SA[Inf, 3+1im, 2+2im, 1+2im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2 => -0.25, 4 => -0.5, 5 => 1),
                ),
                Polygon(
                    SA[3, Inf, 2+1im, 1+2im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => -0.5, 2 => 1, 3=>-0.5, 4=>-0.25, 5=>1),
                ),
                Polygon(
                    SA[3, Inf, 2+1im, 1+1im, 1+2im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => -0.5, 2 => 1, 3=>-0.25, 5=>-0.25, 6=>1.5),
                ),
                Polygon(
                    SA[3, 2+1im, Inf, 1+2im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(2 => -0.5, 3=>1, 4=>-0.75, 5=>1),
                ),
                Polygon(
                    SA[2, 1+1im, Inf, 1+2.5im, 0.5+2im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(
                        2 => -0.5,
                        3 => 1.2,
                        4 => -0.2,
                        5 => -0.75,
                        6 => 1,
                        7 => -0.75,
                    ),
                ),
                Polygon(
                    SA[3, 2+1im, Inf, 1+2im, 0.5+2.5im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(2 => -0.5, 3=>1, 4=>-0.5, 5=>-0.25, 6=>1),
                ),
                Polygon(
                    SA[3, 2+1im, 3+2im, Inf, 1+4im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(3 => -0.25, 4=>1, 5=>-0.5, 6=>1),
                ),
                Polygon(
                    SA[3, 2+1im, 3+2im, Inf, 1+4im, 0.5+4.5im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(3 => -0.25, 4=>1, 5=>-0.25, 6=>-0.25, 7=>1),
                ),
                Polygon(
                    SA[Inf, 2+1im, Inf, 1+2im, 2im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.75, 3=>1, 4=>-0.25),
                ),
                Polygon(
                    SA[Inf, 2+1im, Inf, 1+2im, 1+3im, 3im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.75, 3=>1, 4=>-0.75),
                ),
                Polygon(
                    SA[Inf, 2+1im, 1+2im, Inf, 2im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.25, 3=>-0.5, 4=>1, 5=>-0.5),
                ),
                Polygon(
                    SA[Inf, 2+1im, 1+2im, Inf, 1+3im, 3im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.25, 3=>-0.5, 4=>1, 5=>-0.25),
                ),
                Polygon(
                    SA[Inf, 4+1im, Inf, 3+2im, 2+3im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.5, 3 => 1, 4=>-0.75, 5=>-0.25, 6=>1),
                ),
                Polygon(
                    SA[Inf, 4+1im, Inf, 3+3im, 3+2im, 1+2im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.75, 3 => 1, 4=>0.25, 6=>-0.5, 7=>1),
                ),
                Polygon(
                    SA[Inf, 4+1im, 3+2im, Inf, 1+4im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.25, 3=>-0.5, 4=>1, 5=>-0.75, 6=>1),
                ),
                Polygon(
                    SA[Inf, 4+1im, 3+2im, Inf, 2+4im, 1+5im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.25, 3=>-0.5, 4=>1, 5=>-0.5, 6=>-0.25, 7=>1),
                ),
                Polygon(
                    SA[Inf, 4+1im, 3+2im, 3+3im, Inf, 1+5im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.25, 4=>-0.25, 5=>1, 6=>-0.75, 7=>1),
                ),
                Polygon(
                    SA[Inf, 4+1im, 3+2im, 3+3im, Inf, 2+4im, 1+5im, Inf],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1 => 1, 2=>-0.25, 4=>-0.25, 5=>1, 6=>-0.5, 7=>-0.25, 8=>1),
                ),
                Polygon(
                    SA[0.9, Inf, 2+1im, 1.5im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1=>-0.8, 3=>-0.5),
                ),
                Polygon(
                    SA[0.9, Inf, 2+2im, 1+1im, 1.5im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(1=>-0.8, 3=>-0.25),
                ),
                Polygon(
                    SA[2, 1+1im, Inf, 2.5im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(2=>-0.5, 4=>-0.5),
                ),
                Polygon(
                    SA[2, 1+1im, Inf, 2.5im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(2 => -0.5, 4 => -0.9),
                ),
                Polygon(
                    SA[2, 1+1im, Inf, 1+2im, 2.5im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(2=>-0.5, 4=>-0.5),
                ),
                Polygon(
                    SA[2, 1+1im, 2+2im, Inf, 2.5im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(3=>-0.25, 5=>-0.5),
                ),
                Polygon(
                    SA[2, 1+1im, 2+2im, Inf, 1+2im, 2.5im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(3=>-0.25, 5=>-0.5),
                ),
                Polygon(
                    SA[Inf, 2+0.5im, 2+1.5im, Inf],
                    DihedralSymmetry{2,2}(2im),
                    Dict(1 => 1, 2 => -0.5, 3=>-0.4),
                ),
                Polygon(
                    SA[1, 2+1im, Inf, 0.8+1.2im, Inf, 1+2im, 3im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(2 => -0.25, 3 => 1, 4 => -0.75, 5 => 1, 6 => -0.5),
                ),
                Polygon(
                    SA[1, 2+1im, 0.8+1.2im, Inf, 1+2im, 3im],
                    DihedralSymmetry{2,2}(1im),
                    Dict(3 => -0.75, 5 => -0.5),
                ),
            ])
        end
    end
end

# @testset "InverseTransformation" begin
#     f = SchwarzChristoffel(SA[0, π/2, π], SA[1, 1/2, 1/2])
#     z = -0.876 + 0.229
#     w = sc_trafo(f, z)
#     @test sc_inv(f, w) ≈ z
# end

end  # module
