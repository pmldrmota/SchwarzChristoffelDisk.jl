export Polygon

using StaticArrays

struct Polygon{N,W,B,L}
    w::SVector{N,W}  # vertices
    β::SVector{N,B}  # left-turn angles
    ℓ::SVector{N,L}  # length of edges [i,i+1]
end

function Polygon(w::SVector{N,W}) where {N,W}
    @assert count(isinf, w) == 0 "must specify angles if there are infinities"
    # preallocate output
    β = zeros(N)
    ℓ = zeros(N)
    for (i, wi) ∈ enumerate(w)
        post = w[mod1(i + 1, N)] - wi
        pre = wi - w[mod1(i - 1, N)]
        β[i] = angle(pre' * post) / π
        ℓ[i] = abs(post)
    end
    Polygon(w, SVector{N}(β), SVector{N}(ℓ))
end

function Polygon(w::SVector{N,W}, β::SVector{N,B}) where {N,W,B}
    ∑β = sum(β)
    @assert ∑β ≈ 2 "wrong angles (∑β=$∑β)"
    k = 0
    # preallocate output
    ℓ = zeros(N)
    for (i, wi) ∈ enumerate(w)
        post = mod1(i + 1, N)
        @assert !(isinf(wi) && isinf(w[post])) "remove consecutive infinity at k=$post"
        if isinf(w[mod1(i - 2, N)]) && !isinf(w[mod1(i - 1, N)]) && !isinf(wi)
            # find circshift such that w[N-1] is an infinity and w[1] and w[N] are finite
            k = i - 1
        end
        ℓ[i] = abs(w[post] - wi)
    end
    if k != 0
        w = circshift(w, -k)
        ℓ = circshift(ℓ, -k)
        β = circshift(β, -k)
    end
    Polygon(SVector{N}(w), SVector{N}(β), SVector{N}(ℓ))
end
