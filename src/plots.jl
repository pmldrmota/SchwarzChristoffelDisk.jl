export draw_polygon, sc_plot

using StaticArrays, PyPlot

function draw_polygon(p::Polygon, ax)
    N = Size(p.w)[1]
    γ = angle(p.w[begin] - p.w[end]) .+ π .* cumsum(p.β)
    v = [p.w[end]]
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    R = max(xlim[2] - xlim[1], ylim[2] - ylim[1])
    for (k, wk) ∈ enumerate(p.w)
        if isinf(wk)
            ax.plot(real.(v), imag.(v); color=:black)
            v = [v[end], v[end] + R * cis(γ[mod1(k - 1, N)])]
            ax.plot(real.(v), imag.(v); color=:black, linestyle=":")
            v = [p.w[mod1(k + 1, N)] - R * cis(γ[k]), p.w[mod1(k + 1, N)]]
            ax.plot(real.(v), imag.(v); color=:black, linestyle=":")
            empty!(v)
        else
            push!(v, wk)
        end
    end
    ax.plot(real.(v), imag.(v); color=:black)
end

function sc_draw_polygon!(f, ax)
    N = Size(f.z)[1]
    ws = SVector{N,ComplexF64}([sc_trafo(f, zk) for zk ∈ f.z])
    poly = Polygon(ws, f.β)
    draw_polygon(poly, ax)

    # fix plot range to a square
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    maxrange = max(xlim[2] - xlim[1], ylim[2] - ylim[1])

    function rescale(lims)
        μ = (lims[2] + lims[1]) / 2
        (μ - maxrange / 2, μ + maxrange / 2)
    end

    ax.set_xlim(rescale(xlim))
    ax.set_ylim(rescale(ylim))
end

function sc_plot(f, rpoints, θpoints)
    trafo(z) = sc_trafo(f, z)

    fig, ax = subplots(1, 2)
    ax[1].set_aspect("equal")
    ax[2].set_aspect("equal")

    # plot the transformation
    for r ∈ range(0, 1, rpoints + 2)[begin+1:end-1]
        θ = range(0, 2π, ceil(Int64, sqrt(r) * θpoints))
        zs = cis.(θ) .* r
        ws = trafo.(zs)
        ax[1].plot(real.(ws), imag.(ws))
    end

    # polygon
    sc_draw_polygon!(f, ax[1])

    # center point
    ax[1].scatter([0], [0])

    # prevertices
    ax[2].scatter(real.(f.z), imag.(f.z))
    # unit circle
    zcirc = cispi.(range(0, 2, θpoints))
    ax[2].plot(real.(zcirc), imag.(zcirc))

    display(fig)
    close(fig)
end
