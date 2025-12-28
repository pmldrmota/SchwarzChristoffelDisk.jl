export draw_polygon, sc_plot

using StaticArrays, PyPlot

function draw_polygon(p::Polygon, ax; kwargs...)
    N = Size(p.w)[1]
    k = findfirst(idx₁ -> begin
        idx₂ = mod1(idx₁ + 1, N)
        isfinite(p.w[idx₁]) && isfinite(p.w[idx₂])
    end, 1:N)
    w = circshift(p.w, -k)
    β = circshift(p.β, -k)
    γ = angle(w[begin] - w[end]) .+ π .* cumsum(β)
    v = [w[end]]
    R = 10 * maximum(filter(!isinf, abs.(w)))
    for (k, wk) ∈ enumerate(w)
        if isinf(wk)
            ax.plot(real.(v), imag.(v); kwargs...)
            v = [v[end], v[end] + R * cis(γ[mod1(k - 1, N)])]
            ax.plot(real.(v), imag.(v); kwargs..., linestyle = ":")
            v = [w[mod1(k + 1, N)] - R * cis(γ[k]), w[mod1(k + 1, N)]]
            ax.plot(real.(v), imag.(v); kwargs..., linestyle = ":")
            empty!(v)
        else
            push!(v, wk)
        end
    end
    ax.plot(real.(v), imag.(v); kwargs...)
end

function fill_polygon(p::Polygon, ax; color = color, kwargs...)
    N = Size(p.w)[1]
    k = findfirst(idx₁ -> begin
        idx₂ = mod1(idx₁ + 1, N)
        isfinite(p.w[idx₁]) && isfinite(p.w[idx₂])
    end, 1:N)
    w = circshift(p.w, -k)
    β = circshift(p.β, -k)
    γ = angle(w[begin] - w[end]) .+ π .* cumsum(β)
    R = 10 * maximum(filter(!isinf, abs.(w)))

    function add_draw(x, y)
        ax.fill([x; x[1]], [y; y[1]], color; kwargs...)
    end

    kinf1 = findfirst(isinf, w)
    if isnothing(kinf1)
        add_draw(real.(w), imag.(w))
    else
        v = ComplexF64[]
        for k = kinf1:(kinf1+N+1)
            this = mod1(k, N)
            next = mod1(k + 1, N)
            prev = mod1(k - 1, N)
            if isinf(w[this])
                if k > kinf1
                    push!(v, v[end] + R * cis(γ[prev]))
                    add_draw(real.(v), imag.(v))
                    empty!(v)
                end
                push!(v, w[next] - R * cis(γ[this]))
            else
                push!(v, w[this])
            end
        end
    end
end

function sc_draw_polygon!(f, ax, poly_fill=false; color = :black, kwargs...)
    N = Size(f.z)[1]
    ws = SVector{N,ComplexF64}(sc_trafo(f, zk) for zk ∈ f.z)
    poly = Polygon(ws, Dict(1:N .=> f.β))
    if poly_fill
        fill_polygon(poly, ax; color = color, kwargs...)
    else
        draw_polygon(poly, ax; color = color, kwargs...)
    end

    max_x = maximum(filter(!isinf, abs.(real.(poly.w))))
    max_y = maximum(filter(!isinf, abs.(imag.(poly.w))))
    maxrange = 1.2 * max(max_x, max_y)
    ax.set_xlim((-maxrange, maxrange))
    ax.set_ylim((-maxrange, maxrange))
end

function sc_plot(f, rpoints, θpoints, cmap = "Spectral", poly_fill = false; kwargs...)
    trafo(z) = sc_trafo(f, z)

    fig, ax = subplots(1, 2)
    ax[1].set_aspect("equal")
    ax[2].set_aspect("equal")

    # plot the transformation
    colormap = get_cmap(cmap)
    num_colors = rpoints + 2
    colors = [colormap((i - 1) / (num_colors - 1)) for i = 1:num_colors]

    # polygon
    sc_draw_polygon!(f, ax[1], poly_fill; kwargs...)

    for (i, r) ∈ enumerate(range(0, 1, rpoints + 2)[begin+1:end-1])
        θ = range(0, 2π, ceil(Int64, sqrt(r) * θpoints))
        zs = cis.(θ) .* r
        ws = trafo.(zs)
        ax[1].plot(real.(ws), imag.(ws), color = colors[i+1])
    end

    # center point
    ax[1].scatter([0], [0], marker = ".", color = colors[begin])

    # prevertices
    ax[2].scatter(real.(f.z), imag.(f.z))
    # unit circle
    zcirc = cispi.(range(0, 2, θpoints))
    ax[2].plot(real.(zcirc), imag.(zcirc))

    display(fig)
    close(fig)
end
