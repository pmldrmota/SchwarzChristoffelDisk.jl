export sc_draw!, sc_plot

struct LineSpec
    xs::Vector{Float64}
    ys::Vector{Float64}
    style::Dict{Symbol,<:Any}
end
LineSpec(zs::Vector{<:Complex}, style::Dict{Symbol,<:Any}) =
    LineSpec(real.(zs), imag.(zs), style)

mutable struct PlotSpec
    lines::Vector{LineSpec}
    extent::Float64  # xlim = ylim = (-extent, extent)
end

get_extent(w; enlarge_extent = 1.2, kwargs...) =
    enlarge_extent * maximum(map(x -> isinf(x) ? 0 : max(abs(real(x)), abs(imag(x))), w))

"""
Valid keyword arguments are
  - enlarge_extent: Factor to enlarge scene by, defaults to 1.2
  - edge_style: Dict{Symbol, Any} specifying the style of finite edges
  - ray_style: Dict{Symbol, Any} specifying the style of infinite rays
"""
function PlotSpec(poly::Polygon{N}; kwargs...) where {N}
    extent = get_extent(poly.w; kwargs...)
    k = findfirst(isfinite, poly.ℓ)
    # If there is no finite edge, don't plot anything.
    # todo: for DihedralSymmetry, infer the orientation from the symmetry axis.
    isnothing(k) && return PlotSpec(LineSpec[], extent)

    # Shift vertices so that we can start with a known orientation.
    w = circshift(poly.w, -k)
    β = circshift(poly.β, -k)
    γ = angle(w[begin] - w[end]) .+ π .* cumsum(β)
    len_inf_ray = 3 * extent
    lines = LineSpec[]
    ray_style = get(kwargs, :ray_style, Dict{Symbol,Any}())
    edge_style = get(kwargs, :edge_style, Dict{Symbol,Any}())
    v = [w[end]]  # Temporary set of vertices contained in the current batch.
    for (k, wk) ∈ enumerate(w)
        if isinf(wk)
            # Commit the current batch.
            push!(lines, LineSpec(v, edge_style))
            # Add the first infinite ray.
            ray_dir = len_inf_ray * cis(γ[mod1(k - 1, N)])
            ray = [v[end], v[end] + ray_dir]
            push!(lines, LineSpec(ray, ray_style))
            # Add the second infinite ray.
            ray_dir = len_inf_ray * cis(γ[k])
            ray = [w[mod1(k + 1, N)] - ray_dir, w[mod1(k + 1, N)]]
            push!(lines, LineSpec(ray, ray_style))
            # Continue with the next finite batch.
            empty!(v)
        else
            push!(v, wk)
        end
    end
    # Commit the last batch.
    push!(lines, LineSpec(v, edge_style))
    PlotSpec(lines, extent)
end

function PlotSpec(
    f::SchwarzChristoffel;
    rpoints = 15,
    θpoints = 500,
    colormap = nothing,
    kwargs...,
)
    N = length(f.z)

    trafo(z) = sc_trafo(f, z)

    num_colors = rpoints + 2
    colors = if isnothing(colormap)
        nothing
    elseif colormap isa Symbol
        fill(colormap, num_colors)
    else
        [colormap((i - 1) / (num_colors - 1)) for i = 1:num_colors]
    end

    lines = LineSpec[]
    for (i, r) ∈ enumerate(range(0, 1, rpoints + 2)[(begin+1):(end-1)])
        θ = range(0, 2π, ceil(Int64, sqrt(r) * θpoints))
        zs = cis.(θ) .* r
        ws = trafo.(zs)
        style = if isnothing(colors)
            Dict{Symbol,Any}()
        else
            Dict(:color => colors[i+1])
        end
        push!(lines, LineSpec(ws, style))
    end

    # todo: avoid Polygon construction
    # → would benefit from evaluating f along rays
    spec = PlotSpec(Polygon(trafo.(f.z), Dict(1:N .=> f.β)))
    prepend!(spec.lines, lines)
    spec
end

"Backend capability probe; extensions override this."
supports_backend(::Val) = false

"Backend fallback; extensions override specific Val methods."
_draw_lines!(::Val, spec::PlotSpec, ax; kwargs...) =
    error("Requested plotting backend is not available.")
_sc_plot(::Val, spec::PlotSpec, prevertices; kwargs...) =
    error("Requested plotting backend is not available.")

function backend_switch(backend_call::Function, backend::Symbol = :auto)
    if backend === :auto
        for b in (:makie, :plots, :pyplot)
            if supports_backend(Val(b))
                return backend_call(Val(b))
            end
        end
        error("No supported plotting backend is loaded.")
    else
        return backend_call(Val(backend))
    end
end

function sc_draw!(obj, ax; backend::Symbol = :auto, kwargs...)
    spec = PlotSpec(obj; kwargs...)
    return backend_switch(x -> _draw_lines!(x, spec, ax; kwargs...), backend)
end

function sc_plot(f::SchwarzChristoffel; backend::Symbol = :auto, kwargs...)
    spec = PlotSpec(f; kwargs...)
    return backend_switch(x -> _sc_plot(x, spec, f.z; kwargs...), backend)
end