export LineStyle, draw!, sc_contour!, sc_plot

using ColorSchemes, Colors

struct LineStyle
    color::Union{Nothing,Any}
    linewidth::Union{Nothing,Real}
    linestyle::Union{Nothing,Symbol}
    backend_options::Union{Nothing,NamedTuple}

    function LineStyle(color, linewidth, linestyle, backend_options)
        supported_linestyles = (:solid, :dash, :dot, :dashdot, :dashdotdot)
        isnothing(linestyle) ||
            linestyle ∈ supported_linestyles ||
            throw(ArgumentError("Linestyle must be one of $supported_linestyles"))
        if !isnothing(color)
            color = parse(RGB, color)
        end
        new(color, linewidth, linestyle, backend_options)
    end
end
LineStyle(;
    color = nothing,
    linewidth = nothing,
    linestyle = nothing,
    backend_options = nothing,
) = LineStyle(color, linewidth, linestyle, backend_options)

struct LineSpec
    xs::Vector{Float64}
    ys::Vector{Float64}
    style::LineStyle
end
LineSpec(zs::Vector{<:Complex}, style::LineStyle) = LineSpec(real.(zs), imag.(zs), style)

mutable struct PlotSpec
    lines::Vector{LineSpec}
    extent::Float64  # xlim = ylim = (-extent, extent)
end

get_extent(w) = maximum(map(x -> isinf(x) ? 0 : max(abs(real(x)), abs(imag(x))), w))

"""
Valid keyword arguments are
  - enlarge_extent: Factor to enlarge scene by, defaults to 1.2
  - ray_style: LineStyle specifying the style of infinite rays
  - edge_style: LineStyle specifying the style of finite edges
"""
function polygon_primitives(
    poly::Polygon{N};
    edge_style::LineStyle = LineStyle(linestyle = :solid),
    ray_style::LineStyle = LineStyle(linestyle = :dot),
    enlarge_extent::Real = 1.2,
) where {N}
    extent = enlarge_extent * get_extent(poly.w)
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

get_levels(levels::Integer) = range(0, 1, levels + 1)[(begin+1):end]
get_levels(levels::AbstractVector) = collect(levels)

"""Recursive refinement of contour line

Samples the polar-angle mid-point until the distance between adjacent mapped
points is below `min_distance`.
"""
function refine_segment(
    f::SchwarzChristoffel,
    start::Pair{<:Real,<:Complex},
    finish::Pair{<:Real,<:Complex},
    min_distance::Real,
    radius::Real,
)
    if abs(finish.second - start.second) ≤ min_distance
        return [start, finish]
    end
    θ_mid = (start.first + finish.first) / 2
    mid = Pair(θ_mid, sc_trafo(f, cis(θ_mid) * radius))
    return [
        refine_segment(f, start, mid, min_distance, radius);
        refine_segment(f, mid, finish, min_distance, radius)[2:end]
    ]
end

function sample_contour(f::SchwarzChristoffel, min_distance::Real, radius::Real)
    num_init = 8
    points = map(θ -> Pair(θ, sc_trafo(f, cis(θ) * radius)), range(0, 2π, num_init))
    vcat(
        points[1],
        [
            refine_segment(f, points[i], points[i+1], min_distance, radius)[2:end] for
            i ∈ 1:(num_init-1)
        ]...,
    )
end

function contour_primitives(
    f::SchwarzChristoffel;
    levels::Union{Integer,AbstractVector} = 15,
    color = nothing,
    colorscheme::Symbol = :Spectral,
    enlarge_extent::Real = 1.2,
    backend_options::Union{Nothing,NamedTuple} = nothing,
)
    poly = Polygon(f)
    rlevels = get_levels(levels)
    num_levels = length(rlevels)
    cscheme = colorschemes[colorscheme]
    colors = if isnothing(color)
        map(x -> get(cscheme, x), rlevels)
    else
        fill(parse(RGB, color), num_levels)
    end
    lines = LineSpec[]
    min_distance = get_extent(poly.w) / 50
    for (i, r) ∈ enumerate(filter(<(1), rlevels))
        contour = sample_contour(f, min_distance, r)
        style = LineStyle(color = colors[i], backend_options = backend_options)
        push!(lines, LineSpec(map(p -> p.second, contour), style))
    end
    # todo: avoid Polygon construction
    # → would benefit from evaluating f along rays (especially in polygons without
    # finite edges)
    style = LineStyle(color = colors[end], backend_options = backend_options)
    spec = polygon_primitives(
        poly;
        edge_style = style,
        ray_style = style,
        enlarge_extent = enlarge_extent,
    )
    if rlevels[end] == 1
        prepend!(spec.lines, lines)
    else
        # hacky way of preserving the extent as if the outline was drawn.
        empty!(spec.lines)
        append!(spec.lines, lines)
    end
    spec
end

"Backend capability probe; extensions override this."
supports_backend(::Val) = false

"Backend fallback; extensions override specific Val methods."
_draw_lines!(::Val, spec::PlotSpec, ax; kwargs...) =
    error("Requested plotting backend is not available.")
_sc_plot(::Val, spec::PlotSpec, prevertices; kwargs...) =
    error("Requested plotting backend is not available.")

function auto_backend(backend::Symbol)
    supported_backends = (:makie, :plots, :pyplot)
    if backend === :auto
        for b in supported_backends
            if supports_backend(Val(b))
                return b
            end
        end
        error("No supported plotting backend is loaded.")
    elseif backend ∈ supported_backends
        return backend
    else
        throw(ArgumentError("Backend $b is not supported."))
    end
end

draw!(poly::Polygon, ax; backend::Symbol = :auto, kwargs...) =
    _draw_lines!(Val(auto_backend(backend)), polygon_primitives(poly; kwargs...), ax)

sc_contour!(f::SchwarzChristoffel, ax; backend::Symbol = :auto, kwargs...) =
    _draw_lines!(Val(auto_backend(backend)), contour_primitives(f; kwargs...), ax)

sc_plot(f::SchwarzChristoffel; backend::Symbol = :auto, kwargs...) =
    _sc_plot(Val(auto_backend(backend)), contour_primitives(f; kwargs...), f.z)