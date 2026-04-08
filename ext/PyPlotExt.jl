module PyPlotExt

import SchwarzChristoffelDisk: supports_backend, _draw_lines!, _sc_plot
using SchwarzChristoffelDisk: PlotSpec, LineStyle
using PyPlot, Colors

supports_backend(::Val{:pyplot}) = true

const _to_pyplot_linestyle = Dict(
    nothing => "solid",
    :solid => "solid",
    :dash => "dashed",
    :dot => "dotted",
    :dashdot => "dashdotted",
    :dashdotdot => "dashdotdotted",
)

function _backend_kwargs(s::LineStyle)
    rgb = nothing
    if !isnothing(s.color)
        c = RGB(s.color)
        rgb = (c.r, c.g, c.b)
    end
    kwargs = Dict(
        :color => rgb,
        :linewidth => s.linewidth,
        :linestyle => _to_pyplot_linestyle[s.linestyle],
    )
    if !isnothing(s.backend_options)
        merge!(kwargs, pairs(s.backend_options))
    end
    kwargs
end

function _draw_lines!(::Val{:pyplot}, spec::PlotSpec, ax)
    ax.set_xlim((-spec.extent, spec.extent))
    ax.set_ylim((-spec.extent, spec.extent))
    ax.set_aspect("equal")
    for line ∈ spec.lines
        style = _backend_kwargs(line.style)
        ax.plot(line.xs, line.ys; style...)
    end
end

function _sc_plot(::Val{:pyplot}, spec::PlotSpec, prevertices)
    fig, ax = subplots(1, 2)
    ax[1].set_title("Polygon")
    ax[2].set_aspect("equal")

    _draw_lines!(Val(:pyplot), spec, ax[1])

    # prevertices
    ax[2].scatter(real.(prevertices), imag.(prevertices); color = :black)
    # unit circle
    zcirc = cispi.(range(0, 2, 100))
    ax[2].plot(real.(zcirc), imag.(zcirc); color = :black)
    ax[2].set_title("Prevertices")
    return (fig, ax)
end

end