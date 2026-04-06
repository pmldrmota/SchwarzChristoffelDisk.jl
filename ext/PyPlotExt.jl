module PyPlotExt

import SchwarzChristoffelDisk: supports_backend, _draw_lines!, _sc_plot, PlotSpec
using PyPlot

supports_backend(::Val{:pyplot}) = true

function _draw_lines!(::Val{:pyplot}, spec::PlotSpec, ax; kwargs...)
    ax.set_xlim((-spec.extent, spec.extent))
    ax.set_ylim((-spec.extent, spec.extent))
    ax.set_aspect("equal")
    for line ∈ spec.lines
        # Override general kwargs with specific line style.
        style = Dict(kwargs..., line.style...)
        # Remove kwargs which are only used internally.
        for key ∈ (:edge_style, :ray_style, :enlarge_extent, :colormap, :rpoints, :θpoints)
            delete!(style, key)
        end
        ax.plot(line.xs, line.ys; style...)
    end
end

function _sc_plot(x::Val{:pyplot}, spec::PlotSpec, prevertices; kwargs...)
    fig, ax = subplots(1, 2)
    ax[1].set_title("Polygon")
    ax[2].set_aspect("equal")

    _draw_lines!(x, spec, ax[1]; kwargs...)

    # prevertices
    ax[2].scatter(real.(prevertices), imag.(prevertices); color = :black)
    # unit circle
    zcirc = cispi.(range(0, 2, 100))
    ax[2].plot(real.(zcirc), imag.(zcirc); color = :black)
    ax[2].set_title("Prevertices")
    display(fig)
    close(fig)
end

end