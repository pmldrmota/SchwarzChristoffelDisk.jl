export sc_inv

using NLsolve
using StaticArrays

"""
Calculate the inverse of the Schwarz-Christoffel transformation `f`
"""
function sc_inv(f::SchwarzChristoffel, w::Number; ftol = 1e-10)

    function cost_function!(F, x)
        z = complex(x...)
        w_eval = sc_trafo(f, z)
        F .= reim(w_eval - w)
    end

    sol = nlsolve(cost_function!, @MVector zeros(2); ftol = ftol)
    complex(sol.zero...)
end
