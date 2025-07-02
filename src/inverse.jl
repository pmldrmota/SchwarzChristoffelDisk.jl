export sc_inv

using NLsolve
using StaticArrays


function sc_inv(f, w)

    function cost_function!(F, x)
        z = complex(x...)
        w_eval = sc_trafo(f, z)
        F .= reim(w_eval - w)
    end

    sol = nlsolve(cost_function!, @MVector zeros(2); ftol = 1e-10)
    complex(sol.zero...)
end