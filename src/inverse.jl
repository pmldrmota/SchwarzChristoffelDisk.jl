export sc_inv

using NLsolve
using StaticArrays


function sc_inv(f, w)

    function cost_function!(F, x)
        z = complex(x...)
        w_eval = sc_trafo(f, z)
        F .= reim(w_eval - w)
    end

    sol = nlsolve(cost_function!, @MVector zeros(2))
    z = complex(sol.zero...)
    @assert sc_trafo(f, z) ≈ w "sc_inv failed"
    z
end