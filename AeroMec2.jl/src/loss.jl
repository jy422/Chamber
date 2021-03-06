# This file is licensed under the MIT "Expat" License:

# Copyright (c) 2020: Matthew Ozon.

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



# wall deposition: OK
function WallDeposition!(dx_wall::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble,1},xi::Array{Cdouble,1})
    ws.LR[:] = xi[:]
    dx_wall[:] = -ws.gamma0*ws.t0*xi.*x # I cannot explain this behavior #TODO: report this to the community
end


# jacobian of the wall deposition term: OK
function jacobian_wall_deposition!(F_ev::Array{Cdouble,2},ws::AeroSys,x::Array{Cdouble,1},xi::Array{Cdouble,1})
    ws.LR[:] = xi[:]
    # may replace the following loop using diagind  idx = diagind(F_ev); F_ev[idx] = F_ev[idx] - (ws.gamma0*ws.t0)*ws.LR # the list of index should be known beforehand
    for i in 1:ws.nbin
        # derivation w.r.t. the concentration distribution
        F_ev[i,i]   = F_ev[i,i]   - (ws.gamma0*ws.t0)*xi[i]
    end
end


function WallDepositionRate(ws::AeroSys)
    rr = ws.cst_v^-0.5
    sx = 15.0e-9
    dia0 = 110.0e-9
    (3.0e-16/sqrt(pi*(ws.d0)^3/6.0))*rr.^(collect(1:ws.nbin)) .+ 5.0e-5./(1.0.+exp.(-(ws.d.-dia0)/sx))
end
