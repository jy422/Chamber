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



# this file encodes the condensation mechanism
# it is assumed that the size dependence of the condensational growth rate is known and stored in ws.scale_GR
# the growth rate does not depend on the size (ws.GR could be a scalar instead of an array)
# The mechanism are described dimensionless: to retrieve the values with physical unit, one must multiply the results of CondensationGrowth! or jacobian_condensation! by ws.GR0

# condensation
function CondensationGrowth!(dx_cond::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble,1},zeta::Cdouble)
    ws.GR = (ws.GR0*ws.t0/ws.d0)*zeta # dimensionless growth rate
    # for a growth rate indenpend of the size
    dx_cond[1] = -ws.GR*ws.scale_GR[1]*x[1]
    # if ws.logS
        # dx_cond[2:end] = ws.GR*ws.scale_GR[2:end].*(ws.cst_r*x[1:end-1]-x[2:end])  # it a different type of discretization (growth estimated at the boundary of the bin range, it's not bad, but difficult to expalin)
        dx_cond[2:end] = ws.GR*(ws.scale_GR[1:end-1].*ws.cst_r.*x[1:end-1] - ws.scale_GR[2:end].*x[2:end])
    # else # cst_r = 1.0 if logS==false
    #     dx_cond[2:end] = ws.GR*(ws.scale_GR[1:end-1].*x[1:end-1] - ws.scale_GR[2:end].*x[2:end])
    # end
end

function CondensationGrowth!(dx_cond::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble,1},zeta::Array{Cdouble,1})
    GR_array = (ws.GR0*ws.t0/ws.d0)*zeta[:] # dimensionless growth rate
    ws.GR = GR_array[1]
    # for a growth rate indenpend of the size
    dx_cond[1] = -GR_array[1]*ws.scale_GR[1]*x[1]
    # dx_cond[2:end] = GR_array[2:end].*ws.scale_GR[2:end].*(ws.cst_r*x[1:end-1]-x[2:end]) # it a different type of discretization (growth estimated at the boundary of the bin range, it's not bad, but difficult to expalin)
    # if ws.logS
        dx_cond[2:end] = GR_array[1:end-1].*ws.scale_GR[1:end-1].*ws.cst_r.*x[1:end-1] - GR_array[2:end].*ws.scale_GR[2:end].*x[2:end]
    # else # cst_r = 1.0 if logS==false
    #     dx_cond[2:end] = GR_array[1:end-1].*ws.scale_GR[1:end-1].*x[1:end-1] - GR_array[2:end].*ws.scale_GR[2:end].*x[2:end]
    # end
end


# jacobian of the condensation term: OK
function jacobian_condensation!(F_ev_::Array{Cdouble,2},ws::AeroSys,x::Array{Cdouble,1},zeta::Cdouble) # assume that every element of the matrix is set to zero #LATER it might be improved by not computing the addition in the function
    ws.GR = (ws.GR0*ws.t0/ws.d0)*zeta # CGR(zeta) # dimensionless growth rate
    # for a growth rate indenpend of the size (known and stored in ws)
    F_ev_[1,1]   =     F_ev_[1,1]   - ws.GR*ws.scale_GR[1]
    for i in 2:ws.nbin
        # derivation w.r.t. the concentration distribution
        F_ev_[i,i]   = F_ev_[i,i]   - ws.GR*ws.scale_GR[i]
        # F_ev_[i,i-1] = F_ev_[i,i-1] + ws.GR*ws.scale_GR[i]*ws.cst_r
        # if ws.logS
            F_ev_[i,i-1] = F_ev_[i,i-1] + ws.GR*ws.scale_GR[i-1]*ws.cst_r
        # else # cst_r = 1.0 if logS==false
        #     F_ev_[i,i-1] = F_ev_[i,i-1] + ws.GR*ws.scale_GR[i-1]
        # end
    end
end

function jacobian_condensation!(F_ev_::Array{Cdouble,2},ws::AeroSys,x::Array{Cdouble,1},zeta::Array{Cdouble,1}) # assume that every element of the matrix is set to zero #LATER it might be improved by not computing the addition in the function
    GR_array = (ws.GR0*ws.t0/ws.d0)*zeta[:] # CGR(zeta) # dimensionless growth rate
    ws.GR = GR_array[1]
    # for a growth rate that denpends on the size
    F_ev_[1,1]   =     F_ev_[1,1]   - GR_array[1]*ws.scale_GR[1]
    for i in 2:ws.nbin
        # derivation w.r.t. the concentration distribution
        F_ev_[i,i]   = F_ev_[i,i]   - GR_array[i]*ws.scale_GR[i]
        # F_ev_[i,i-1] = F_ev_[i,i-1] + GR_array[i]*ws.scale_GR[i]*ws.cst_r # cf CondensationGrowth! comments
        # if ws.logS
            F_ev_[i,i-1] = F_ev_[i,i-1] + GR_array[i-1]*ws.scale_GR[i-1]*ws.cst_r
        # else # cst_r = 1.0 if logS==false
        #     F_ev_[i,i-1] = F_ev_[i,i-1] + GR_array[i-1]*ws.scale_GR[i-1]
        # end
    end
end
