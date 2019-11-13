"""
SDP to SOC relaxation of type 2, applied to real-value matrix,  as described in:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_real(m, mx)
    @assert size(mx,1) == size(mx,2)
    n_elements = size(mx,1)
    for i in 1:n_elements-1
        for j in i+1:n_elements
            JuMP.@constraint(m, mx[i,j]^2 <= mx[i,i]*mx[j,j])
        end
    end
end


"""
SDP to SOC relaxation of type 2, applied to complex-value matrix,  as described in:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_complex(m, mxreal, mximag)
    @assert size(mxreal) == size(mximag)
    n_elements = size(mxreal,1)
    for i in 1:n_elements-1
        for j in i+1:n_elements
            JuMP.@constraint(m, mxreal[i,j]^2 + mximag[i,j]^2 <= mxreal[i,i]*mxreal[j,j])
        end
    end
end


"""
SDP to SOC relaxation of type 2, applied to real-value matrix,  as described in:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_real_conic(m, mx)
    @assert size(mx,1) == size(mx,2)
    n_elements = size(mx,1)
    for i in 1:n_elements-1
        for j in i+1:n_elements
            JuMP.@constraint(m, [mx[i,i]+mx[j,j], 2*mx[i,j], mx[i,i]-mx[j,j]] in JuMP.SecondOrderCone())
        end
    end
end


"""
SDP to SOC relaxation of type 2, applied to complex-value matrix,  as described in:
```
@article{Kim2003,
author = {Kim, S and Kojima, M and Yamashita, M},
title = {{Second order cone programming relaxation of a positive semidefinite constraint}},
doi = {10.1080/1055678031000148696},
journal = {Optimization Methods and Software},
number = {5},
pages = {535--541},
volume = {18},
year = {2003}
}
```
"""
function relaxation_psd_to_soc_complex_conic(m, mxreal, mximag)
    @assert size(mxreal) == size(mximag)
    n_elements = size(mxreal,1)
    for i in 1:n_elements-1
        for j in i+1:n_elements
            JuMP.@constraint(m, [mxreal[i,i]+mxreal[j,j], 2*mxreal[i,j], 2*mximag[i,j], mxreal[i,i]-mxreal[j,j]] in JuMP.SecondOrderCone())
        end
    end
end


"""
See section 4.3 for complex to real PSD constraint transformation:
@article{Fazel2001,
author = {Fazel, M. and Hindi, H. and Boyd, S.P.},
title = {{A rank minimization heuristic with application to minimum order system approximation}},
doi = {10.1109/ACC.2001.945730},
journal = {Proc. American Control Conf.},
number = {2},
pages = {4734--4739},
url = {http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=945730},
volume = {6},
year = {2001}
}
"""
function relaxation_psd_to_soc(m, mxreal, mximag; complex=true)
    if complex==false
        @assert size(mxreal) == size(mximag)
        mx =
            [
            mxreal -mximag;
            mximag  mxreal
            ]

        relaxation_psd_to_soc_real(m, mx)
    else
        relaxation_psd_to_soc_complex(m, mxreal, mximag)
    end
end


"""
See section 4.3 for complex to real PSD constraint transformation:
@article{Fazel2001,
author = {Fazel, M. and Hindi, H. and Boyd, S.P.},
title = {{A rank minimization heuristic with application to minimum order system approximation}},
doi = {10.1109/ACC.2001.945730},
journal = {Proc. American Control Conf.},
number = {2},
pages = {4734--4739},
url = {http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=945730},
volume = {6},
year = {2001}
}
"""
function relaxation_psd_to_soc_conic(m, mxreal, mximag; complex=true)
    if complex==false
        @assert size(mxreal) == size(mximag)
        mx =
            [
            mxreal -mximag;
            mximag  mxreal
            ]

        relaxation_psd_to_soc_real_conic(m, mx)
    else
        relaxation_psd_to_soc_complex_conic(m, mxreal, mximag)
    end
end


"""
For debugging / exploration: real-valued SDP to SDP relaxation based on PSDness of principal minors, default is 3x3 SDP relaxation
"""
function relaxation_psd_to_psd_real(m, mxreal, mximag; ndim=3)
    @assert size(mxreal) == size(mximag)
    @assert size(mxreal,1) >= ndim
    n_elements = size(mxreal,1)
    for i in 1:n_elements-(ndim-1)
        j = i+(ndim-1)
        mr = mxreal[i:j, i:j]
        mi = mximag[i:j, i:j]
        JuMP.@constraint(m, [mr -mi; mi mr] in JuMP.PSDCone())
    end
end

"""
This constraints models
M = [ Ar+im*Ai    Cr+im*Ci
     (Cr+im*Ci)'  Br+im*Bi]
M in PSDCone, rank(M) = 1
as a set of polynomial matrix equations
(Cr+im*Ci) (Cr+im*Ci)' = tr(Br)(Ar+im*Ai)
(Cr+im*Ci)'(Cr+im*Ci)  = tr(Ar)(Br+im*Bi)
"""
function reformulation_psd_rank1_to_quadratic(model, Ar, Ai, Br, Bi, Cr, Ci)
    @assert size(Ar) == size(Ai)
    @assert size(Br) == size(Bi)
    @assert size(Cr) == size(Ci)
    @assert size(Cr,1) == size(Ar,1)
    @assert size(Cr,2) == size(Br,2)

    # n = size(Ar,1)
    # m = size(Br,1)
    #
    # utridiagA   = [(i,j) for i=1:n, j=1:n if i<=j] # upper triangle + diagonal elements of A
    # utriA       = [(i,j) for i=1:n, j=1:n if i< j] # upper triangle elements of A
    #
    # utridiagB   = [(i,j) for i=1:m, j=1:m if i<=j] # upper triangle + diagonal elements of B
    # utriB       = [(i,j) for i=1:m, j=1:m if i< j] # upper triangle elements of B
    #
    # for (a,b) in utridiagA
    #     JuMP.@constraint(model, sum(Cr[a,j]*Cr[b,j] + Ci[a,j]*Ci[b,j] for j in 1:m) == Ar[a,b] * sum(Br[j,j] for j in 1:m))
    # end
    #
    # for (a,b) in utriA
    #     JuMP.@constraint(model, sum(Ci[a,j]*Cr[b,j] - Cr[a,j]*Ci[b,j] for j in 1:m) == Ai[a,b] * sum(Br[j,j] for j in 1:m))
    # end
    #
    # for (a,b) in utridiagB
    #     JuMP.@constraint(model, sum(Cr[i,a]*Cr[i,b] + Ci[i,a]*Ci[i,b] for i in 1:n) == Br[a,b] * sum(Ar[i,i] for i in 1:n))
    # end
    #
    # for (a,b) in utriB
    #     JuMP.@constraint(model, sum(Ci[i,a]*Cr[i,b] - Cr[i,a]*Ci[i,b] for i in 1:n) == Bi[a,b] * sum(Ar[i,i] for i in 1:n))
    # end


    # (Cr+im*Ci) (Cr+im*Ci)' = tr(Br)(Ar+im*Ai)
    matrix_product_real(model, Ar, Br, Cr, Ci)
    matrix_product_imag(model, Ai, Br, Cr, Ci)

    # equations are redundant when size(Cr) == (1,1)
    if size(Cr) != (1, 1)
        # (Cr+im*Ci)'(Cr+im*Ci)  = tr(Ar)(Br+im*Bi)
        # swap A and B, take complex conjugate transpose of C
        matrix_product_real(model, Br, Ar, Cr', -Ci')
        matrix_product_imag(model, Bi, Ar, Cr', -Ci')
    end
end

"""
This constraints models
M = Mr + im*Mi in PSDCone, rank(M) = 1
as a set of quadratic matrix equations by partitioning M in four equal quadrants
Ar+im*Ai the upper diagonal block
Br+im*Bi the lower diagonal block
Cr+im*Ci the upper right block

and then calls reformulation_psd_rank1_to_quadratic(model, Ar, Ai, Br, Bi, Cr, Ci)
"""
function reformulation_psd_rank1_to_quadratic(model, Mr, Mi)
    @assert size(Mr) == size(Mi)
    @assert size(Mr,1) == size(Mr,2)
    @assert iseven(size(Mr,1))

    n = Int(size(Mr,1)/2)

    Ar = Mr[1:n,    1:n]
    Ai = Mi[1:n,    1:n]
    Br = Mr[n+1:2n, n+1:2n]
    Bi = Mi[n+1:2n, n+1:2n]
    Cr = Mr[1:n,    n+1:2n]
    Ci = Mi[1:n,    n+1:2n]

    reformulation_psd_rank1_to_quadratic(model, Ar, Ai, Br, Bi, Cr, Ci)
end

"""
This constraints models the matrix equation
(Cr)(Cr)' + (Ci)(Ci)' = tr(Br)(Ar)
as a set of scalar equations.
The lower triangular elements generate redundant constraints due to symmetry, and are therefore not constructed
"""
function matrix_product_real(model, Ar, Br, Cr, Ci)
    n = size(Cr,1)
    m = size(Cr,2)
    utridiag   = [(i,j) for i=1:n, j=1:n if i<=j] # upper triangle + diagonal elements of A

    for (a,b) in utridiag
        JuMP.@constraint(model, sum(Cr[a,j]*Cr[b,j] + Ci[a,j]*Ci[b,j] for j in 1:m) == Ar[a,b] * sum(Br[j,j] for j in 1:m))
    end
end

"""
This constraints models the matrix equation
(Ci)(Cr)' - (Cr)(Ci)' = tr(Br)(Ai)
as a set of scalar equations.
The lower triangular and diagonal elements generate redundant constraints due to symmetry, and are therefore not constructed
"""
function matrix_product_imag(model, Ai, Br, Cr, Ci)
    n = size(Cr,1)
    m = size(Cr,2)
    utri       = [(i,j) for i=1:n, j=1:n if i< j] # upper triangle elements of A

    for (a,b) in utri
        JuMP.@constraint(model, sum(Ci[a,j]*Cr[b,j] - Cr[a,j]*Ci[b,j] for j in 1:m) == Ai[a,b] * sum(Br[j,j] for j in 1:m))
    end
end


function sdp_to_soc_kim_kojima(model, Are, Aim, are, aim, alphare, Ure, Uim; tol=1e-8)
    @assert size(Are) == size(Aim)
    @assert size(are) == size(aim)

    @assert size(Are)[1] == size(are)[1]

    Ure[abs.(Ure).<=tol] .=0
    Uim[abs.(Uim).<=tol] .=0

    U = Ure+im*Uim
    C = U*U'
    Cre = real(C)
    Cim = imag(C)
    Cre[abs.(Cre).<=tol] .=0
    Cim[abs.(Cim).<=tol] .=0

    @assert size(Cre) == size(Are)

    rhs_1 = alphare
    rhs_2 = sum(Cre.*Are) + sum(Cim.*Aim)

    lhs_re = Ure'* are + Uim'* aim
    lhs_im = Ure'* aim - Uim'* are
    @show rhs_1, rhs_2
    @show lhs_re, lhs_im

    JuMP.@constraint(model, rhs_2 >= 0)
    JuMP.@constraint(model,     [rhs_1+rhs_2;
                                 rhs_1-rhs_2;
                                 2*lhs_re;
                                 2*lhs_im] in JuMP.SecondOrderCone())

end

#
# function sdp_to_soc_kim_kojima(model, Mre, Mim, theta)
#     @assert size(Mre) == size(Mim) == (3,3)
#
#     Are = Mre[[1;2;3], [1;2;3]]
#     Aim = Mim[[1;2;3], [1;2;3]]
#
#     # sdp_to_soc_kim_kojima(model, Are, Aim, are, aim, alphare, Ure, Uim)
#
#     CA = Are[1,1] + Are[2,2] + 2*(cos(theta)*Are[1,2] - sin(theta)*Aim[1,2] )
#     UHare = Are[1,3] + cos(theta)*Are[2,3] + sin(theta)*Aim[2,3]
#     UHaim = Aim[1,3] + cos(theta)*Aim[2,3] - sin(theta)*Are[2,3]
#
#     JuMP.@constraint(model,     [Are[3,3] + CA;
#                                  Are[3,3] - CA;
#                                  2*UHare;
#                                  2*UHaim] in JuMP.SecondOrderCone())
#
#
#     Are = Mre[[2;3;1], [2;3;1]]
#     Aim = Mim[[2;3;1], [2;3;1]]
#     # sdp_to_soc_kim_kojima(model, Are, Aim, are, aim, alphare, Ure, Uim)
#
#     CA = Are[1,1] + Are[2,2] + 2*(cos(theta)*Are[1,2] - sin(theta)*Aim[1,2] )
#     UHare = Are[1,3] + cos(theta)*Are[2,3] + sin(theta)*Aim[2,3]
#     UHaim = Aim[1,3] + cos(theta)*Aim[2,3] - sin(theta)*Are[2,3]
#
#     JuMP.@constraint(model,     [Are[3,3] + CA;
#                                  Are[3,3] - CA;
#                                  2*UHare;
#                                  2*UHaim] in JuMP.SecondOrderCone())
#
#
#     Are = Mre[[3;1;2], [3;1;2]]
#     Aim = Mim[[3;1;2], [3;1;2]]
#     # sdp_to_soc_kim_kojima(model, Are, Aim, are, aim, alphare, Ure, Uim)
#     CA = Are[1,1] + Are[2,2] + 2*(cos(theta)*Are[1,2] - sin(theta)*Aim[1,2] )
#     UHare = Are[1,3] + cos(theta)*Are[2,3] + sin(theta)*Aim[2,3]
#     UHaim = Aim[1,3] + cos(theta)*Aim[2,3] - sin(theta)*Are[2,3]
#
#     JuMP.@constraint(model,     [Are[3,3] + CA;
#                                  Are[3,3] - CA;
#                                  2*UHare;
#                                  2*UHaim] in JuMP.SecondOrderCone())
# end
