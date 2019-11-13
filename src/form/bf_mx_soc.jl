"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_tp_model_current(pm::_PMs.GenericPowerModel{T}, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr) where T <: SOCUBFForm
    p_fr = _PMs.var(pm, n, :P_mx)[f_idx]
    q_fr = _PMs.var(pm, n, :Q_mx)[f_idx]

    w_fr_re = _PMs.var(pm, n, :W_re)[f_bus]
    w_fr_im = _PMs.var(pm, n, :W_im)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CC_re)[i]
    ccm_im =  _PMs.var(pm, n, :CC_im)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc(pm.model, mat_real, mat_imag, complex=true)

    # reformulation_psd_rank1_to_quadratic(pm.model, w_fr_re, w_fr_im, ccm_re, ccm_im, p_s_fr, q_s_fr)

    # reformulation_psd_rank1_to_quadratic(pm.model, mat_real, mat_imag)
    # code below useful for debugging: valid inequality equired to make the SOC-NLP formulation more accurate
    # (l,i,j) = f_idx
    # t_idx = (l,j,i)
    # p_to = _PMs.var(pm, n, :P_mx)[t_idx]
    # total losses are positive when g_fr, g_to and r are positive
    # not guaranteed for individual phases though when matrix obtained through Kron's reduction
    # JuMP.@constraint(pm.model, sum(p_fr[i,i] for i in 1:size(p_fr,1)) + sum(p_to[i,i] for i in 1:size(p_to,1)) >= 0)

end


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_tp_model_current(pm::_PMs.GenericPowerModel{T}, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr) where T <: SOCConicUBFForm
    p_fr = _PMs.var(pm, n, :P_mx)[f_idx]
    q_fr = _PMs.var(pm, n, :Q_mx)[f_idx]

    w_fr_re = _PMs.var(pm, n, :W_re)[f_bus]
    w_fr_im = _PMs.var(pm, n, :W_im)[f_bus]

    ccm_re =  _PMs.var(pm, n, :CC_re)[i]
    ccm_im =  _PMs.var(pm, n, :CC_im)[i]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    mat_real = [
    w_fr_re     p_s_fr  ;
    p_s_fr'    ccm_re  ;
    ]

    mat_imag = [
    w_fr_im     q_s_fr  ;
    -q_s_fr'    ccm_im  ;
    ]

    relaxation_psd_to_soc_conic(pm.model, mat_real, mat_imag, complex=true)


    # JuMP.@constraint(pm.model,
    # [
    # w_fr_re  -w_fr_im;
    # w_fr_im   w_fr_re
    # ] in JuMP.PSDCone())
    #
    # JuMP.@constraint(pm.model,
    # [
    # ccm_re  -ccm_im;
    # ccm_im   ccm_re
    # ] in JuMP.PSDCone())

    # JuMP.@constraint(pm.model,
    # [
    # mat_real[2:4, 2:4]  -mat_imag[2:4, 2:4];
    # mat_imag[2:4, 2:4]   mat_real[2:4, 2:4]
    # ] in JuMP.PSDCone())
    #
    # JuMP.@constraint(pm.model,
    # [
    # mat_real[3:5, 3:5]  -mat_imag[3:5, 3:5];
    # mat_imag[3:5, 3:5]   mat_real[3:5, 3:5]
    # ] in JuMP.PSDCone())
    #
    # JuMP.@constraint(pm.model,
    # [
    # mat_real[[5;6;1], [5;6;1]]  -mat_imag[[5;6;1], [5;6;1]];
    # mat_imag[[5;6;1], [5;6;1]]   mat_real[[5;6;1], [5;6;1]]
    # ] in JuMP.PSDCone())

    mixre = mat_real[2:4, 2:4]
    mixim = mat_imag[2:4, 2:4]

    mixre2 = mat_real[3:5, 3:5]
    mixim2 = mat_imag[3:5, 3:5]

    mixre3 = mat_real[[5;6;1], [5;6;1]]
    mixim3 = mat_imag[[5;6;1], [5;6;1]]


    steps = pm.setting["n_steps_KK"]
    for theta= pm.setting["init_KK"]:(2*π/steps):(pm.setting["init_KK"]+2*π-0.001)
        U = [1; pm.setting["tv"]* exp(im*theta)]
        Ure = real(U)
        Uim = imag(U)
        @show U

        Are = mat_real[1:5, 1:5]
        Aim = mat_imag[1:5, 1:5]
        are = mat_real[1:5,6]
        aim = mat_imag[1:5,6]
        alphare = mat_real[6,6]


        #
        # Ure = [1 1 1 1 1]'
        # Uim = 0*[1 1 1 1 1]'
        #
        # using Combinatorics
        # n = 6
        # k = 5
        # collect(combinations(1:n,k))
        #
        Mre = w_fr_re[[1;2;3], [1;2;3]]
        Mim = w_fr_im[[1;2;3], [1;2;3]]
        Are = Mre[1:2, 1:2]
        Aim = Mim[1:2, 1:2]
        are = Mre[1:2,3]
        aim = Mim[1:2,3]
        alphare = Mre[3,3]
        sdp_to_soc_kim_kojima(pm.model, Are, Aim, are, aim, alphare, Ure, Uim)

        Mre = w_fr_re[[2;3;1], [2;3;1]]
        Mim = w_fr_im[[2;3;1], [2;3;1]]
        Are = Mre[1:2, 1:2]
        Aim = Mim[1:2, 1:2]
        are = Mre[1:2,3]
        aim = Mim[1:2,3]
        alphare = Mre[3,3]
        sdp_to_soc_kim_kojima(pm.model, Are, Aim, are, aim, alphare, Ure, Uim)


        Mre = w_fr_re[[3;1;2], [3;1;2]]
        Mim = w_fr_im[[3;1;2], [3;1;2]]
        Are = Mre[1:2, 1:2]
        Aim = Mim[1:2, 1:2]
        are = Mre[1:2,3]
        aim = Mim[1:2,3]
        alphare = Mre[3,3]
        sdp_to_soc_kim_kojima(pm.model, Are, Aim, are, aim, alphare, Ure, Uim)

        Mre = ccm_re[[1;2;3], [1;2;3]]
        Mim = ccm_im[[1;2;3], [1;2;3]]
        Are = Mre[1:2, 1:2]
        Aim = Mim[1:2, 1:2]
        are = Mre[1:2,3]
        aim = Mim[1:2,3]
        alphare = Mre[3,3]
        sdp_to_soc_kim_kojima(pm.model, Are, Aim, are, aim, alphare, Ure, Uim)

        Mre = ccm_re[[2;3;1], [2;3;1]]
        Mim = ccm_im[[2;3;1], [2;3;1]]
        Are = Mre[1:2, 1:2]
        Aim = Mim[1:2, 1:2]
        are = Mre[1:2,3]
        aim = Mim[1:2,3]
        alphare = Mre[3,3]
        sdp_to_soc_kim_kojima(pm.model, Are, Aim, are, aim, alphare, Ure, Uim)


        Mre = ccm_re[[3;1;2], [3;1;2]]
        Mim = ccm_im[[3;1;2], [3;1;2]]
        Are = Mre[1:2, 1:2]
        Aim = Mim[1:2, 1:2]
        are = Mre[1:2,3]
        aim = Mim[1:2,3]
        alphare = Mre[3,3]
        sdp_to_soc_kim_kojima(pm.model, Are, Aim, are, aim, alphare, Ure, Uim)
        # Ure = [1 -1 1 -1 1]'
        # Uim = 0*[1 1 1 1 1]'
        # sdp_to_soc_kim_kojima(pm.model, w_fr_re,    w_fr_im,    theta)
        #
        # sdp_to_soc_kim_kojima(pm.model, mixre,    mixim,    theta)
        # sdp_to_soc_kim_kojima(pm.model, mixre2,    mixim2,    theta)
        # sdp_to_soc_kim_kojima(pm.model, mixre3,    mixim3,    theta)

        # sdp_to_soc_kim_kojima(pm.model, ccm_re,     ccm_im,     theta)
    end

end
