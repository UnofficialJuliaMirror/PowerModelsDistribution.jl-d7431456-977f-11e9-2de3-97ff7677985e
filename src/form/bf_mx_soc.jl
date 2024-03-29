"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCUBFModels, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
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

    # code below useful for debugging: valid inequality equired to make the SOC-NLP formulation more accurate
    # (l,i,j) = f_idx
    # t_idx = (l,j,i)
    # p_to = _PMs.var(pm, n, :P_mx)[t_idx]
    # total losses are positive when g_fr, g_to and r are positive
    # not guaranteed for individual phases though when matrix obtained through Kron's reduction
    # JuMP.@constraint(pm.model, tr(p_fr) + tr(p_to) >= 0)
end


"Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude"
function constraint_mc_model_current(pm::SOCConicUBFModel, n::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr)
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
end
