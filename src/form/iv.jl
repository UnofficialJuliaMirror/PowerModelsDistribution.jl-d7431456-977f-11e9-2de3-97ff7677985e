""
function variable_mc_voltage(pm::_PMs.AbstractIVRModel; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_voltage(pm, cnd=c; nw=nw, kwargs...)
    end
end

""
function variable_mc_branch_current(pm::_PMs.AbstractIVRModel; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_branch_current(pm, cnd=c; nw=nw, kwargs...)
    end
end

""
function variable_mc_dcline(pm::_PMs.AbstractIVRModel; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_dcline(pm, cnd=c; nw=nw, kwargs...)
    end
end

""
function variable_mc_load(pm::_PMs.AbstractIVRModel; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_load(pm, cnd=c; nw=nw, kwargs...)
    end
end

""
function variable_mc_gen(pm::_PMs.AbstractIVRModel; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.variable_gen(pm, cnd=c; nw=nw, kwargs...)
    end
end

""
function constraint_mc_load_power_setpoint(pm::_PMs.AbstractIVRModel, i::Int; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.constraint_load_power_setpoint(pm, i, cnd=c; nw=nw, kwargs...)
    end
end

""
function constraint_mc_gen_power_limits(pm::_PMs.AbstractIVRModel, i::Int; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.constraint_gen_power_limits(pm, i, cnd=c; nw=nw, kwargs...)
    end
end

""
function constraint_mc_current_balance(pm::_PMs.AbstractIVRModel, i::Int; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.constraint_current_balance(pm, i, cnd=c; nw=nw, kwargs...)
    end
end

""
function constraint_mc_voltage_magnitude(pm::_PMs.AbstractIVRModel, i::Int; nw=pm.cnw, kwargs...)
    for c in _PMs.conductor_ids(pm)
        _PMs.constraint_voltage_magnitude(pm, i, cnd=c; nw=nw, kwargs...)
    end
end

""
function constraint_mc_voltage_drop(pm::_PMs.AbstractIVRModel, i::Int; nw=pm.cnw, kwargs...)

    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    r = branch["br_r"].values
    x = branch["br_x"].values
    constraint_mc_voltage_drop(pm, nw, i, f_bus, t_bus, f_idx, r, x)
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::_PMs.AbstractIVRModel, n::Int, c::Int, d)
    vr = _PMs.var(pm, n, c, :vr, d)
    vi = _PMs.var(pm, n, c, :vi, d)
    @show "hallo"
    nconductors = length(_PMs.conductor_ids(pm))
    theta = _wrap_to_pi(2 * pi / nconductors * (1-c))
    # deal with cases first where tan(theta)==Inf or tan(theta)==0
    if theta == pi/2
        JuMP.@constraint(pm.model, vr == 0)
        JuMP.@constraint(pm.model, vi >= 0)
    elseif theta == -pi/2
        JuMP.@constraint(pm.model, vr == 0)
        JuMP.@constraint(pm.model, vi <= 0)
    elseif theta == 0
        JuMP.@constraint(pm.model, vr >= 0)
        JuMP.@constraint(pm.model, vi == 0)
    elseif theta == pi
        JuMP.@constraint(pm.model, vr >= 0)
        JuMP.@constraint(pm.model, vi == 0)
    else
        JuMP.@constraint(pm.model, vi == tan(theta)*vr)
        # theta also implies a sign for vi
        if 0<theta && theta < pi
            JuMP.@constraint(pm.model, vi >= 0)
        else
            JuMP.@constraint(pm.model, vi <= 0)
        end
        # theta also implies a sign for vr
        if -π/2<theta && theta < π/2
            JuMP.@constraint(pm.model, vr >= 0)
        else
            JuMP.@constraint(pm.model, vr <= 0)
        end
    end
end


"""
Defines voltage drop over a branch, linking from and to side complex voltage
"""
function constraint_mc_voltage_drop(pm::_PMs.AbstractIVRModel, n::Int, i, f_bus, t_bus, f_idx, r, x)
    Vrfr = [_PMs.var(pm, n, c, :vr, f_bus) for c in _PMs.conductor_ids(pm)]
    Vifr = [_PMs.var(pm, n, c, :vi, f_bus) for c in _PMs.conductor_ids(pm)]

    Vrto = [_PMs.var(pm, n, c, :vr, t_bus) for c in _PMs.conductor_ids(pm)]
    Vito = [_PMs.var(pm, n, c, :vi, t_bus) for c in _PMs.conductor_ids(pm)]

    Csr =  [_PMs.var(pm, n, c, :csr, f_idx) for c in _PMs.conductor_ids(pm)]
    Csi =  [_PMs.var(pm, n, c, :csi, f_idx) for c in _PMs.conductor_ids(pm)]


    JuMP.@constraint(pm.model, Vrto .== Vrfr - r*Csr + x*Csi)
    JuMP.@constraint(pm.model, Vito .== Vifr - r*Csi - x*Csr)
end

""
function solution_mc_opf_iv!(pm::_PMs.AbstractPowerModel, sol::Dict{String,<:Any})
    _PMs.add_setpoint_bus_voltage!(sol, pm)
    _PMs.add_setpoint_branch_current!(sol, pm)
    _PMs.add_setpoint_load_current!(sol, pm)
    _PMs.add_setpoint_generator_current!(sol, pm)

    _PMs.add_setpoint_dcline_flow!(sol, pm)

    # _PMs.add_dual_kcl!(sol, pm)
    # _PMs.add_dual_sm!(sol, pm) # Adds the duals of the transmission lines' thermal limits.
end
