function orbital_dynamics(::CartesianConfig, ::ImpulsiveConfig, ::VelocityConfig)
    ode_func!(dx, x, u, t) = begin
        dx[1:3] .= x[4:6]
        dx[4:6] .= -x[1:3] * μ ./ sum(x[1:3].^2)^(3/2)
    end

    return ode_func!
end


function orbital_dynamics(::CartesianConfig, ::ZeroOrderHoldConfig, ::VelocityConfig)
    ode_func!(dx, x, u, t) = begin
        dx[1:3] .= x[4:6]
        dx[4:6] .= -x[1:3] * μ ./ sum(x[1:3].^2)^(3/2) .+ u[1:3] .* thrust/maximum_mass
    end

    return ode_func!
end


function orbital_dynamics(::CartesianConfig, ::ZeroOrderHoldConfig, ::LoggedMassConfig)
    ode_func!(dx, x, u, t) = begin
        dx[1:3] .= x[4:6]
        dx[4:6] .= -x[1:3] * μ ./ sum(x[1:3].^2)^(3/2) .+ u[1:3] .* thrust
        dx[7] = -u[4]*thrust/g0_isp
    end

    return ode_func!
end


function orbital_dynamics(::CartesianConfig, ::ZeroOrderHoldConfig, ::MassConfig)
    ode_func!(dx, x, u, t) = begin
        dx[1:3] .= x[4:6]
        dx[4:6] .= -x[1:3] * μ ./ sum(x[1:3].^2)^(3/2) .+ thrust .* u[1:3] ./ x[7]
        dx[7] = -u[4]*thrust/g0_isp
    end

    return ode_func!
end


function apply_control!(u, p, Δt, index, u_nodes, ::CartesianConfig, ::ImpulsiveConfig)
    u[4:6] .+= u_nodes[:, index] * Δt
end


function apply_control!(u, p, Δt, index, u_nodes, ::CartesianConfig, ::ZeroOrderHoldConfig)
    p[1:end] .= u_nodes[:, index]
end


function integrate_arc(
    x0::Vector{T}, 
    u0::Vector{T},
    tspan;
    coordinate_config::CoordinateConfig = CartesianConfig(), 
    thrust_config::ThrustConfig = ZeroOrderHoldConfig(),
    objective_config::ObjectiveConfig = VelocityConfig(),
) where T
    x0 = copy(x0)
    p = zeros(T, 4)

    callback_index = 1

    apply_control!(
        x0, 
        p, 
        tspan[2] - tspan[1], 
        callback_index, 
        u0, 
        coordinate_config, 
        thrust_config
    )

    prob = ODEProblem(
        orbital_dynamics(coordinate_config, thrust_config, objective_config), 
        x0, 
        tspan, 
        p
    )

    return DifferentialEquations.solve(
        prob, 
        Vern9();
        abstol=1e-12, 
        reltol=1e-12, 
        save_everystep=false
    )[end]
end


function integrate_trajectory(
    x0::Vector{T}, 
    t_out; 
    t_nodes = Float64[], 
    u_nodes = zeros(0, 0), 
    coordinate_config::CoordinateConfig = CartesianConfig(), 
    thrust_config::ThrustConfig = ZeroOrderHoldConfig(),
    objective_config::ObjectiveConfig = VelocityConfig(),
) where T
    x0 = copy(x0)
    p = zeros(T, 4)

    callback_index = 1

    # Add initial control if it lines up with first node
    if size(u_nodes, 1) > 0 && t_nodes[1] ≈ 0.0 
        apply_control!(
            x0, 
            p, 
            t_nodes[callback_index + 1] - t_nodes[callback_index], 
            callback_index, 
            u_nodes, 
            coordinate_config, 
            thrust_config
        )

        callback_index += 1
    end

    affect!(integrator) = begin
        apply_control!(
            integrator.u, 
            integrator.p, 
            t_nodes[callback_index + 1] - t_nodes[callback_index], 
            callback_index, 
            u_nodes, 
            coordinate_config, 
            thrust_config
        )

        callback_index += 1
    end

    cb = PresetTimeCallback(t_nodes[callback_index:(end-1)], affect!)

    prob = ODEProblem(
        orbital_dynamics(coordinate_config, thrust_config, objective_config), 
        x0, 
        (0.0, t_out[end]), 
        p
    )

    x_nodes = stack(DifferentialEquations.solve(
        prob, 
        Vern9();
        abstol=1e-12, 
        reltol=1e-12, 
        callback=cb,
        saveat = t_out
    )(t_out))

    return x_nodes
end



function get_arc_stms(
    t_nodes, 
    x_nodes, 
    u_nodes;
    coordinate_config::CartesianConfig = CartesianConfig(),
    thrust_config::ThrustConfig = ImpulsiveConfig(),
    objective_config::ObjectiveConfig = VelocityConfig()
)
    state_size, control_size = if typeof(objective_config) == VelocityConfig
        6, 4
    else
        7, 4
    end
 
    stm_function(val, Δt) = begin
        return integrate_arc(
            val[1:state_size], 
            val[(state_size + 1):end],
            [0.0, Δt];
            coordinate_config, 
            thrust_config,
            objective_config
        )
    end

    arc_stms = zeros(state_size, state_size + control_size, size(u_nodes, 2))
    
    for i in 1:size(u_nodes, 2)
        arc_stms[:, :, i] = ForwardDiff.jacobian(
            val -> stm_function(val, t_nodes[i+1] - t_nodes[i]), 
            vcat(x_nodes[:, i], u_nodes[:, i])
        )
    end

    return arc_stms
end


function get_arc_stms_adaptive_time(
    Δt_nodes, 
    x_nodes, 
    u_nodes;
    coordinate_config::CartesianConfig = CartesianConfig(),
    thrust_config::ThrustConfig = ImpulsiveConfig(),
    objective_config::ObjectiveConfig = VelocityConfig()
)
    state_size, control_size = if typeof(objective_config) == VelocityConfig
        6, 4
    else
        7, 4
    end
 
    stm_function(val) = begin
        return integrate_arc(
            val[1:state_size], 
            val[(state_size + 1):(end-1)],
            [0.0, val[end]];
            coordinate_config, 
            thrust_config,
            objective_config
        )
    end

    arc_stms = zeros(state_size, state_size + control_size + 1, size(u_nodes, 2))
    
    for i in 1:size(u_nodes, 2)
        arc_stms[:, :, i] = ForwardDiff.jacobian(
            val -> stm_function(val), 
            vcat(x_nodes[:, i], u_nodes[:, i], Δt_nodes[i])
        )
    end

    return arc_stms
end
