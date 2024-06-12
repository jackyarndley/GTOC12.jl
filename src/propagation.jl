function orbital_dynamics(::CartesianConfig, ::ImpulsiveConfig, ::VelocityConfig)
    ode_func!(dx, x, u, t) = begin
        dx[1:3] .= x[4:6]
        dx[4:6] .= -x[1:3] * μ ./ sum(x[1:3].^2)^(3/2)
    end

    return ode_func!
end

function orbital_dynamics(::ModifiedEquinoctialConfig, ::ImpulsiveConfig, ::VelocityConfig)
    ode_func!(dx, x, u, t) = begin
        p, f, g, h, k, L = x
        w = 1.0 + f*cos(L) + g*sin(L)

        dx[1] = 0.0
        dx[2] = 0.0
        dx[3] = 0.0
        dx[4] = 0.0
        dx[5] = 0.0
        dx[6] = sqrt(μ*p)*(w/p)^2
    end

    return ode_func!
end

function orbital_dynamics(::ModifiedOrbitalConfig, ::ImpulsiveConfig, ::VelocityConfig)
    ode_func!(dx, x, u, t) = begin
        Λ, η, s, γ, κ, β, t = x
    
        dx[1] = -η
        dx[2] = Λ
        dx[3] = γ
        dx[4] = -s
        dx[5] = 0
        dx[6] = 0
        dx[7] = 1 / (κ*(κ + Λ)^2)
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


function orbital_dynamics(::CartesianConfig, ::ZeroOrderHoldConfig, ::FinalConfig)
    ode_func!(dx, x, u, t) = begin
        dx[1:3] .= x[4:6]
        dx[4:6] .= -x[1:3] * μ ./ sum(x[1:3].^2)^(3/2) .+ thrust .* u[1:3] ./ x[7]
        dx[7] = -norm(u[1:3])*thrust/g0_isp
    end

    return ode_func!
end

function apply_control(x, u, Δt, index, u_nodes, coordinate_config::CartesianConfig, thrust_config::ImpulsiveConfig)
    x[4:6] .+= u_nodes[:, index] * Δt
    return x, u
end

function apply_control(x, u, Δt, index, u_nodes, coordinate_config::CartesianConfig, thrust_config::ZeroOrderHoldConfig)
    u[1:size(u_nodes, 1)] .= u_nodes[:, index]
    return x, u
end


"""
Propagate in dynamics with manoeuvres specified
by the nodes. If t out matches with t nodes, the out velocity is provided.
Manoeuvre Δv are at each node but not the last one.
Integrated time starts from t = 0.
"""
function propagate_spacecraft(
    x_start::Vector{T}, 
    t_out; 
    t_nodes = Float64[], 
    u_nodes = zeros(0, 0), 
    coordinate_config::CoordinateConfig = CartesianConfig(), 
    thrust_config::ThrustConfig = ImpulsiveConfig(),
    objective_config::ObjectiveConfig = VelocityConfig(),
    extra_callback = nothing
) where T
    # Check if gain needs to be applied
    apply_initial_control = size(u_nodes, 1) > 0 && t_nodes[1] ≈ 0.0 
    x_start = copy(x_start)
    u_start = zeros(T, 7)
    callback_index = 1

    Δt = if size(t_nodes, 1) > 1
        t_nodes[2] - t_nodes[1]
    else
        t_out[1]
    end

    affect!(integrator) = begin
        integrator.u, integrator.p = apply_control(
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

    # Add initial Δv if it lines up with first node
    cb = if apply_initial_control
        x_start, u_start = apply_control(
            x_start, 
            u_start, 
            Δt, 
            callback_index, 
            u_nodes, 
            coordinate_config, 
            thrust_config
        )

        callback_index += 1
        PresetTimeCallback(t_nodes[2:(end-1)], affect!)
    else
        PresetTimeCallback(t_nodes[1:(end-1)], affect!)
    end

    if !isnothing(extra_callback)
        cb = CallbackSet(cb, extra_callback)
    end

    prob = ODEProblem(
        orbital_dynamics(coordinate_config, thrust_config, objective_config), 
        x_start, 
        (0.0, t_out[end]), 
        u_start
    )

    x_nodes = stack(DifferentialEquations.solve(
        prob, 
        Vern9();
        abstol=1e-12, 
        reltol=1e-12, 
        callback=cb
    )(t_out))

    # Remove the initial Δv that was added if it lines up with the first node
    if apply_initial_control && t_out[1] ≈ 0.0
        x_nodes[:, 1], _ = apply_control(
            x_start, 
            u_start, 
            Δt, 
            1, 
            -u_nodes[:, 1:1], 
            coordinate_config, 
            thrust_config
        )
    end

    return x_nodes
end

"""
Get the spacecraft STM for an arc with Δt specified.
"""
function get_arc_stm(
    t_nodes, 
    x_nodes, 
    u_nodes,
    coordinate_config::CartesianConfig = CartesianConfig(),
    thrust_config::ThrustConfig = ImpulsiveConfig(),
    objective_config::ObjectiveConfig = VelocityConfig()
)
    thrust_nodes = size(u_nodes, 2)

    state_size, control_size = if typeof(objective_config) == VelocityConfig
        6, 4
    else
        7, 4
    end
 
    stm_function(values, Δt) = begin
        x_start = values[1:state_size]
        u_start = values[(state_size + 1):end][:, :]

        return propagate_spacecraft(
            x_start, 
            [Δt]; 
            t_nodes = [0.0], 
            u_nodes = u_start, 
            coordinate_config, 
            thrust_config,
            objective_config
        )
    end

    arc_stms = zeros(state_size, state_size + control_size, thrust_nodes)
    
    for i in 1:thrust_nodes
        arc_stms[:, :, i] = ForwardDiff.jacobian(
            value -> stm_function(value, t_nodes[i+1] - t_nodes[i]), 
            vcat(x_nodes[:, i], u_nodes[:, i])
        )
    end

    return arc_stms
end

"""
Get the spacecraft STM for an arc with Δt specified.
"""
function get_arc_stm_adaptive_time(
    Δt_nodes, 
    x_nodes, 
    u_nodes,
    coordinate_config::CartesianConfig = CartesianConfig(),
    thrust_config::ThrustConfig = ImpulsiveConfig(),
    objective_config::ObjectiveConfig = VelocityConfig()
)
    thrust_nodes = size(u_nodes, 2)

    state_size, control_size = if typeof(objective_config) == VelocityConfig
        6, 4
    else
        7, 4
    end
 
    stm_function(values) = begin
        x_start = values[1:state_size]
        u_start = values[(state_size + 1):(end - 1)][:, :]

        return propagate_spacecraft(
            x_start, 
            [values[end]]; 
            t_nodes = [0.0], 
            u_nodes = u_start, 
            coordinate_config, 
            thrust_config,
            objective_config
        )
    end

    arc_stms = zeros(state_size, state_size + control_size + 1, thrust_nodes)
    
    for i in 1:thrust_nodes
        arc_stms[:, :, i] = ForwardDiff.jacobian(
            value -> stm_function(value), 
            vcat(x_nodes[:, i], u_nodes[:, i], Δt_nodes[i])
        )
    end

    return arc_stms
end
