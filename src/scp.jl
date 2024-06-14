

# Inputs:
# Locations, times.

# Outputs. 




struct SequentialConvexProblem{T <: Real}
    id_journey::Vector{Vector{Int64}}
    times_journey::Vector{Vector{T}}
    t_nodes::Vector{Vector{Vector{T}}}
    u_nodes::Vector{Vector{Matrix{T}}}
    x_nodes::Vector{Vector{Matrix{T}}}
    x0::Vector{Vector{Vector{T}}}
    xf::Vector{Vector{Vector{T}}}
    Δv0::Vector{Vector{Vector{T}}}
    Δvf::Vector{Vector{Vector{T}}}
    Δv0_limit::Vector{Vector{T}}
    Δvf_limit::Vector{Vector{T}}
    Δm0::Vector{Vector{T}}
    objective_config::Union{MassConfig, LoggedMassConfig}
    dynamical_error::T
    trust_region_factor::T
    optimizer
end



abstract type SequentialConvexAlgorithm end

struct SegmentedTimeFixed <: SequentialConvexAlgorithm end
struct UnifiedTimeFixed <: SequentialConvexAlgorithm end
struct MixedTimeAdaptive <: SequentialConvexAlgorithm end




function SequentialConvexProblem(
    id_journey,
    times_journey;
    objective_config = LoggedMassConfig(),
    dynamical_error = 1e-6,
    trust_region_factor = 0.1
)
    t_nodes, u_nodes, x_nodes, x0, xf, Δv0, Δvf, Δv0_limit, Δvf_limit, Δm0 = get_lambert_guess_for_scp(id_journey, times_journey)

    display(x0[1][1])

    return SequentialConvexProblem(
        id_journey,
        times_journey,
        t_nodes,
        u_nodes,
        x_nodes,
        x0,
        xf,
        Δv0,
        Δvf,
        Δv0_limit,
        Δvf_limit,
        Δm0,
        objective_config,
        dynamical_error,
        trust_region_factor,
        Mosek.Optimizer
    )
end





prob = SequentialConvexProblem(id_journey, times_journey)













function solve!(p::SequentialConvexProblem, ::MixedTimeAdaptive)

    mixing = length(p.t_nodes)

    model = Model(p.optimizer)
    set_silent(model)

    s_nodes = [[fill(1.0, length(t_nodes) - 1) for t_nodes in t_nodes] for t_nodes in p.t_nodes]

    # State variables
    x = [[] for _ in 1:mixing]

    # Control variables
    u = [[] for _ in 1:mixing]

    # Time-scaling variables
    s = [[] for _ in 1:mixing]

    # Δv variables
    Δv0 = [[] for _ in 1:mixing]
    Δvf = [[] for _ in 1:mixing]

    # State violation at start and end of each segment
    Δx0 = [[] for _ in 1:mixing]
    Δxf = [[] for _ in 1:mixing]

    # Maximum violation variable
    x_violation = [[] for _ in 1:mixing]

    # Terminal mass violation
    m_violation = [@variable(model, lower_bound = 0.0, upper_bound = 1000.0) for _ in 1:mixing]
    Δt_start = [@variable(model) for _ in 1:mixing]

    for n in 1:mixing
        segments = length(p.x0[n])

        # For each segment in each trajectory
        for k in 1:segments
            # Number of nodes in current segment
            nodes = length(p.t_nodes[n][k])

            # Variables
            push!(x[n], @variable(model, [i=1:7, j=1:nodes], start = p.x_nodes[n][k][i, j]))
            push!(u[n], @variable(model, [i=1:4, j=1:(nodes-1)], start = p.u_nodes[n][k][i, j]))
            push!(s[n], @variable(model, [i=1:(nodes-1)], start = 1.0))
            push!(Δx0[n], @variable(model, [i=1:6]))
            push!(Δxf[n], @variable(model, [i=1:6]))
            push!(x_violation[n], @variable(model, [i=1:6]))
            push!(Δv0[n], @variable(model, [i=1:3], start = p.Δv0[n][k][i]))
            push!(Δvf[n], @variable(model, [i=1:3], start = p.Δvf[n][k][i]))

            # SOC constraint for control
            @constraint(model, [i=1:(nodes-1)], [u[n][k][4, i]; u[n][k][1:3, i]] in SecondOrderCone())

            # SOC constraint for injection Δv
            @constraint(model, [p.Δv0_limit[n][k]; Δv0[n][k][1:3]] in SecondOrderCone())
            @constraint(model, [p.Δvf_limit[n][k]; Δvf[n][k][1:3]] in SecondOrderCone())

            # Ensure positivity of the violation
            @constraint(model, x_violation[n][k] .>= Δx0[n][k])
            @constraint(model, x_violation[n][k] .>= -Δx0[n][k])
            @constraint(model, x_violation[n][k] .>= Δxf[n][k])
            @constraint(model, x_violation[n][k] .>= -Δxf[n][k])

            @constraint(model, s[n][k] .== 1.0)
        end

        # Limit maximum mass at start
        if typeof(p.objective_config) == LoggedMassConfig
            @constraint(model, x[n][1][7, 1] <= log(3000/m_scale))
        else
            @constraint(model, x[n][1][7, 1] <= 3000.0/m_scale)
        end

        
    end

    # Variables for SCP convergence
    dynamic_con = [fill([], length(x)) for x in p.x0]
    delation_dynamic_con = [fill([], length(x)) for x in p.x0]
    trust_region_dynamic_con = [fill([], length(x)) for x in p.x0]
    segment_time_con = [fill([], length(x)) for x in p.x0]
    segment_start_con = [fill([], length(x)) for x in p.x0]
    segment_end_con = [fill([], length(x)) for x in p.x0]
    mass_end_con = [[] for x in p.x0]
    time_start_movement_con = [[] for x in p.x0]
    actual_time_con = [[] for x in p.x0]

    thrust_mass_con = [fill([], length(x)) for x in p.x0]
    mass_link_con = [fill([], length(x)) for x in p.x0]

    # Variables to track which segments need to be relinearised
    bad_dynamic_segments = [collect(1:length(x)) for x in p.x0]
    bad_thrust_mass_segments = [collect(1:length(x)) for x in p.x0]
    bad_mass_link_segments = [collect(2:length(x)) for x in p.x0]

    # active_parts = collect(1:mixing)

    initial_trust_region_factor = p.trust_region_factor

    for i in 1:scp_iterations
        Δt_nodes = [[t_nodes[2:end] .- t_nodes[1:(end-1)] for t_nodes in t_nodes] for t_nodes in p.t_nodes]

        objective = 0.0

        for n in 1:mixing
            segments = length(p.x0[n])

            Δt_segments = @expression(model, [i=1:segments],
                sum(s[n][i] .* Δt_nodes[n][i])
            )

            actual_time = @expression(model, [i=1:length(p.times_journey[n])], 
                p.times_journey[n][1] + sum(vcat([0.0], Δt_segments[1:(i-1)])) + Δt_start[n]
            )

            for k in bad_dynamic_segments[n]
                nodes = size(p.t_nodes[n][k], 1)

                for con in dynamic_con[n][k]
                    delete(model, con)
                end
                
                for con in delation_dynamic_con[n][k]
                    delete(model, con)
                end

                for con in trust_region_dynamic_con[n][k]
                    delete(model, con)
                end

                for con in segment_time_con[n][k]
                    delete(model, con)
                end

                for con in segment_start_con[n][k]
                    delete(model, con)
                end

                for con in segment_end_con[n][k]
                    delete(model, con)
                end

                # Get the state transition matrices from automatic differentiation
                stms = get_arc_stm_adaptive_time(Δt_nodes[n][k], p.x_nodes[n][k], p.u_nodes[n][k], CartesianConfig(), ZeroOrderHoldConfig(), p.objective_config)

                # Derived matrices from the STMs
                As = stms[:, 1:7, 1:(nodes-1)]
                Bs = stms[:, (7+1):(end-1), 1:(nodes-1)]
                Ts = stms[:, end, 1:(nodes-1)]

                # Offset vector
                cs = p.x_nodes[n][k][:, 2:end] .- stack(
                    As[:, :, j] * p.x_nodes[n][k][:, j] + 
                    Bs[:, :, j] * p.u_nodes[n][k][:, j] +
                    Ts[:, j] * Δt_nodes[n][k][j]
                    for j in 1:(nodes-1))

                dynamic_con[n][k] = @constraint(model,
                    [i=1:(nodes-1)],
                    x[n][k][:, i+1] .== As[:, :, i] * x[n][k][:, i] .+ Bs[:, :, i] * u[n][k][:, i] .+ Ts[:, i] * (Δt_nodes[n][k][i] .* s[n][k][i]) .+ cs[:, i]
                )

                delation_dynamic_con[n][k] = @constraint(model,
                    1.0 - p.trust_region_factor .<= s[n][k] .<= 1.0 + p.trust_region_factor
                )

                trust_region_dynamic_con[n][k] = @constraint(model,
                    [i=1:(nodes-1)],
                    -5e1*p.trust_region_factor - 1e-3 .<= x[n][k][1:6, i] .- p.x_nodes[n][k][1:6, i] .<= 5e1*p.trust_region_factor + 1e-3
                )

                start_position_function(time) = begin
                    ephermeris_cartesian_from_id(id_journey[n][k], time)[:]
                end 

                final_position_function(time) = begin
                    ephermeris_cartesian_from_id(id_journey[n][k+1], time)[:]
                end 

                start_jac = FiniteDiff.finite_difference_jacobian(start_position_function, [p.times_journey[n][k]])
                final_jac = FiniteDiff.finite_difference_jacobian(final_position_function, [p.times_journey[n][k + 1]])

                segment_start_con[n][k] = @constraint(model, 
                    x[n][k][1:6, 1] .== 
                        (p.x0[n][k][1:6] .+ vcat([0, 0, 0], Δv0[n][k]) .+ 
                        start_jac * (actual_time[k] - p.times_journey[n][k]) .+ Δx0[n][k])[:]
                )      

                segment_end_con[n][k] = @constraint(model, 
                    x[n][k][1:6, nodes] .== 
                        (p.xf[n][k][1:6] .- vcat([0, 0, 0], Δvf[n][k][1:3]) .+ 
                        final_jac * (actual_time[k+1] - p.times_journey[n][k + 1]) .+ Δxf[n][k])[:]
                )    
            end
            
            # Control constraints due to mass
            for k in bad_thrust_mass_segments[n]
                nodes = size(p.t_nodes[n][k], 1)

                for con in thrust_mass_con[n][k]
                    delete(model, con)
                end

                thrust_mass_con[n][k] = if typeof(p.objective_config) == LoggedMassConfig
                    @constraint(model,
                        [i=1:(nodes-1)],
                        u[n][k][4, i] <= exp.(-p.x_nodes[n][k][7, i]).*(1.0 .- x[n][k][7, i] .+ p.x_nodes[n][k][7, i])
                    )
                else
                    @constraint(model,
                        [i=1:(nodes-1)],
                        u[n][k][4, i] <= 1.0
                    )
                end
            end

            # Mass linkage constraints
            for k in bad_mass_link_segments[n]
                for con in mass_link_con[n][k]
                    delete(model, con)
                end

                mass_link_con[n][k] = if typeof(p.objective_config) == LoggedMassConfig 
                    [@constraint(model,
                        x[n][k][7, 1] == log(exp(p.x_nodes[n][k-1][7, end]) + p.Δm0[n][k]) + (1.0/(1.0 + p.Δm0[n][k]/exp(p.x_nodes[n][k-1][7, end])))*(x[n][k-1][7, end] - p.x_nodes[n][k-1][7, end])
                    )]
                else
                    [@constraint(model,
                        x[n][k][7, 1] == x[n][k-1][7, end] + p.Δm0[n][k]
                    )]
                end
            end

            for con in mass_end_con[n]
                delete(model, con)
            end

            mass_end_con[n] = if typeof(p.objective_config) == LoggedMassConfig
                [@constraint(model, x[n][end][7, end] + m_violation[n] >= log(500/m_scale - p.Δm0[n][end]))]
            else
                [@constraint(model, x[n][end][7, end] + m_violation[n] >= 500/m_scale - p.Δm0[n][end])]
            end


            for con in time_start_movement_con[n]
                delete(model[n], con)
            end

            time_start_movement_con[n] = [@constraint(model, -2e-1*p.trust_region_factor <= Δt_start[n] <= 2e-1*p.trust_region_factor)]
            
            for con in actual_time_con[n]
                delete(model[n], con)
            end

            actual_time_con[n] = [
                @constraint(model, actual_time[end] <= maximum_time - 10*p.trust_region_factor*day_scale),
                @constraint(model, actual_time[1] >= 0.0 + 10*p.trust_region_factor*day_scale)
            ]

            dropoff = p.Δm0[n] .≈ -40/m_scale
            pickup = p.Δm0[n] .> 0.0

            objective += sum(actual_time[pickup]) - sum(actual_time[dropoff]) - 1e4*sum(sum.(x_violation[n])) - 5e3*m_violation[n]
        end


        @objective(model, 
            Max, 
            objective
        )

        JuMP.optimize!(model)

        solution_summary(model)



        JuMP.value.(x[1][1])

        JuMP.value.(u[1][1])

        JuMP.value.(s[1][1])

        JuMP.value.(s[1][1])


        # Things to do

        # update times_journey
        # Update mass changes
        # Update u_nodes from control
        # Update x_nodes from control
        # Update Δt_nodes from control
        # Update t_nodes from control

        # Update Δv0, Δvf
        # Update x0, xf


        # Calculate errors 


        # Display output



        # Check for convergence


        # Update trust region

        if i >= 10
            p.trust_region_factor = initial_trust_region_factor * ((scp_iterations - i) / (scp_iterations - 10))^2.0 + 5e-5
        end 

















        




















    end















end


















function process_solution_for_scp(
    x_start_cartesian, 
    x_final_cartesian,
    Δv_start_injection,
    t_nodes, 
    u_nodes, 
    coordinate_config::CartesianConfig,
    thrust_config::ThrustConfig,
    objective_config::ObjectiveConfig
)
    x_start = copy(x_start_cartesian)
    x_final = copy(x_final_cartesian)

    x_start_injection = copy(x_start)
    x_start_injection[4:6] .+= Δv_start_injection

    x_nodes = propagate_spacecraft(x_start_injection, t_nodes; t_nodes, u_nodes, coordinate_config, thrust_config, objective_config)

    x_final_current = x_nodes[:, end]
    x_final_current_cartesian = copy(x_final_current)

    return x_start, x_final, x_nodes, x_final_current, x_final_current_cartesian
end



function solve_scp_full_segmented(
    x_departure_cartesian::Vector{Vector{T}},
    x_arrival_cartesian::Vector{Vector{T}},
    t_nodes::Vector{Vector{T}},
    u_nodes::Vector{Matrix{T}},
    Δv_departure_injection::Vector{Vector{T}},
    Δv_departure_injection_limit::Vector{T},
    Δv_arrival_injection::Vector{Vector{T}},
    Δv_arrival_injection_limit::Vector{T},
    mass_change::Vector{T},
    id_segments::Vector{Int64},
    objective_config::Union{MassConfig, LoggedMassConfig};
    coordinate_config::CoordinateConfig = CartesianConfig(),
    thrust_config::ThrustConfig = ZeroOrderHoldConfig(),
    optimizer = Mosek.Optimizer,
    linearization_error = 1e-4,
    display_output = false,
    reverse_mass = false,
    iteration_callback = nothing,
) where {T <: Real}
    x_nodes_out = []
    u_nodes_out = copy(u_nodes)
    Δv_departure_injection_out = copy(Δv_departure_injection)
    Δv_arrival_injection_out = copy(Δv_arrival_injection)
    termination = nothing
    maximum_error = nothing

    i = 1

    while i <= length(x_departure_cartesian)
        # Find the index of the next non-GA
        j = i

        while j < length(x_departure_cartesian) && id_segments[j + 1] < 0
            j += 1
        end

        if id_segments[i] == id_segments[i+1]
            i = j + 1
            continue
        end

        x_nodes_temp, u_nodes_temp, Δv_departure_injection_temp, Δv_arrival_injection_temp, termination, maximum_error = solve_scp_full(
            x_departure_cartesian[i:j],
            x_arrival_cartesian[i:j],
            t_nodes[i:j],
            u_nodes[i:j],
            Δv_departure_injection[i:j],
            Δv_departure_injection_limit[i:j],
            Δv_arrival_injection[i:j],
            Δv_arrival_injection_limit[i:j],
            mass_change[i:j],
            id_segments[i:j],
            objective_config;
            coordinate_config,
            thrust_config,
            optimizer,
            linearization_error,
            display_output,
            reverse_mass,
            fix_start_mass = true,
            iteration_callback,
        )

        if !((termination == OPTIMAL) || (termination == SLOW_PROGRESS))
            if display_output
                println("\nLEG: $(id_segments[i])->$(id_segments[j+1]) $(termination)")
            end

            break
        end

        # Set the correct starting mass for joining the segments
        if j < length(x_departure_cartesian)
            # Get the mass difference between the guess and the actual
            temp = x_nodes_temp[end][7, end] - x_departure_cartesian[j + 1][7] + mass_change[j + 1]

            for k in (j+1):length(x_departure_cartesian)
                x_departure_cartesian[k][7] += temp
            end
        end

        x_nodes_out = vcat(x_nodes_out, x_nodes_temp)
        u_nodes_out[i:j] .= u_nodes_temp
        Δv_departure_injection_out[i:j] .= Δv_departure_injection_temp
        Δv_arrival_injection_out[i:j] .= Δv_arrival_injection_temp

        i = j + 1
    end

    return x_nodes_out, u_nodes_out, Δv_departure_injection_out, Δv_arrival_injection_out, termination, maximum_error
end


function solve_scp_full(
    x_departure_cartesian::Vector{Vector{T}},
    x_arrival_cartesian::Vector{Vector{T}},
    t_nodes::Vector{Vector{T}},
    u_nodes::Vector{Matrix{T}},
    Δv_departure_injection::Vector{Vector{T}},
    Δv_departure_injection_limit::Vector{T},
    Δv_arrival_injection::Vector{Vector{T}},
    Δv_arrival_injection_limit::Vector{T},
    mass_change::Vector{T},
    id_segments::Vector{Int64},
    objective_config::Union{MassConfig, LoggedMassConfig};
    coordinate_config::CoordinateConfig = CartesianConfig(),
    thrust_config::ThrustConfig = ZeroOrderHoldConfig(),
    optimizer = Mosek.Optimizer,
    linearization_error = 1e-4,
    display_output = false,
    reverse_mass = false,
    fix_start_mass = false,
    iteration_callback = nothing,
) where {T <: Real}
    x_departure_cartesian = deepcopy(x_departure_cartesian)
    x_arrival_cartesian = deepcopy(x_arrival_cartesian)
    u_nodes = deepcopy(u_nodes)
    Δv_departure_injection = deepcopy(Δv_departure_injection)
    Δv_arrival_injection = deepcopy(Δv_arrival_injection)

    # Create the model for optimisation
    model = Model(optimizer)
    set_silent(model)
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1e-10)

    mass_total_reference = []
    x_nodes = []

    x = []
    u = []
    Δ_start = []
    Δ_end = []
    γ = []

    Δt_nodes = [t_nodes[2:end] .- t_nodes[1:(end-1)] for t_nodes in t_nodes]

    Δv_departure_injection_segments = []
    Δv_arrival_injection_segments = []
    gravity_assists = collect(1:length(x_departure_cartesian))[id_segments[1:length(x_departure_cartesian)] .< 0]
    μ_gravity_assists = μ_planets[-id_segments[gravity_assists] .- 1]
    rp_min_gravity_assists = rp_min_planets[-id_segments[gravity_assists] .- 1]

    state_size = size(x_departure_cartesian[1], 1)
    control_size = size(u_nodes[1], 1)

    for k in 1:length(x_departure_cartesian)
        # Get number of nodes for the problem
        nodes = size(t_nodes[k], 1)

        if typeof(objective_config) == LoggedMassConfig
            # Convert into the log mass for convex formulation
            x_departure_cartesian[k][7] = log(x_departure_cartesian[k][7])
            x_arrival_cartesian[k][7] = log(x_arrival_cartesian[k][7])
        end

        x_departure_segment, x_arrival_segment, x_nodes_segment, x_arrival_current_segment, x_arrival_current_cartesian_segment = process_solution_for_scp(
            x_departure_cartesian[k], 
            x_arrival_cartesian[k], 
            Δv_departure_injection[k], 
            t_nodes[k], 
            u_nodes[k], 
            coordinate_config, 
            thrust_config,
            objective_config
        )

        push!(x_nodes, x_nodes_segment)

        # State at each node
        push!(x, @variable(model, [i=1:state_size, j=1:nodes], start = x_nodes_segment[i, j]))

        # Control at each node
        push!(u, @variable(model, [i=1:control_size, j=1:(nodes-1)], start = u_nodes[k][i, j]))

        push!(Δ_start, @variable(model, [i=1:6]))
        push!(Δ_end, @variable(model, [i=1:6]))
        push!(γ, @variable(model, [i=1:6]))

        # Injection Δv for start and end nodes
        push!(Δv_departure_injection_segments, @variable(model, [i=1:3], start = Δv_departure_injection[k][i]))
        push!(Δv_arrival_injection_segments, @variable(model, [i=1:3], start = Δv_arrival_injection[k][i]))

        # SOC constraint for injection Δv
        @constraint(model, [Δv_departure_injection_limit[k]; Δv_departure_injection_segments[k][1:3]] in SecondOrderCone())
        @constraint(model, [Δv_arrival_injection_limit[k]; Δv_arrival_injection_segments[k][1:3]] in SecondOrderCone())

        # Ensure start and end node states match
        # @constraint(model, x[k][1:6, 1] .== x_departure_segment[1:6] .+ vcat([0, 0, 0], Δv_departure_injection_segments[k][1:3]) .+ Δ_start[k])
        # @constraint(model, x[k][1:6, nodes] .== x_arrival_segment[1:6] .- vcat([0, 0, 0], Δv_arrival_injection_segments[k][1:3]) .+ Δ_end[k])

        @constraint(model, x[k][1:6, 1] .== x_departure_segment[1:6] .+ vcat([0, 0, 0], Δv_departure_injection_segments[k][1:3]))
        @constraint(model, x[k][1:6, nodes] .== x_arrival_segment[1:6] .- vcat([0, 0, 0], Δv_arrival_injection_segments[k][1:3]))

        # Ensure the mass is at all times greater than the ship mass - could also add the bit here for the pickup mass
        # @constraint(model, [i=1:nodes], x[k][7, i] >= log(dry_mass))

        # Ensure positivity of the violation parameter
        @constraint(model, γ[k] .>= Δ_start[k])
        @constraint(model, γ[k] .>= -Δ_start[k])
        @constraint(model, γ[k] .>= Δ_end[k])
        @constraint(model, γ[k] .>= -Δ_end[k])

        push!(mass_total_reference, x_nodes_segment[7, :])

        # SOC constraint for control
        @constraint(model, [i=1:(nodes-1)], [u[k][4, i]; u[k][1:3, i]] in SecondOrderCone())
    end

    if reverse_mass
        @constraint(model, x[end][7, end] == x_arrival_cartesian[end][7])
        @objective(model, Min, x[1][7, 1])
    else
        if typeof(objective_config) == LoggedMassConfig
            if fix_start_mass
                @constraint(model, x[1][7, 1] == x_departure_cartesian[1][7])
            else
                @constraint(model, x[1][7, 1] <= 0.0)
            end
            
            # 
            # @constraint(model, x[1][7, 1] == log(2000/m_scale))
            # @constraint(model, x[end][7, end] >= log(500.0/m_scale))
            # @constraint(model, x[end][7, end] >= log(1200.0 / m_scale))
            # @constraint(model, x[end][7, end] >= log(1/6 - mass_change[end] + 0.05/m_scale))
        else
            if fix_start_mass
                @constraint(model, x[1][7, 1] == x_departure_cartesian[1][7])
            else
                @constraint(model, x[1][7, 1] <= 1.0)
            end

            # @constraint(model, x[1][7, 1] == x_departure_cartesian[1][7])
            # @constraint(model, x[1][7, 1] == 2000/m_scale)
            # @constraint(model, x[end][7, end] >= 500.0/m_scale)
            # @constraint(model, x[end][7, end] >= 1200.0 / m_scale)
            # @constraint(model, x[end][7, end] >= 1/6 - mass_change[end] + 0.05/m_scale)
        end

        # @objective(model, Max, x[end][7, end] - 1e3*sum(sum.(γ)) - 5e-2*sum([sum(u[4, :])/size(u, 2) for u in u])/length(u))
        # @objective(model, Max, x[end][7, end] - 1e1*sum(sum.(γ)))
        # @objective(model, Max, x[end][7, end] - 1e-3*sum([sum(u[i][4, :].*Δt_nodes[i]) for i in 1:length(u)]) - 5e1*sum(sum.(γ)))
        # @objective(model, Max, x[end][7, end] - 1e-3*sum([sum(u[i][4, :].*Δt_nodes[i]) for i in 1:length(u)]))
        # @objective(model, Max, x[end][7, end])
        @objective(model, Max, x[end][7, end] - 1e-3*sum([sum(u[i][4, :].*Δt_nodes[i]) for i in 1:length(u)])/(t_nodes[end][end] - t_nodes[1][1]))
        # @objective(model, Max, x[end][7, end] - 5e-2*sum([sum(u[4, :])/size(u, 2) for u in u])/length(u))
    end

    # Variables for SCP convergence
    dynamic_constraint = fill([], length(x_departure_cartesian))
    trust_region_dynamic_constraint = fill([], length(x_departure_cartesian))
    thrust_mass_constraint = fill([], length(x_departure_cartesian))
    mass_link_constraint = fill([], length(x_departure_cartesian))
    gravity_assist_constraint = fill([], length(gravity_assists))
    bad_dynamic_segments = collect(1:length(x_departure_cartesian))
    bad_thrust_mass_segments = collect(1:length(x_departure_cartesian))
    bad_mass_link_segments = collect(2:length(x_departure_cartesian))
    bad_gravity_assist_segments = collect(1:length(gravity_assists))

    maximum_error = 1e6

    # SCP core loop
    for i in 1:scp_iterations
        # Timer for the setup of the convex program
        t1 = time()

        # Change those which need to be relinearised
        for k in bad_dynamic_segments
            nodes = size(t_nodes[k], 1)

            for con in dynamic_constraint[k]
                delete(model, con)
            end

            # Get the state transition matrices from automatic differentiation
            stms = get_arc_state_transition(t_nodes[k], x_nodes[k], u_nodes[k], coordinate_config, thrust_config, objective_config)

            # Derived matrices from the STMs
            As = stms[:, 1:state_size, 1:(nodes-1)]
            Bs = stms[:, (state_size+1):end, 1:(nodes-1)]
            cs = x_nodes[k][:, 2:end] .- stack(As[:, :, j] * x_nodes[k][:, j] + Bs[:, :, j] * u_nodes[k][:, j] for j in 1:(nodes-1))

            # Dynamics propagated using STMs
            dynamic_constraint[k] = @constraint(model,
                [i=1:(nodes-1)],
                x[k][:, i+1] .== As[:, :, i] * x[k][:, i] .+ Bs[:, :, i] * u[k][:, i] .+ cs[:, i]
            )

            # trust_region_dynamic_constraint[k] = @constraint(model,
            #     [i=1:(nodes-1)],
            #     -0.2 .<= x[k][1:3, i] .- x_nodes[k][1:3, i] .<= 0.2
            # )
        end

        for k in bad_thrust_mass_segments
            nodes = size(t_nodes[k], 1)

            for con in thrust_mass_constraint[k]
                delete(model, con)
            end

            thrust_mass_constraint[k] = if typeof(objective_config) == LoggedMassConfig
                @constraint(model,
                    [i=1:(nodes-1)],
                    u[k][4, i] <= exp.(-mass_total_reference[k][i]).*(1.0 .- x[k][7, i] .+ mass_total_reference[k][i])
                    # u[k][4, i] <= 1.0
                )
            else
                @constraint(model,
                    [i=1:(nodes-1)],
                    u[k][4, i] <= 1.0
                )
            end
        end

        # For the mass linkage - drop off or pickup mass at a segment end
        for k in bad_mass_link_segments
            for con in mass_link_constraint[k]
                delete(model, con)
            end

            mass_link_constraint[k] = if typeof(objective_config) == LoggedMassConfig 
                [@constraint(model,
                    # x[k][7, 1] == log(max(1e-6, exp(mass_total_reference[k-1][end]) + mass_change[k])) + (1.0/(1.0 + mass_change[k]/exp(mass_total_reference[k-1][end])))*(x[k-1][7, end] - mass_total_reference[k-1][end])
                    x[k][7, 1] == log(max(exp(mass_total_reference[k-1][end]) + mass_change[k], 0.001)) + (1.0/(1.0 + mass_change[k]/exp(mass_total_reference[k-1][end])))*(x[k-1][7, end] - mass_total_reference[k-1][end])
                )]
            else
                [@constraint(model,
                    x[k][7, 1] == x[k-1][7, end] + mass_change[k]
                )]
            end
        end

        for k in bad_gravity_assist_segments
            for con in gravity_assist_constraint[k]
                delete(model, con)
            end

            maximum_deviation = get_maximum_gravity_assist_deviation(Δv_arrival_injection[gravity_assists[k] - 1], μ_gravity_assists[k], rp_min_gravity_assists[k])
            maximum_deviation = min(maximum_deviation, get_maximum_gravity_assist_deviation(Δv_departure_injection[gravity_assists[k]], μ_gravity_assists[k], rp_min_gravity_assists[k]))

            norm_arrival = norm(Δv_arrival_injection[gravity_assists[k] - 1])
            norm_departure = norm(Δv_departure_injection[gravity_assists[k]])

            norm_arrival_stm = get_norm_linearisation(Δv_arrival_injection[gravity_assists[k] - 1])
            norm_departure_stm = get_norm_linearisation(Δv_departure_injection[gravity_assists[k]])

            gravity_assist_constraint[k] = [
                @constraint(model,
                    [maximum_deviation; Δv_arrival_injection_segments[gravity_assists[k] - 1] .+ Δv_departure_injection_segments[gravity_assists[k]]] in SecondOrderCone()
                ),
                # @constraint(model,
                #     [max(norm_arrival, norm_departure); Δv_departure_injection_segments[gravity_assists[k]]] in SecondOrderCone()
                # ),

                @constraint(model,
                    norm_arrival .+ norm_arrival_stm * (Δv_arrival_injection_segments[gravity_assists[k] - 1] .- Δv_arrival_injection[gravity_assists[k] - 1]) .== norm_departure .+ norm_departure_stm * (Δv_departure_injection_segments[gravity_assists[k]] .- Δv_departure_injection[gravity_assists[k]])
                ),
            ]
        end

        # Time for convex program solution
        t2 = time()

        # Optimize the model
        JuMP.optimize!(model)

        # Retrieve the optimal solution from the model
        u_nodes = [value.(u) for u in u]
        Δv_departure_injection = [value.(Δv_departure_injection) for Δv_departure_injection in Δv_departure_injection_segments]
        Δv_arrival_injection = [value.(Δv_arrival_injection) for Δv_arrival_injection in Δv_arrival_injection_segments]

        segment_dynamics_error = zeros(Float64, length(x_departure_cartesian))
        segment_thrust_mass_error = zeros(Float64, length(x_departure_cartesian))
        segment_mass_link_error = zeros(Float64, length(x_departure_cartesian))
        segment_gravity_assist_error = zeros(Float64, length(gravity_assists))

        if !isnothing(iteration_callback)
            x_initial_departure_cartesian = locations_journey[:, 1]
            x_initial_departure_cartesian[4:6] .+= Δv_departure_injection[1]
            x_initial_departure_cartesian = vcat(x_initial_departure_cartesian, mass_starting)

            plot_ship_trajectory_low_thrust(
                t_nodes,
                u_nodes,
                x_initial_departure_cartesian,
                times_journey,
                id_journey,
                mass_changes,
                planets_classical,
                LoggedMassConfig(),
                save_location = @sprintf("figures/convex_animation/%06d.png", i),
                zoom = 1.2
            )
        end

        for k in 1:length(x_departure_cartesian)
            x_departure_cartesian[k][7] = value.(x[k][7, 1])

            # Use the optimized values to get the current solution
            x_departure_segment, x_arrival_segment, x_nodes_segment, x_arrival_current_segment, x_arrival_current_cartesian_segment = process_solution_for_scp(
                x_departure_cartesian[k], 
                x_arrival_cartesian[k], 
                Δv_departure_injection[k], 
                t_nodes[k], 
                u_nodes[k], 
                coordinate_config, 
                thrust_config,
                objective_config
            )

            # Check the maximum state error for the SCP propagation
            segment_dynamics_error[k] = maximum(abs.(x_arrival_current_cartesian_segment[1:6] .- x_arrival_cartesian[k][1:6] .+ vcat([0, 0, 0], Δv_arrival_injection[k])))
            
            if typeof(objective_config) == LoggedMassConfig
                if k > 1
                    # Check that the previous end mass 
                    segment_mass_link_error[k] = abs(exp(x_nodes_segment[7, 1]) - exp(mass_total_reference[k - 1][end]) - mass_change[k])
                end

                segment_thrust_mass_error[k] = abs(exp(mass_total_reference[k][end]) - exp(x_nodes_segment[7, end]))
                # mass_total_reference[k] = clamp.(x_nodes_segment[7, :], log(dry_mass), log(maximum_mass))
            else
                segment_thrust_mass_error[k] = abs(mass_total_reference[k][end] - x_nodes_segment[7, end])
            end

            mass_total_reference[k] = x_nodes_segment[7, :]
            x_nodes[k] = x_nodes_segment
        end

        for k in 1:length(gravity_assists)
            segment_gravity_assist_error[k] = norm(Δv_arrival_injection[gravity_assists[k] - 1]) - norm(Δv_departure_injection[gravity_assists[k]])

            # TODO better checking
            segment_gravity_assist_error[k] = max(segment_gravity_assist_error[k], 1.0)
        end

        bad_dynamic_segments = collect(1:length(x_departure_cartesian))[segment_dynamics_error .>= linearization_error]
        bad_thrust_mass_segments = collect(1:length(x_departure_cartesian))[segment_thrust_mass_error .>= linearization_error]
        bad_mass_link_segments = collect(1:length(x_departure_cartesian))[segment_mass_link_error .>= linearization_error]
        bad_gravity_assist_segments = collect(1:length(gravity_assists))[segment_gravity_assist_error .>= linearization_error]

        # TODO add gravity assist error
        maximum_error = maximum(vcat(segment_dynamics_error, segment_thrust_mass_error, segment_mass_link_error))

        end_mass = if typeof(objective_config) == LoggedMassConfig
            exp(value.(x[end][7, end]))
        else
            value.(x[end][7, end])
        end

        start_mass = if typeof(objective_config) == LoggedMassConfig
            exp(value.(x[1][7, 1]))
        else
            value.(x[1][7, 1])
        end

        # global ad_runtime += t2-t1
        # global solve_runtime += time()-t2

        if display_output
            @printf "\n%3i  " i
            
            temp = @sprintf "%10.6e  %10.5f  %10.5f  %7.6f  %7.6f  %s" maximum_error start_mass*m_scale end_mass*m_scale t2-t1 time()-t2 termination_status(model)

            if maximum_error < linearization_error
                printstyled(temp; color=:green)
            else
                printstyled(temp; color=:red)
            end
        end

        if maximum_error < linearization_error || termination_status(model) == INFEASIBLE
            break
        end
    end

    if typeof(objective_config) == LoggedMassConfig
        for k in eachindex(x_nodes)
            # Convert back into actual mass
            x_nodes[k][7, :] = exp.(x_nodes[k][7, :])
        end
    end

    u_nodes = [value.(u) for u in u_nodes]
    Δv_departure_injection = [value.(Δv_departure_injection) for Δv_departure_injection in Δv_departure_injection_segments]
    Δv_arrival_injection = [value.(Δv_arrival_injection) for Δv_arrival_injection in Δv_arrival_injection_segments]

    return x_nodes, u_nodes, Δv_departure_injection, Δv_arrival_injection, termination_status(model), maximum_error
end


function solve_scp_full_mixed_adaptive_time(
    x_departure_cartesian::Vector{Vector{Vector{T}}},
    x_arrival_cartesian::Vector{Vector{Vector{T}}},
    t_nodes::Vector{Vector{Vector{T}}},
    u_nodes::Vector{Vector{Matrix{T}}},
    Δv_departure_injection::Vector{Vector{Vector{T}}},
    Δv_departure_injection_limit::Vector{Vector{T}},
    Δv_arrival_injection::Vector{Vector{Vector{T}}},
    Δv_arrival_injection_limit::Vector{Vector{T}},
    mass_change::Vector{Vector{T}},
    id_segments::Vector{Vector{Int64}},
    objective_config::Union{MassConfig, LoggedMassConfig};
    coordinate_config::CoordinateConfig = CartesianConfig(),
    thrust_config::ThrustConfig = ZeroOrderHoldConfig(),
    optimizer = Mosek.Optimizer,
    linearization_error = 1e-4,
    display_output = false,
    reverse_mass = false,
    initial_scaling_factor = 0.1,
    time_change = false,
    iteration_callback = nothing
) where {T <: Real}
    x_departure_cartesian = deepcopy(x_departure_cartesian)
    x_arrival_cartesian = deepcopy(x_arrival_cartesian)
    u_nodes = deepcopy(u_nodes)
    t_nodes = deepcopy(t_nodes)
    Δv_departure_injection = deepcopy(Δv_departure_injection)
    Δv_arrival_injection = deepcopy(Δv_arrival_injection)

    state_size = size(x_departure_cartesian[1][1], 1)
    control_size = size(u_nodes[1][1], 1)

    model = []

    for n in 1:length(t_nodes)
        # Create the model for optimisation
        push!(model, Model(optimizer))
        set_silent(model[n])
        set_attribute(model[n], "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1e-6)
        # set_attribute(model[n], "MSK_DPAR_INTPNT_CO_TOL_PFEAS", 1e-6)
    end

    mass_total_reference = [[] for _ in t_nodes]
    x_nodes = [[] for _ in t_nodes]

    Δt_nodes = [[t_nodes[2:end] .- t_nodes[1:(end-1)] for t_nodes in t_nodes] for t_nodes in t_nodes]

    s_nodes = [[fill(1.0, length(t_nodes) - 1) for t_nodes in t_nodes] for t_nodes in t_nodes]


    x = [[] for _ in t_nodes]
    u = [[] for _ in t_nodes]
    s = [[] for _ in t_nodes]
    Δ_start = [[] for _ in t_nodes]
    Δ_end = [[] for _ in t_nodes]
    segment_time = [[] for _ in t_nodes]
    γ = [[] for _ in t_nodes]

    Δv_departure_injection_segments = [[] for _ in t_nodes]
    Δv_arrival_injection_segments = [[] for _ in t_nodes]

    start_offset = [@variable(model[n]) for n in 1:length(t_nodes)]
    mass_end_violation = [@variable(model[n]) for n in 1:length(t_nodes)]
    actual_time = []


    for n in 1:length(t_nodes)
        for k in 1:length(x_departure_cartesian[n])
            # Get number of nodes for the problem
            nodes = size(t_nodes[n][k], 1)

            if typeof(objective_config) == LoggedMassConfig
                # Convert into the log mass for convex formulation
                x_departure_cartesian[n][k][7] = log(x_departure_cartesian[n][k][7])
                x_arrival_cartesian[n][k][7] = log(x_arrival_cartesian[n][k][7])
            end

            x_departure_segment, x_arrival_segment, x_nodes_segment, x_arrival_current_segment, x_arrival_current_cartesian_segment = process_solution_for_scp(
                x_departure_cartesian[n][k], 
                x_arrival_cartesian[n][k], 
                Δv_departure_injection[n][k], 
                t_nodes[n][k], 
                u_nodes[n][k], 
                coordinate_config, 
                thrust_config,
                objective_config
            )

            push!(x_nodes[n], x_nodes_segment)

            # State at each node
            push!(x[n], @variable(model[n], [i=1:state_size, j=1:nodes], start = x_nodes_segment[i, j]))

            # Control at each node
            push!(u[n], @variable(model[n], [i=1:control_size, j=1:(nodes-1)], start = u_nodes[n][k][i, j]))

            # Time delation parameter
            push!(s[n], @variable(model[n], [i=1:(nodes-1)], start = 1.0))

            push!(Δ_start[n], @variable(model[n], [i=1:6]))
            push!(Δ_end[n], @variable(model[n], [i=1:6]))
            push!(γ[n], @variable(model[n], [i=1:6]))
            push!(segment_time[n], @variable(model[n]))

            @constraint(model[n], 0.01 >= mass_end_violation[n] >= 0.0)

            # Injection Δv for start and end nodes
            push!(Δv_departure_injection_segments[n], @variable(model[n], [i=1:3], start = Δv_departure_injection[n][k][i]))
            push!(Δv_arrival_injection_segments[n], @variable(model[n], [i=1:3], start = Δv_arrival_injection[n][k][i]))

            # SOC constraint for injection Δv
            @constraint(model[n], [Δv_departure_injection_limit[n][k]; Δv_departure_injection_segments[n][k][1:3]] in SecondOrderCone())
            @constraint(model[n], [Δv_arrival_injection_limit[n][k]; Δv_arrival_injection_segments[n][k][1:3]] in SecondOrderCone())

            # Ensure positivity of the violation parameter
            @constraint(model[n], γ[n][k] .>= Δ_start[n][k])
            @constraint(model[n], γ[n][k] .>= -Δ_start[n][k])
            @constraint(model[n], γ[n][k] .>= Δ_end[n][k])
            @constraint(model[n], γ[n][k] .>= -Δ_end[n][k])

            push!(mass_total_reference[n], x_nodes_segment[7, :])

            # SOC constraint for control
            @constraint(model[n], [i=1:(nodes-1)], [u[n][k][4, i]; u[n][k][1:3, i]] in SecondOrderCone())
        end

        if reverse_mass
            @constraint(model[n], x[n][end][7, end] == x_arrival_cartesian[n][end][7])
        else
            if typeof(objective_config) == LoggedMassConfig
                @constraint(model[n], x[n][1][7, 1] <= 0.0)
            else
                @constraint(model[n], x[n][1][7, 1] <= 1.0)
            end
            # @constraint(model[n], x[n][1][7, 1] <= x_departure_cartesian[n][1][7])
        end

        # Restrict the maximum time of the transfer
        push!(actual_time, @expression(model[n], [i=1:length(id_segments[n])], 
            times_journey_mixed[n][1] + sum(vcat([0.0], segment_time[n][1:(i-1)])) + start_offset[n]
        ))

        # @constraint(model[n], start_offset[n] == 0.0)
    end

    # Variables for SCP convergence
    dynamic_constraint = [fill([], length(x)) for x in x_departure_cartesian]
    delation_dynamic_constraint = [fill([], length(x)) for x in x_departure_cartesian]
    trust_region_dynamic_constraint = [fill([], length(x)) for x in x_departure_cartesian]
    time_total_dynamic_constraint = [fill([], length(x)) for x in x_departure_cartesian]
    segment_time_constraint = [fill([], length(x)) for x in x_departure_cartesian]
    segment_start_constraint = [fill([], length(x)) for x in x_departure_cartesian]
    segment_end_constraint = [fill([], length(x)) for x in x_departure_cartesian]
    mass_end_constraint = [[] for x in x_departure_cartesian]
    time_start_movement_constraint = [[] for x in x_departure_cartesian]
    actual_time_constraint = [[] for x in x_departure_cartesian]

    thrust_mass_constraint = [fill([], length(x)) for x in x_departure_cartesian]
    mass_link_constraint = [fill([], length(x)) for x in x_departure_cartesian]

    # Variables to track which segments need to be relinearised
    bad_dynamic_segments = [collect(1:length(x)) for x in x_departure_cartesian]
    bad_thrust_mass_segments = [collect(1:length(x)) for x in x_departure_cartesian]
    bad_mass_link_segments = [collect(2:length(x)) for x in x_departure_cartesian]

    
    start_position = [fill([], length(x)) for x in x_departure_cartesian]
    final_position = [fill([], length(x)) for x in x_departure_cartesian]

    active_parts = collect(1:length(t_nodes))
    current_maximum_error = [0.0 for _ in t_nodes]
    current_time_convergence = [0.0 for _ in t_nodes]

    maximum_error = 1e6
    scaling_factor = initial_scaling_factor

    ref_offset = [0.0 for _ in t_nodes]

    # Sequential convex programming core loop
    for i in 1:scp_iterations
        # Timer for the setup of the convex program
        t1 = time()

        for n in active_parts
            # Change those which need to be relinearised
            for k in bad_dynamic_segments[n]
                nodes = size(t_nodes[n][k], 1)

                for con in dynamic_constraint[n][k]
                    delete(model[n], con)
                end
                
                for con in delation_dynamic_constraint[n][k]
                    delete(model[n], con)
                end

                for con in trust_region_dynamic_constraint[n][k]
                    delete(model[n], con)
                end

                for con in segment_time_constraint[n][k]
                    delete(model[n], con)
                end

                for con in segment_start_constraint[n][k]
                    delete(model[n], con)
                end

                for con in segment_end_constraint[n][k]
                    delete(model[n], con)
                end

                # Get the state transition matrices from automatic differentiation
                stms = get_arc_state_transition_tmod(Δt_nodes[n][k], x_nodes[n][k], u_nodes[n][k], coordinate_config, thrust_config, objective_config)

                # Derived matrices from the STMs
                As = stms[:, 1:state_size, 1:(nodes-1)]
                Bs = stms[:, (state_size+1):(end-1), 1:(nodes-1)]
                Ts = stms[:, end, 1:(nodes-1)]

                # Offset vector
                cs = x_nodes[n][k][:, 2:end] .- stack(
                    As[:, :, j] * x_nodes[n][k][:, j] + 
                    Bs[:, :, j] * u_nodes[n][k][:, j] +
                    Ts[:, j] * Δt_nodes[n][k][j]
                    for j in 1:(nodes-1))

                # Dynamics propagated using STMs
                dynamic_constraint[n][k] = @constraint(model[n],
                    [i=1:(nodes-1)],
                    x[n][k][:, i+1] .== As[:, :, i] * x[n][k][:, i] .+ Bs[:, :, i] * u[n][k][:, i] .+ Ts[:, i] * (Δt_nodes[n][k][i] .* s[n][k][i]) .+ cs[:, i]
                )

                delation_dynamic_constraint[n][k] = @constraint(model[n],
                    1.0 - scaling_factor .<= s[n][k] .<= 1.0 + scaling_factor
                )

                trust_region_dynamic_constraint[n][k] = @constraint(model[n],
                    [i=1:(nodes-1)],
                    -5e1*scaling_factor - 1e-3 .<= x[n][k][1:6, i] .- x_nodes[n][k][1:6, i] .<= 5e1*scaling_factor + 1e-3
                )

                segment_time_constraint[n][k] = [
                    @constraint(model[n], 
                        sum(s[n][k] .* Δt_nodes[n][k]) == segment_time[n][k]
                    )
                ]

                if !time_change
                    time_total_dynamic_constraint[n][k] = [@constraint(model[n], 
                        segment_time[n][k] == t_nodes[n][k][end] - t_nodes[n][k][1]
                    )]
                end

                start_position_function(time) = begin
                    ephermeris_cartesian_from_id(id_segments[n][k], time)[:]
                end 

                final_position_function(time) = begin
                    ephermeris_cartesian_from_id(id_segments[n][k+1], time)[:]
                end 

                start_time = times_journey_mixed[n][1] + sum(vcat([0.0], [t_nodes[n][j][end] for j in 1:(k-1)])) + ref_offset[n]
                final_time = times_journey_mixed[n][1] + sum([t_nodes[n][k][end] for k in 1:k]) + ref_offset[n]

                start_position[n][k] = start_position_function(start_time)
                final_position[n][k] = final_position_function(final_time)

                start_jac = FiniteDiff.finite_difference_jacobian(start_position_function, [start_time])
                final_jac = FiniteDiff.finite_difference_jacobian(final_position_function, [final_time])

                segment_start_constraint[n][k] = @constraint(model[n], 
                    x[n][k][1:6, 1] .== 
                        (start_position[n][k] .+ vcat([0, 0, 0], Δv_departure_injection_segments[n][k][1:3]) .+ 
                        start_jac * (actual_time[n][k] - start_time) .+ Δ_start[n][k])[:]
                )      

                segment_end_constraint[n][k] = @constraint(model[n], 
                    x[n][k][1:6, nodes] .== 
                        (final_position[n][k] .- vcat([0, 0, 0], Δv_arrival_injection_segments[n][k][1:3]) .+ 
                        final_jac * (actual_time[n][k+1] - final_time) .+ Δ_end[n][k])[:]
                )    
            end


            for k in bad_thrust_mass_segments[n]
                nodes = size(t_nodes[n][k], 1)

                for con in thrust_mass_constraint[n][k]
                    delete(model[n], con)
                end

                thrust_mass_constraint[n][k] = if typeof(objective_config) == LoggedMassConfig
                    @constraint(model[n],
                        [i=1:(nodes-1)],
                        u[n][k][4, i] <= exp.(-mass_total_reference[n][k][i]).*(1.0 .- x[n][k][7, i] .+ mass_total_reference[n][k][i])
                    )
                else
                    @constraint(model[n],
                        [i=1:(nodes-1)],
                        u[n][k][4, i] <= 1.0
                    )
                end
            end

            # For the mass linkage - drop off or pickup mass at a segment end
            for k in bad_mass_link_segments[n]
                for con in mass_link_constraint[n][k]
                    delete(model[n], con)
                end

                mass_link_constraint[n][k] = if typeof(objective_config) == LoggedMassConfig 
                    [@constraint(model[n],
                        x[n][k][7, 1] == log(exp(mass_total_reference[n][k-1][end]) + mass_change[n][k]) + (1.0/(1.0 + mass_change[n][k]/exp(mass_total_reference[n][k-1][end])))*(x[n][k-1][7, end] - mass_total_reference[n][k-1][end])
                    )]
                else
                    [@constraint(model[n],
                        x[n][k][7, 1] == x[n][k-1][7, end] + mass_change[n][k]
                    )]
                end
            end

            for con in mass_end_constraint[n]
                delete(model[n], con)
            end

            if time_change
                mass_end_constraint[n] = if typeof(objective_config) == LoggedMassConfig
                    [@constraint(model[n], x[n][end][7, end] + mass_end_violation[n] >= log(1/6 - mass_change[n][end]))]
                else
                    [@constraint(model[n], x[n][end][7, end] + mass_end_violation[n] >= 1/6 - mass_change[n][end])]
                end
            end

            for con in time_start_movement_constraint[n]
                delete(model[n], con)
            end

            if time_change
                time_start_movement_constraint[n] = [@constraint(model[n], -2e-1*scaling_factor <= start_offset[n] - ref_offset[n] <= 2e-1*scaling_factor)]
            else
                time_start_movement_constraint[n] = [@constraint(model[n], start_offset[n] - ref_offset[n] == 0.0)]
            end

            for con in actual_time_constraint[n]
                delete(model[n], con)
            end

            actual_time_constraint[n] = [
                @constraint(model[n], actual_time[n][end] <= maximum_time - 10*scaling_factor*day_scale),
                @constraint(model[n], actual_time[n][1] >= 0.0 + 10*scaling_factor*day_scale)
            ]

            if !time_change
                @objective(model[n], Max, sum(x[n][end][7, end] - 100.0*sum(sum.(γ[n]))))
            else
                dropoff = mass_change[n] .≈ -0.013333333333333334
                pickup = mass_change[n] .> 0.0

                @objective(model[n], 
                    Max, 
                    sum(actual_time[n][pickup]) - sum(actual_time[n][dropoff]) - 1e4*sum(sum.(γ[n])) - 5e3*mass_end_violation[n]
                )
            end
        end

        # Time for convex program solution
        t2 = time()

        for n in active_parts
            # Optimize the model
            JuMP.optimize!(model[n])
        end

        if i >= 10
            scaling_factor = initial_scaling_factor * ((scp_iterations - i) / (scp_iterations - 10))^2.0 + 5e-5
        end

        # temp_factor = i / scp_iterations

        # scaling_factor = initial_scaling_factor * 1/(1+(temp_factor/(1-temp_factor))^(1.8)) + 5e-5

        temp = [value.(actual_time[n]) for n in 1:length(t_nodes)]

        for n in active_parts
            mass_change[n] = get_mass_change_at_ids_mixed(
                id_segments[n], 
                temp[n],
                id_segments,
                temp
            )
        end

        ref_offset = value.(start_offset)



        for n in active_parts
            current_time_convergence[n] = maximum([maximum(abs.(value.(s) .- 1.0)) for (s, s_nodes) in zip(s[n], s_nodes[n])])

            # Retrieve the optimal solution from the model
            u_nodes[n] = [value.(u) for u in u[n]]
            s_nodes[n] = [value.(s) for s in s[n]]

            t_nodes[n] = [vcat(t_nodes[1], t_nodes[1] .+ cumsum(s_nodes .* Δt_nodes)) for (t_nodes, s_nodes, Δt_nodes) in zip(t_nodes[n], s_nodes[n], Δt_nodes[n])]

            Δt_nodes[n] = [t_nodes[2:end] .- t_nodes[1:(end-1)] for t_nodes in t_nodes[n]]

            Δv_departure_injection[n] = [value.(Δv_departure_injection) for Δv_departure_injection in Δv_departure_injection_segments[n]]
            Δv_arrival_injection[n] = [value.(Δv_arrival_injection) for Δv_arrival_injection in Δv_arrival_injection_segments[n]]
            
            segment_dynamics_error = zeros(Float64, length(x_departure_cartesian[n]))
            segment_thrust_mass_error = zeros(Float64, length(x_departure_cartesian[n]))
            segment_mass_link_error = zeros(Float64, length(x_departure_cartesian[n]))

            for k in 1:length(x_departure_cartesian[n])
                x_departure_cartesian[n][k][1:6] = ephermeris_cartesian_from_id(id_segments[n][k], value(actual_time[n][k]))[:]

                x_departure_cartesian[n][k][7] = value.(x[n][k][7, 1])

                x_arrival_cartesian[n][k][1:6] = ephermeris_cartesian_from_id(id_segments[n][k+1], value(actual_time[n][k+1]))[:]

                # Use the optimized values to get the current solution
                _, _, x_nodes_segment, _, x_final_current_cartesian = process_solution_for_scp(
                    x_departure_cartesian[n][k], 
                    x_arrival_cartesian[n][k], 
                    Δv_departure_injection[n][k], 
                    t_nodes[n][k], 
                    u_nodes[n][k], 
                    coordinate_config, 
                    thrust_config,
                    objective_config
                )

                # Check the maximum state error for the SCP propagation
                segment_dynamics_error[k] = maximum(abs.(x_final_current_cartesian[1:6] .- value.(x[n][k][1:6, end])))
                
                if typeof(objective_config) == LoggedMassConfig
                    if k > 1
                        # Check that the previous end mass 
                        segment_mass_link_error[k] = abs(exp(x_nodes_segment[7, 1]) - exp(mass_total_reference[n][k - 1][end]) - mass_change[n][k])
                    end

                    segment_thrust_mass_error[k] = abs(exp(mass_total_reference[n][k][end]) - exp(x_nodes_segment[7, end]))
                    # mass_total_reference[k] = clamp.(x_nodes_segment[7, :], log(dry_mass), log(maximum_mass))
                else
                    if k > 1
                        # Check that the previous end mass 
                        segment_mass_link_error[k] = abs(x_nodes_segment[7, 1] - mass_total_reference[n][k - 1][end] - mass_change[n][k])
                    end

                    segment_thrust_mass_error[k] = abs(mass_total_reference[n][k][end] - x_nodes_segment[7, end])
                end

                mass_total_reference[n][k] = x_nodes_segment[7, :]
                x_nodes[n][k] = x_nodes_segment
            end

            # bad_dynamic_segments = collect(1:length(x_departure_cartesian))[segment_dynamics_error .>= linearization_error]
            # bad_thrust_mass_segments[n] = collect(1:length(x_departure_cartesian[n]))[segment_thrust_mass_error .>= linearization_error]
            # bad_mass_link_segments[n] = collect(1:length(x_departure_cartesian[n]))[segment_mass_link_error .>= linearization_error]

            maximum_error = maximum(vcat(segment_dynamics_error, segment_thrust_mass_error, segment_mass_link_error))

            current_maximum_error[n] = maximum_error

            if (maximum_error < linearization_error) && is_self_cleaning(id_segments[n])
                active_parts = setdiff(active_parts, [n])
            end
        end

        if !isnothing(iteration_callback)
            x_initial_departure_cartesian = ephermeris_cartesian_from_id(id_segments[1][1], temp[1][1])[:]
            x_initial_departure_cartesian[4:6] .+= Δv_departure_injection[1][1]
            x_initial_departure_cartesian = vcat(x_initial_departure_cartesian, exp.(x_nodes[1][1][7, 1]))

            repeat = if i <= 10
                5
            else
                1
            end

            for _ in 1:repeat
                plot_ship_trajectory_low_thrust(
                    t_nodes[1],
                    u_nodes[1],
                    x_initial_departure_cartesian,
                    temp[1],
                    id_journey_mixed[1],
                    mass_change[1],
                    planets_classical,
                    LoggedMassConfig(),
                    save_location = @sprintf("figures/convex_animation2/%06d.png", length(readdir("figures/convex_animation2"))),
                    zoom = 1.2,
                    iteration = i
                )
            end
        end

        for n in 1:length(t_nodes)
            end_mass = if typeof(objective_config) == LoggedMassConfig
                exp(value.(x[n][end][7, end]))
            else
                value.(x[n][end][7, end])
            end

            start_mass = if typeof(objective_config) == LoggedMassConfig
                exp(value.(x[n][1][7, 1]))
            else
                value.(x[n][1][7, 1])
            end

            if display_output
                if n == 1
                    @printf "\n%3i  " i
                else
                    @printf "\n     "
                end
                
                temp = @sprintf "%10.6e  %10.5f  %10.5f  %7.6f  %9.6f  %s" current_maximum_error[n] start_mass*m_scale end_mass*m_scale current_time_convergence[n] ref_offset[n] termination_status(model[n])

                if current_maximum_error[n] < linearization_error
                    printstyled(temp; color=:green)
                else
                    printstyled(temp; color=:red)
                end
            end
        end

        if (maximum(current_maximum_error) < linearization_error) || any(termination_status.(model) .== INFEASIBLE)
            break
        end
    end

    if typeof(objective_config) == LoggedMassConfig
        for n in 1:length(t_nodes)
            for k in eachindex(x_nodes[n])
                # Convert back into actual mass
                x_nodes[n][k][7, :] = exp.(x_nodes[n][k][7, :])
            end
        end
    end

    return x_nodes, u_nodes, t_nodes, [value.(actual_time[n]) for n in 1:length(t_nodes)], Δv_departure_injection, Δv_arrival_injection, termination_status.(model), maximum_error
end

