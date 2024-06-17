

mutable struct SequentialConvexProblem{T <: Real}
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
    trust_region_factor = 0.1,
)
    t_nodes, u_nodes, x_nodes, x0, xf, Δv0, Δvf, Δv0_limit, Δvf_limit, Δm0 = get_lambert_guess_for_scp(id_journey, times_journey; objective_config)

    return SequentialConvexProblem(
        deepcopy(id_journey),
        deepcopy(times_journey),
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


function solve!(
    p::SequentialConvexProblem, 
    ::MixedTimeAdaptive;
    adaptive_time = true
)
    mixing = length(p.t_nodes)

    model = Model(p.optimizer)
    set_silent(model)
    set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1e-4)

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

            if !adaptive_time
                @constraint(model, s[n][k] .== 1.0)
            end
        end

        if !adaptive_time
            @constraint(model, Δt_start[n] .== 0.0)
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

        Δt_segments = [@expression(model, [i=1:length(p.x0[n])],
            sum(s[n][i] .* Δt_nodes[n][i])
        ) for n in 1:mixing]

        actual_time = [@expression(model, [i=1:length(p.times_journey[n])], 
            p.times_journey[n][1] + sum(vcat([0.0], Δt_segments[n][1:(i-1)])) + Δt_start[n]
        ) for n in 1:mixing]

        objective = 0.0

        for n in 1:mixing
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
                        start_jac * (actual_time[n][k] - p.times_journey[n][k]) .+ Δx0[n][k])[:]
                )      

                segment_end_con[n][k] = @constraint(model, 
                    x[n][k][1:6, nodes] .== 
                        (p.xf[n][k][1:6] .- vcat([0, 0, 0], Δvf[n][k][1:3]) .+ 
                        final_jac * (actual_time[n][k+1] - p.times_journey[n][k + 1]) .+ Δxf[n][k])[:]
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
                delete(model, con)
            end

            time_start_movement_con[n] = [@constraint(model, -2e-1*p.trust_region_factor <= Δt_start[n] <= 2e-1*p.trust_region_factor)]
            
            for con in actual_time_con[n]
                delete(model, con)
            end

            actual_time_con[n] = [
                @constraint(model, actual_time[n][end] <= maximum_time - 10*p.trust_region_factor*day_scale),
                @constraint(model, actual_time[n][1] >= 0.0 + 10*p.trust_region_factor*day_scale)
            ]

            dropoff = p.Δm0[n] .≈ -40/m_scale
            pickup = p.Δm0[n] .> 0.0

            if !adaptive_time
                # objective += x[n][end][7, end]
                objective -= x[n][1][7, 1]
            else
                objective += sum(actual_time[n][pickup]) - sum(actual_time[n][dropoff])
            end

            objective -= 5e3*m_violation[n]
            objective -= 1e4*sum(sum.(x_violation[n]))
            # objective -= 1e-2*sum([sum(u[n][i][4, :].*Δt_nodes[n][i]) for i in 1:length(u[n])])
        end

        @objective(model, 
            Max, 
            objective
        )
        
        JuMP.optimize!(model)

        for n in 1:mixing
            p.times_journey[n] = JuMP.value.(actual_time[n])
            p.u_nodes[n] = [JuMP.value.(u) for u in u[n]]
            p.Δv0[n] = [JuMP.value.(Δv0) for Δv0 in Δv0[n]]
            p.Δvf[n] = [JuMP.value.(Δvf) for Δvf in Δvf[n]]

            p.t_nodes[n] = [vcat(0.0, cumsum(JuMP.value.(s[n][k]) .* Δt_nodes[n][k])) for k in 1:length(p.x0[n])]
        end

        for n in 1:mixing
            dynamical_errors = []
            mass_link_errors = []

            p.Δm0[n] = get_mass_change_at_ids_mixed(
                p.id_journey[n],
                p.times_journey[n],
                p.id_journey,
                p.times_journey
            )

            for k in 1:length(p.x0[n])
                p.x0[n][k][1:6] = ephermeris_cartesian_from_id(
                    p.id_journey[n][k], 
                    p.times_journey[n][k]
                )[:]
                p.x0[n][k][7] = value.(x[n][k][7, 1])

                p.xf[n][k][1:6] = ephermeris_cartesian_from_id(
                    p.id_journey[n][k+1], 
                    p.times_journey[n][k+1]
                )[:]
                p.xf[n][k][7] = value.(x[n][k][7, end])

                p.x_nodes[n][k] = propagate_spacecraft(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    p.t_nodes[n][k];
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                push!(
                    dynamical_errors, 
                    maximum(abs.(p.xf[n][k] .- p.x_nodes[n][k][:, end] .- vcat([0.0, 0.0, 0.0], p.Δvf[n][k], [0.0])))
                )

                if k > 1
                    push!(
                        mass_link_errors,
                        abs(exp(p.x_nodes[n][k-1][7, end]) + p.Δm0[n][k] - exp(p.x_nodes[n][k][7, 1]))
                    )
                end
            end

            current_maximum_error = maximum(vcat(dynamical_errors, mass_link_errors))

            start_mass = if typeof(p.objective_config) == LoggedMassConfig
                exp(p.x_nodes[n][1][7, 1])
            else
                p.x_nodes[n][1][7, 1]
            end

            end_mass = if typeof(p.objective_config) == LoggedMassConfig
                exp(p.x_nodes[n][end][7, end])
            else
                p.x_nodes[n][end][7, end]
            end
            
            if n == 1
                @printf "\n%3i  " i
            else
                @printf "\n     "
            end
            
            temp = @sprintf "%10.6e  %10.5f  %10.5f  %10.5f  %7.6f  %s" current_maximum_error start_mass*m_scale end_mass*m_scale -m_scale*p.Δm0[n][end] p.trust_region_factor termination_status(model)

            if current_maximum_error < p.dynamical_error
                printstyled(temp; color=:green)
            else
                printstyled(temp; color=:red)
            end
        end

        if i >= 10
            p.trust_region_factor = initial_trust_region_factor * ((scp_iterations - i) / (scp_iterations - 10))^2.0 + 5e-5
        end 
    end
end