abstract type SequentialConvexGuess end

struct BallisticGuess <: SequentialConvexGuess end 
struct LambertGuess <: SequentialConvexGuess end 


abstract type SequentialConvexDiscretization end

struct FixedTimeDiscretization <: SequentialConvexDiscretization end
struct FixedNumberDiscretization <: SequentialConvexDiscretization end



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
    mixing_number::Int64
    segment_number::Vector{Int64}
    objective_config::Union{MassConfig, LoggedMassConfig}
    dynamical_error::T
    mass_overhead::T
    trust_region_factor::T
    optimizer
end

Base.show(io::IO, p::SequentialConvexProblem) = begin
    println(io, "A sequential convex problem for GTOC12 using ", p.objective_config)
    println(io)

    total_deployments = 0
    total_collections = 0

    for n in 1:p.mixing_number
        deployments = sum(p.Δm0[n] .≈ -40/m_scale)
        collections = sum(p.Δm0[n] .> 0.0)

        total_deployments += deployments
        total_collections += collections

        println(io, "Ship $n ($deployments deployments, $collections collections)")
        println(io, "Sequence: ", join(get_id_label.(p.id_journey[n]), " → "))
        println(io, "Times: ", join([@sprintf("%.2f", convert_time_to_mjd(val)) for val in p.times_journey[n]], " → "))

        if typeof(p.objective_config) == LoggedMassConfig
            println(io, @sprintf("Fuel starting: %.2fkg", m_scale*exp(p.x_nodes[n][1][7, 1]) - 500.0 - 40.0*deployments))
            println(io, @sprintf("Fuel remaining: %.2fkg", m_scale*exp(p.x_nodes[n][end][7, end]) - 500.0 + m_scale*p.Δm0[n][end]))
        else
            println(io, @sprintf("Fuel starting: %.2fkg", m_scale*p.x_nodes[n][1][7, 1] - 500.0 - 40.0*deployments))
            println(io, @sprintf("Fuel remaining: %.2fkg", m_scale*p.x_nodes[n][end][7, end] - 500.0 + m_scale*p.Δm0[n][end]))
        end

        println(io, @sprintf("Mass returned: %.2fkg", -m_scale*p.Δm0[n][end]))
        println(io)
    end

    println(io, "Total ($total_deployments deployments, $total_collections collections)")
    println(io, @sprintf("Mass Returned: %.2fkg", -m_scale*sum([val[end] for val in p.Δm0])))

    println(io, 
        @sprintf("Average: %.2fkg (%.2f .. %.2f)", 
            -m_scale*mean([val[end] for val in p.Δm0]),
            -m_scale*maximum([val[end] for val in p.Δm0]),
            -m_scale*minimum([val[end] for val in p.Δm0]),
        )
    )
end


function SequentialConvexProblem(
    id_journey,
    times_journey;
    objective_config = LoggedMassConfig(),
    dynamical_error = 1e-6,
    trust_region_factor = 0.1,
    mass_overhead = 0.01/m_scale,
    optimizer = Mosek.Optimizer
)
    mixing_number = length(id_journey)
    segment_number = [length(id_journey[n]) - 1 for n in 1:mixing_number]

    t_nodes = Vector{Vector{Float64}}[]
    u_nodes = Vector{Matrix{Float64}}[]
    x_nodes = Vector{Matrix{Float64}}[]
    x0 = Vector{Vector{Float64}}[]
    xf = Vector{Vector{Float64}}[]
    Δv0 = Vector{Vector{Float64}}[]
    Δvf = Vector{Vector{Float64}}[]
    Δv0_limit = Vector{Float64}[]
    Δvf_limit = Vector{Float64}[]
    Δm0 = Vector{Float64}[]

    p = SequentialConvexProblem(
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
        mixing_number,
        segment_number,
        objective_config,
        dynamical_error,
        mass_overhead,
        trust_region_factor,
        optimizer
    )

    for n in 1:p.mixing_number
        push!(p.Δm0, get_mass_change_at_ids_mixed(
            p.id_journey[n],
            p.times_journey[n],
            p.id_journey,
            p.times_journey
        ))

        push!(p.x0, [ephermeris_cartesian_from_id(
            p.id_journey[n][i], 
            p.times_journey[n][i]
        )[:] for i in 1:(p.segment_number[n])])

        push!(p.xf, [ephermeris_cartesian_from_id(
            p.id_journey[n][i], 
            p.times_journey[n][i]
        )[:] for i in 2:(p.segment_number[n] + 1)])

        push!(p.t_nodes, fill(Float64[], p.segment_number[n]))
        push!(p.x_nodes, fill(Matrix{Float64}(undef, 0, 0), p.segment_number[n]))
        push!(p.u_nodes, fill(Matrix{Float64}(undef, 0, 0), p.segment_number[n]))

        push!(p.Δv0, fill(Float64[], p.segment_number[n]))
        push!(p.Δvf, fill(Float64[], p.segment_number[n]))

        push!(p.Δv0_limit, fill(0.0, p.segment_number[n]))
        push!(p.Δvf_limit, fill(0.0, p.segment_number[n]))
    end

    initialize_scp_discretization!(p, FixedTimeDiscretization())

    # initialize_scp_guess!(p, BallisticGuess())
    initialize_scp_guess!(p, LambertGuess())

    return p
end

function initialize_scp_discretization!(
    p::SequentialConvexProblem,
    ::FixedTimeDiscretization
)
    for n in 1:p.mixing_number, i in 1:p.segment_number[n]
        Δt_segment = p.times_journey[n][i+1] - p.times_journey[n][i]

        t_nodes_segment = collect(0.0:node_time_spacing:(Δt_segment - Δt_segment % node_time_spacing))

        if !(t_nodes_segment[end] ≈ Δt_segment)
            t_nodes_segment = vcat(t_nodes_segment, Δt_segment)
        end
    
        if p.id_journey[n][i] == p.id_journey[n][i + 1]
            t_nodes_segment = [0.0, Δt_segment/2, Δt_segment]
        end

        p.t_nodes[n][i] = t_nodes_segment
    end

    return
end



function initialize_scp_guess!(
    p::SequentialConvexProblem,
    ::BallisticGuess
)
    for n in 1:p.mixing_number
        m0 = 3000.0 / m_scale

        for k in 1:p.segment_number[n]
            p.Δv0[n][k] = [0.0, 0.0, 0.0]
            p.Δv0_limit[n][k] = 0.0

            p.Δvf[n][k] = [0.0, 0.0, 0.0]
            p.Δvf_limit[n][k] = 0.0

            mf = m0

            if typeof(p.objective_config) == LoggedMassConfig
                p.x0[n][k] = vcat(p.x0[n][k], log(m0))
                p.xf[n][k] = vcat(p.xf[n][k], log(mf))
            else
                p.x0[n][k] = vcat(p.x0[n][k], m0)
                p.xf[n][k] = vcat(p.xf[n][k], mf)
            end

            p.u_nodes[n][k] = zeros(4, length(p.t_nodes[n][k]) - 1)

            temp = deepcopy(p.x0[n][k])
            temp[4:6] .+= p.Δv0[n][k]

            p.x_nodes[n][k] = integrate_trajectory(
                temp,
                p.t_nodes[n][k];
                t_nodes = p.t_nodes[n][k],
                u_nodes = p.u_nodes[n][k],
                objective_config = p.objective_config,
            )

            m0 = mf + p.Δm0[n][k + 1]
        end
    end
end

function initialize_scp_guess!(
    p::SequentialConvexProblem,
    ::LambertGuess
)
    for n in 1:p.mixing_number
        m0 = 3000.0 / m_scale

        for k in 1:p.segment_number[n]
            Δv0, Δvf = if p.id_journey[n][k] == p.id_journey[n][k + 1]
                zeros(Float64, 3), zeros(Float64, 3)
            else
                temp = deepcopy(p.xf[n][k])

                tof = p.times_journey[n][k + 1] - p.times_journey[n][k]

                # Prevent near 180 degree transfers
                for _ in 1:10
                    if norm(cross(p.x0[n][k][1:3], temp[1:3])) > 0.5
                        break
                    else
                        tof += 1*day_scale

                        temp = ephermeris_cartesian_from_id(
                            p.id_journey[n][k + 1], 
                            p.times_journey[n][k] + tof
                        )
                    end
                end

                find_best_lambert_transfer(p.x0[n][k], temp, tof)
            end

            p.Δv0[n][k] = Δv0[:]
            p.Δv0_limit[n][k] = (p.id_journey[n][k] == 0) ? 6.0/v_scale : 0.0

            p.Δvf[n][k] = Δvf[:]
            p.Δvf_limit[n][k] = ((k == p.segment_number[n]) && p.id_journey[n][k+1] == -3) ? 6.0/v_scale : 0.0

            Δv_use = max(norm(p.Δv0[n][k]) - p.Δv0_limit[n][k], 0.0) + max(norm(p.Δvf[n][k]) - p.Δvf_limit[n][k], 0.0) 

            mf = m0/exp(Δv_use / g0_isp)

            if typeof(p.objective_config) == LoggedMassConfig
                p.x0[n][k] = vcat(p.x0[n][k], log(m0))
                p.xf[n][k] = vcat(p.xf[n][k], log(mf))
            else
                p.x0[n][k] = vcat(p.x0[n][k], m0)
                p.xf[n][k] = vcat(p.xf[n][k], mf)
            end

            temp = deepcopy(p.x0[n][k])
            temp[4:6] .+= Δv0

            p.u_nodes[n][k] = zeros(4, length(p.t_nodes[n][k]) - 1)

            # Use Lambert transfer as guess for x_nodes
            p.x_nodes[n][k] = integrate_trajectory(
                temp,
                p.t_nodes[n][k];
                objective_config = p.objective_config,
            )

            m0 = mf + p.Δm0[n][k + 1]
        end
    end
end


function convert_logged_mass_to_mass!(
    p::SequentialConvexProblem
)
    if typeof(p.objective_config) != LoggedMassConfig
        error("Conversion to mass problem attempted on non logged-mass problem.")
        return
    end

    p.objective_config = MassConfig()

    for n in 1:p.mixing_number, k in 1:p.segment_number[n]
        mass = exp.(p.x_nodes[n][k][7, :])

        p.x0[n][k][7] = exp(p.x0[n][k][7])
        p.xf[n][k][7] = exp(p.xf[n][k][7])

        p.u_nodes[n][k][1, :] .*= mass[1:end-1]
        p.u_nodes[n][k][2, :] .*= mass[1:end-1]
        p.u_nodes[n][k][3, :] .*= mass[1:end-1]
        p.u_nodes[n][k][4, :] .*= mass[1:end-1]

        p.x_nodes[n][k] = integrate_trajectory(
            p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
            p.t_nodes[n][k];
            t_nodes = p.t_nodes[n][k],
            u_nodes = p.u_nodes[n][k],
            p.objective_config,
        )
    end

    return
end




function solve!(
    p::SequentialConvexProblem;
    fixed_segments = false,
    fixed_rendezvous = false,
    maximum_mass = false
)
    if fixed_segments && !fixed_rendezvous
        error("Must have non-fixed segments for adaptive rendezvous")
    end

    mixing = length(p.t_nodes)

    id_groups = get_journey_groups(p.id_journey)

    models = []

    for _ in 1:maximum(id_groups)
        model = Model(p.optimizer)
        set_silent(model)

        push!(models, model)
    end

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
    m_violation = [@variable(models[id_groups[n]], lower_bound = 0.0, upper_bound = 1000.0) for n in 1:mixing]
    Δt_start = [@variable(models[id_groups[n]]) for n in 1:mixing]

    for n in 1:mixing
        segments = length(p.x0[n])

        # For each segment in each trajectory
        for k in 1:segments
            # Number of nodes in current segment
            nodes = length(p.t_nodes[n][k])

            # Variables
            push!(x[n], @variable(models[id_groups[n]], [i=1:7, j=1:nodes], start = p.x_nodes[n][k][i, j]))
            push!(u[n], @variable(models[id_groups[n]], [i=1:4, j=1:(nodes-1)], start = p.u_nodes[n][k][i, j]))
            push!(s[n], @variable(models[id_groups[n]], [i=1:(nodes-1)], start = 1.0))
            push!(Δx0[n], @variable(models[id_groups[n]], [i=1:6]))
            push!(Δxf[n], @variable(models[id_groups[n]], [i=1:6]))
            push!(x_violation[n], @variable(models[id_groups[n]], [i=1:6]))
            push!(Δv0[n], @variable(models[id_groups[n]], [i=1:3], start = p.Δv0[n][k][i]))
            push!(Δvf[n], @variable(models[id_groups[n]], [i=1:3], start = p.Δvf[n][k][i]))

            # SOC constraint for control
            @constraint(models[id_groups[n]], [i=1:(nodes-1)], [u[n][k][4, i]; u[n][k][1:3, i]] in SecondOrderCone())

            # SOC constraint for injection Δv
            @constraint(models[id_groups[n]], [p.Δv0_limit[n][k]; Δv0[n][k][1:3]] in SecondOrderCone())
            @constraint(models[id_groups[n]], [p.Δvf_limit[n][k]; Δvf[n][k][1:3]] in SecondOrderCone())

            # Ensure positivity of the violation
            @constraint(models[id_groups[n]], x_violation[n][k] .>= Δx0[n][k])
            @constraint(models[id_groups[n]], x_violation[n][k] .>= -Δx0[n][k])
            @constraint(models[id_groups[n]], x_violation[n][k] .>= Δxf[n][k])
            @constraint(models[id_groups[n]], x_violation[n][k] .>= -Δxf[n][k])

            if fixed_segments
                @constraint(models[id_groups[n]], s[n][k] .== 1.0)
            end
        end

        if fixed_segments || fixed_rendezvous
            @constraint(models[id_groups[n]], Δt_start[n] .== 0.0)
        end

        # Limit maximum mass at start
        mass_limit = if !maximum_mass 
            3000.0/m_scale
        else
            25000.0/m_scale
        end

        if typeof(p.objective_config) == LoggedMassConfig
            # @constraint(model, x[n][1][7, 1] <= log(1000/m_scale))
            @constraint(models[id_groups[n]], x[n][1][7, 1] <= log(mass_limit))
        else
            # @constraint(model, x[n][1][7, 1] <= 1000/m_scale)
            @constraint(models[id_groups[n]], x[n][1][7, 1] <= mass_limit)
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

    active_models = collect(1:length(models))

    current_trust_region_factor = p.trust_region_factor

    for i in 1:scp_iterations
        Δt_nodes = [[t_nodes[2:end] .- t_nodes[1:(end-1)] for t_nodes in t_nodes] for t_nodes in p.t_nodes]

        Δt_segments = [@expression(models[id_groups[n]], [i=1:length(p.x0[n])],
            sum(s[n][i] .* Δt_nodes[n][i])
        ) for n in 1:mixing]

        actual_time = [@expression(models[id_groups[n]], [i=1:length(p.times_journey[n])], 
            p.times_journey[n][1] + sum(vcat([0.0], Δt_segments[n][1:(i-1)])) + Δt_start[n]
        ) for n in 1:mixing]

        objective = fill(AffExpr(0.0), length(models))

        for n in 1:mixing
            if id_groups[n] ∉ active_models
                continue
            end

            for k in bad_dynamic_segments[n]
                nodes = size(p.t_nodes[n][k], 1)

                for con in dynamic_con[n][k]
                    delete(models[id_groups[n]], con)
                end
                
                for con in delation_dynamic_con[n][k]
                    delete(models[id_groups[n]], con)
                end

                for con in trust_region_dynamic_con[n][k]
                    delete(models[id_groups[n]], con)
                end

                for con in segment_time_con[n][k]
                    delete(models[id_groups[n]], con)
                end

                for con in segment_start_con[n][k]
                    delete(models[id_groups[n]], con)
                end

                for con in segment_end_con[n][k]
                    delete(models[id_groups[n]], con)
                end

                # Get the state transition matrices from automatic differentiation
                stms = get_arc_stms_adaptive_time(
                    Δt_nodes[n][k], 
                    p.x_nodes[n][k], 
                    p.u_nodes[n][k]; 
                    coordinate_config = CartesianConfig(), 
                    thrust_config = ZeroOrderHoldConfig(), 
                    objective_config = p.objective_config
                )

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

                dynamic_con[n][k] = @constraint(models[id_groups[n]],
                    [i=1:(nodes-1)],
                    x[n][k][:, i+1] .== As[:, :, i] * x[n][k][:, i] .+ Bs[:, :, i] * u[n][k][:, i] .+ Ts[:, i] * (Δt_nodes[n][k][i] .* s[n][k][i]) .+ cs[:, i]
                )

                delation_dynamic_con[n][k] = @constraint(models[id_groups[n]],
                    1.0 - current_trust_region_factor .<= s[n][k] .<= 1.0 + current_trust_region_factor
                )

                trust_region_dynamic_con[n][k] = @constraint(models[id_groups[n]],
                    [i=1:(nodes-1)],
                    # -1e-0 .<= x[n][k][7, i] .- p.x_nodes[n][k][7, i] .<= 1e-0
                    -5e1*current_trust_region_factor - 1e-3 .<= x[n][k][1:6, i] .- p.x_nodes[n][k][1:6, i] .<= 5e1*current_trust_region_factor + 1e-3
                )

                start_position_function(time) = begin
                    ephermeris_cartesian_from_id(p.id_journey[n][k], time)[:]
                end 

                final_position_function(time) = begin
                    ephermeris_cartesian_from_id(p.id_journey[n][k+1], time)[:]
                end 

                start_jac = FiniteDiff.finite_difference_jacobian(start_position_function, [p.times_journey[n][k]])
                final_jac = FiniteDiff.finite_difference_jacobian(final_position_function, [p.times_journey[n][k + 1]])

                segment_start_con[n][k] = @constraint(models[id_groups[n]], 
                    x[n][k][1:6, 1] .== 
                        (p.x0[n][k][1:6] .+ vcat([0, 0, 0], Δv0[n][k]) .+ 
                        start_jac * (actual_time[n][k] - p.times_journey[n][k]) .+ Δx0[n][k])[:]
                )      

                segment_end_con[n][k] = @constraint(models[id_groups[n]], 
                    x[n][k][1:6, nodes] .== 
                        (p.xf[n][k][1:6] .- vcat([0, 0, 0], Δvf[n][k][1:3]) .+ 
                        final_jac * (actual_time[n][k+1] - p.times_journey[n][k + 1]) .+ Δxf[n][k])[:]
                )    
            end
            
            # Control constraints due to mass
            for k in bad_thrust_mass_segments[n]
                nodes = size(p.t_nodes[n][k], 1)

                for con in thrust_mass_con[n][k]
                    delete(models[id_groups[n]], con)
                end

                thrust_mass_con[n][k] = if typeof(p.objective_config) == LoggedMassConfig
                    @constraint(models[id_groups[n]],
                        [i=1:(nodes-1)],
                        u[n][k][4, i] <= exp.(-p.x_nodes[n][k][7, i]).*(1.0 .- x[n][k][7, i] .+ p.x_nodes[n][k][7, i])
                    )
                else
                    @constraint(models[id_groups[n]],
                        [i=1:(nodes-1)],
                        u[n][k][4, i] <= 1.0
                    )
                end
            end

            # Mass linkage constraints
            for k in bad_mass_link_segments[n]
                for con in mass_link_con[n][k]
                    delete(models[id_groups[n]], con)
                end

                time_indices = []

                for l in collect(1:p.mixing_number)[id_groups .== id_groups[n]]
                    for m in 1:length(p.id_journey[l])
                        if p.id_journey[n][k] == p.id_journey[l][m]
                            push!(time_indices, (l, m))
                        end
                    end
                end

                sort!(time_indices; by = x->p.times_journey[x[1]][x[2]])

                Δm = if n == time_indices[1][1] && k == time_indices[1][2]
                    -40/m_scale
                else
                    mining_rate*(actual_time[time_indices[2][1]][time_indices[2][2]] - actual_time[time_indices[1][1]][time_indices[1][2]])
                end

                mass_link_con[n][k] = if typeof(p.objective_config) == LoggedMassConfig 
                    [@constraint(models[id_groups[n]],
                        x[n][k][7, 1] == log(exp(p.x_nodes[n][k-1][7, end]) + p.Δm0[n][k]) + (1.0/(1.0 + p.Δm0[n][k]/exp(p.x_nodes[n][k-1][7, end])))*(x[n][k-1][7, end] - p.x_nodes[n][k-1][7, end])
                    )]
                else
                    [@constraint(models[id_groups[n]],
                        # x[n][k][7, 1] == x[n][k-1][7, end] + Δm
                        x[n][k][7, 1] == x[n][k-1][7, end] + p.Δm0[n][k]
                    )]
                end
            end

            for con in mass_end_con[n]
                delete(models[id_groups[n]], con)
            end

            if !maximum_mass
                mass_end_con[n] = if typeof(p.objective_config) == LoggedMassConfig
                    [@constraint(models[id_groups[n]], x[n][end][7, end] + m_violation[n] >= log(500/m_scale + p.mass_overhead - p.Δm0[n][end]))]
                else
                    [@constraint(models[id_groups[n]], x[n][end][7, end] + m_violation[n] >= 500/m_scale + p.mass_overhead - p.Δm0[n][end])]
                end
            end

            for con in actual_time_con[n]
                delete(models[id_groups[n]], con)
            end

            if !fixed_rendezvous
                for con in time_start_movement_con[n]
                    delete(models[id_groups[n]], con)
                end

                time_start_movement_con[n] = [@constraint(models[id_groups[n]], -2e-1*current_trust_region_factor <= Δt_start[n] <= 2e-1*current_trust_region_factor)]

                actual_time_con[n] = [
                    @constraint(models[id_groups[n]], actual_time[n][end] <= maximum_time - 10*current_trust_region_factor*day_scale),
                    @constraint(models[id_groups[n]], actual_time[n][1] >= 0.0 + 10*current_trust_region_factor*day_scale),
                ]
            else
                actual_time_con[n] = [
                    @constraint(models[id_groups[n]], actual_time[n][i] == p.times_journey[n][i]) for i in 1:length(p.times_journey[n])
                ]
            end

            dropoff = p.Δm0[n] .≈ -40/m_scale
            pickup = p.Δm0[n] .> 0.0

            if fixed_rendezvous
                if maximum_mass 
                    objective[id_groups[n]] += x[n][end][7, end]
                else
                    objective[id_groups[n]] -= x[n][1][7, 1]
                end
            else
                objective[id_groups[n]] += sum(actual_time[n][pickup]) - sum(actual_time[n][dropoff])
            end

            objective[id_groups[n]] -= 5e3*m_violation[n]
            objective[id_groups[n]] -= 1e4*sum(sum.(x_violation[n]))
            # objective -= 1e-1*sum([sum(u[n][i][4, :].*Δt_nodes[n][i]) for i in 1:length(u[n])])
        end

        for i in active_models
            @objective(models[i], 
                Max, 
                objective[i]
            )

            # display(models[i])

            # display(sum(num_constraints(models[i], F, S) for (F, S) in list_of_constraint_types(models[i])))

            JuMP.optimize!(models[i])
        end
        
        for n in 1:mixing
            if id_groups[n] ∉ active_models
                continue
            end

            p.times_journey[n] = JuMP.value.(actual_time[n])
            p.u_nodes[n] = [JuMP.value.(u) for u in u[n]]
            p.Δv0[n] = [JuMP.value.(Δv0) for Δv0 in Δv0[n]]
            p.Δvf[n] = [JuMP.value.(Δvf) for Δvf in Δvf[n]]

            p.t_nodes[n] = [vcat(0.0, cumsum(JuMP.value.(s[n][k]) .* Δt_nodes[n][k])) for k in 1:length(p.x0[n])]
        end

        maximum_error = fill(0.0, mixing)

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

                p.xf[n][k][1:6] = ephermeris_cartesian_from_id(
                    p.id_journey[n][k+1], 
                    p.times_journey[n][k+1]
                )[:]

                if k != 1
                    if typeof(p.objective_config) == LoggedMassConfig
                        p.x0[n][k][7] = log(max(1/6, exp(p.xf[n][k - 1][7]) + p.Δm0[n][k]))
                    else
                        p.x0[n][k][7] = max(1/6, p.xf[n][k - 1][7] + p.Δm0[n][k])
                    end
                else
                    p.x0[n][k][7] = value.(x[n][k][7, 1])
                end

                p.u_nodes[n][k][4, :] = norm.(eachcol(p.u_nodes[n][k][1:3, :]))

                p.x_nodes[n][k] = integrate_trajectory(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    p.t_nodes[n][k];
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                push!(dynamical_errors, abs(p.xf[n][k][7] - p.x_nodes[n][k][7, end]))

                p.xf[n][k][7] = p.x_nodes[n][k][7, end]

                push!(
                    dynamical_errors, 
                    maximum(abs.(p.xf[n][k][1:6] .- p.x_nodes[n][k][1:6, end] .- vcat([0.0, 0.0, 0.0], p.Δvf[n][k])))
                )

                if k > 1
                    mass_link_error = if typeof(p.objective_config) == LoggedMassConfig
                        abs(exp(JuMP.value.(x[n][k-1][7, end])) + p.Δm0[n][k] - exp(JuMP.value.(x[n][k][7, 1])))
                    else
                        abs(JuMP.value.(x[n][k-1][7, end]) + p.Δm0[n][k] - JuMP.value.(x[n][k][7, 1]))
                    end

                    # push!(
                    #     mass_link_errors,
                    #     mass_link_error
                    # )
                end
            end

            maximum_error[n] = maximum(vcat(dynamical_errors, mass_link_errors))

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

            mixed_ship = if count(==(id_groups[n]), id_groups) == 1
                " "
            else
                "*"
            end
            
            temp = @sprintf "%2i%1s  %10.6e  %10.5f  %10.5f  %10.5f  %7.6f  %s" id_groups[n] mixed_ship maximum_error[n] start_mass*m_scale end_mass*m_scale -m_scale*p.Δm0[n][end] current_trust_region_factor termination_status(models[id_groups[n]])

            if maximum_error[n] < p.dynamical_error
                printstyled(temp; color=:green)
            else
                printstyled(temp; color=:red)
            end
        end

        for n in 1:length(models)
            if (maximum(maximum_error[id_groups .== n]) < p.dynamical_error) || (termination_status(models[n]) ∉ [OPTIMAL, SLOW_PROGRESS])
                setdiff!(active_models, n)
            end
        end

        if length(active_models) == 0
            break
        end

        if i >= 10
            current_trust_region_factor = max(p.trust_region_factor * ((scp_iterations - i) / (scp_iterations - 10))^2.0, 1e-4)
        end 
    end
end