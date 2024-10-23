


function get_id_label(
    id
)
    return if id == 0
        "Earth"
    elseif id < 0
        ["Venus", "Earth", "Mars"][-id - 1]
    else
        string(id)
    end
end

function get_event_label_and_color(
    id,
    mass_change
)
    label = get_id_label(id)

    color = :black

    return label, color
end



function plot_location_vertical_markers(
    ax,
    times_journey,
    id_journey,
    mass_changes;
    vertical_position = 0.0
)
    for i in 1:length(times_journey)
        label, color = get_event_label_and_color(id_journey[i], mass_changes[i])

        lines!(ax,
            [times_journey[i], times_journey[i]],
            [-10.0, vertical_position],
            color=color,
            linewidth=1.0,
            alpha=0.5,
            label=false,
        )

        text!(ax,
            [(times_journey[i], vertical_position)],
            text = label,
            align = (:left, :baseline),
            rotation = π/4,
            fontsize = 10
        )
    end
end


function plot_thrust_information(
    p::SequentialConvexProblem;
    solution_indices = nothing
)
    mixing = length(p.x0)

    if isnothing(solution_indices)
        solution_indices = collect(1:mixing)
    end
    
    f = Figure(size = (800, 250*length(solution_indices)), backgroundcolor = :white)

    axs = []

    for (i, n) in enumerate(solution_indices)
        t_nodes_combined = vcat(
            [p.t_nodes[n][k][1:end-1] .+ p.times_journey[n][k] for k in 1:length(p.x0[n])]...
        )
        u_nodes_combined = hcat(p.u_nodes[n]...)
        x_nodes_combined = hcat([val[:, 1:end-1] for val in p.x_nodes[n]]...)

        push!(t_nodes_combined, p.times_journey[end][end])
        u_nodes_combined = hcat(u_nodes_combined, p.u_nodes[end][end][:, end])
        x_nodes_combined = hcat(x_nodes_combined, p.x_nodes[end][end][:, end])

        thrust_force = if p.objective_config == LoggedMassConfig()
            temp = exp.(x_nodes_combined[7, :])
            u_nodes_combined[4, :] .* temp * thrust * m_scale * a_scale * 1e3
        else
            u_nodes_combined[4, :] * thrust * m_scale * a_scale * 1e3
        end
        
        ax = Axis(
            f[i, 1]; 
            xlabel = n == solution_indices[end] ? "t [MJD]" : "", 
            ylabel = "thrust [N]", 
            # xticks = [65000, 66000, 67000, 68000, 69000],
            yticks = [0.0, 0.6],
            xgridvisible = false
        )

        push!(axs, ax)

        ylims!(ax, (-0.05 * maximum(thrust_force), 1.25 * maximum(thrust_force)))

        plot_location_vertical_markers(
            ax, 
            convert_time_to_mjd(p.times_journey[n]), 
            p.id_journey[n], 
            p.Δm0[n]; 
            vertical_position = 1.05*maximum(thrust_force)
        )

        stairs!(
            ax, 
            convert_time_to_mjd(t_nodes_combined), 
            thrust_force;
            step = :post,
            color = :black,
            linewidth = 1.5,
        )
    end

    linkxaxes!(axs...)
    resize_to_layout!(f)
    display(f)

    return
end


function plot_trajectory(
    p::SequentialConvexProblem;
    plot_3d = false,
    solution_indices = nothing,
    rotating = false
)
    f = Figure(size = (900, 800), backgroundcolor = :white)

    if isnothing(solution_indices)
        solution_indices = collect(1:p.mixing_number)
    end

    ax1 = if plot_3d
        Axis3(
            f[1:2, 1:2]; 
            xlabel = "x [AU]", 
            ylabel = "y [AU]", 
            zlabel = "z [AU]", 
            # limits = (-3.0, 3.0, -3.0, 3.0, -3.0, 3.0),
            # aspect = 1,
            # perspectiveness = 1.0,
        )
    else
        Axis(
            f[1:2, 1:2]; 
            xlabel = "x [AU]", 
            ylabel = "y [AU]", 
            limits = (-3.1, 3.1, -3.1, 3.1),
            # limits = (0.5, 1.5, -3.0, -2.0),
            aspect = 1
        )
    end

    rotation_rate = if rotating
        visited_ids = setdiff(unique(stack(p.id_journey[solution_indices])), [0, -3])

        a_visited = asteroids_classical[1, visited_ids]

        period_visited = 2π * sqrt.(a_visited.^3/μ)

        mean(period_visited)
    else
        0.0
    end

    # scatter!(ax1,
    #     [0.0],
    #     [0.0],
    #     [0.0],
    #     color = :black
    # )

    ν_plot = collect(LinRange(0.0, 2π, 400)) 

    for (i, planet_classical) in enumerate(eachcol(planets_classical))
        temp = repeat(planet_classical, 1, 400)
        temp[6, :] .= ν_plot
        temp = hcat(classical_to_cartesian.(eachcol(temp))...)

        lines!(ax1,
            temp[1, :],
            temp[2, :],
            temp[3, :],
            linestyle = :dashdot,
            color = Makie.wong_colors()[[2, 1, 6][i]]
        )
    end

    for n in solution_indices
        for k in 1:p.segment_number[n]
            x0 = p.x0[n][k]
            xf = p.xf[n][k]

            if rotating
                x0 = get_state_rotation_matrix(2π*(p.times_journey[n][k])/rotation_rate)*x0[1:6]
                xf = get_state_rotation_matrix(2π*(p.times_journey[n][k+1])/rotation_rate)*xf[1:6]
            end

            scatter!(ax1,
                x0[1],
                x0[2],
                x0[3],
            )

            if k == p.segment_number[n]
                scatter!(ax1,
                    xf[1],
                    xf[2],
                    xf[3],
                )
            end

            t_fine = collect(p.t_nodes[n][k][1]:1.0*day_scale:p.t_nodes[n][k][end])

            if !(t_fine[end] ≈ p.t_nodes[n][k][end])
                push!(t_fine, p.t_nodes[n][k][end])
            end

            x_fine = integrate_trajectory(
                p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                t_fine;
                t_nodes = p.t_nodes[n][k],
                u_nodes = p.u_nodes[n][k],
                p.objective_config,
            )

            if rotating
                rotation_matrix = get_state_rotation_matrix.(2π*(t_fine .+ p.times_journey[n][k])/rotation_rate)

                x_fine[1:6, :] = stack(rotation_matrix.*eachcol(x_fine[1:6, :]))
            end

            lines!(ax1,
                x_fine[1, :],
                x_fine[2, :],
                x_fine[3, :],
            )

            thrust_force = if p.objective_config == LoggedMassConfig()
                temp = exp.(p.x_nodes[n][k][7, 1:end-1])
                p.u_nodes[n][k][4, :] .* temp * thrust * m_scale * a_scale * 1e3
            else
                p.u_nodes[n][k][4, :] * thrust * m_scale * a_scale * 1e3
            end

            sel = thrust_force .>= 0.01

            thrust_vectors = stack(normalize.(eachcol(p.u_nodes[n][k][1:3, :]))).*thrust_force'

            length_scale = 0.2

            for i in collect(1:length(thrust_force))[sel]
                lines!(ax1,
                    [p.x_nodes[n][k][1, i], p.x_nodes[n][k][1, i] + length_scale*thrust_vectors[1, i]],
                    [p.x_nodes[n][k][2, i], p.x_nodes[n][k][2, i] + length_scale*thrust_vectors[2, i]],
                    [p.x_nodes[n][k][3, i], p.x_nodes[n][k][3, i] + length_scale*thrust_vectors[3, i]],
                    color = :black,
                    alpha = 0.5,
                    linewidth = 1,
                )
            end
        end
    end

    
    # hidedecorations!(
    #     ax1
    # )

    resize_to_layout!(f)
    display(f)

    return
end


function plot_graph(
    p::MixedIntegerProblem;
    plot_optimal_path = false,
    plot_pruned = true,
    output_file = nothing,
    selection_index = 1,
    figure_size = (900, 340)
)

    stage_size = length(p.id_subset)

    stage_Δx = 0.1
    node_Δy = 0.1

    unused_edge_list = []

    used_edge_lists = [[] for _ in 1:p.solutions]
    used_edge_colors = [[] for _ in 1:p.solutions]
    used_edge_styles = [[] for _ in 1:p.solutions]
    used_edge_widths = [[] for _ in 1:p.solutions]

    layout = []

    # color_range = ColorSchemes.tab10
    color_range = ColorSchemes.tableau_10
    # color_range = ColorSchemes.tableau_jewel_bright
    # color_range = ColorSchemes.tableau_miller_stone

    for stage in 1:p.deployment_nums[selection_index]
        for i in 1:stage_size
            push!(layout, ((stage - 1)*stage_Δx, -(i-1)*node_Δy))
        end
    end

    for stage in 1:p.deployment_arcs[selection_index]
        for i in 1:stage_size, j in 1:stage_size
            if !plot_pruned || ((p.deployment_cost[selection_index][i, j, stage] <= p.cost_limit) && (i != j))
                used_check = false

                temp = []

                if plot_optimal_path
                    for k in p.solutions:-1:1
                    # for k in 1:p.solutions
                        if p.id_journey_solutions[k][selection_index][stage + 1] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][stage + 2] == p.id_subset[j]
                            used_check = true

                            for l in temp
                                used_edge_widths[l][end] += 3
                            end

                            push!(used_edge_lists[k], (i + stage_size*(stage-1), j + stage_size*stage))
                            push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                            push!(used_edge_styles[k], :solid)
                            push!(used_edge_widths[k], 4)
                            push!(temp, k)
                        end
                    end
                end

                if !used_check
                    push!(unused_edge_list, (i + stage_size*(stage-1), j + stage_size*stage))
                end
            end
        end
    end

    collection_offset = length(layout)

    for i in 1:stage_size, j in 1:stage_size
        if !plot_pruned || (p.intermediate_cost[selection_index][i, j] <= p.cost_limit)
            used_check = false

            temp = []

            if plot_optimal_path
                for k in p.solutions:-1:1
                # for k in 1:p.solutions
                    if p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 1] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 2] == p.id_subset[j]
                        used_check = true

                        for l in temp
                            used_edge_widths[l][end] += 3
                        end
                    
                        push!(used_edge_lists[k], (collection_offset - stage_size + i, collection_offset + j))
                        push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                        push!(used_edge_styles[k], :solid)
                        push!(used_edge_widths[k], 4)
                        push!(temp, k)
                    end
                end
            end

            if !used_check
                push!(unused_edge_list, (collection_offset - stage_size + i, collection_offset + j))
            end
        end
    end

    for stage in (p.deployment_nums[selection_index] + 1):(p.deployment_nums[selection_index] + p.collection_nums[selection_index])
        for i in 1:stage_size
            push!(layout, ((stage - 1 + 1)*stage_Δx, -(i-1)*node_Δy))
        end
    end

    for stage in 1:p.collection_arcs[selection_index]
        for i in 1:stage_size, j in 1:stage_size
            if !plot_pruned || ((p.collection_cost[selection_index][i, j, stage] <= p.cost_limit) && (i != j))
                used_check = false

                temp = []

                if plot_optimal_path
                    for k in p.solutions:-1:1
                    # for k in 1:p.solutions
                        if p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 1 + stage] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 2 + stage] == p.id_subset[j]
                            used_check = true

                            for l in temp
                                used_edge_widths[l][end] += 3
                            end
                        
                            push!(used_edge_lists[k], (collection_offset + i + stage_size*(stage-1), collection_offset + j + stage_size*stage))
                            push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                            push!(used_edge_styles[k], :solid)
                            push!(used_edge_widths[k], 4)
                            push!(temp, k)
                        end
                    end
                end

                if !used_check
                    push!(unused_edge_list, (collection_offset + i + stage_size*(stage-1), collection_offset + j + stage_size*stage))
                end
            end
        end
    end
    
    unused_graph = SimpleDiGraph(length(layout))

    [add_edge!(unused_graph, e) for e in unused_edge_list]

    fixed_layout(_) = layout

    f = Figure(size = figure_size, backgroundcolor = :white)

    ax = Axis(
        f[1, 1];
        # size = (800, 800),
    )


    graphplot!(ax, unused_graph, 
        layout=fixed_layout,
        # size = (800, 800),
        edge_plottype=:beziersegments,
        edge_color= (:black, 0.1),
        edge_width=2,
        arrow_show = false,
        node_size = 0,
    )

    for k in p.solutions:-1:1
    # for k in 1:p.solutions
        used_graph = SimpleDiGraph(length(layout))

        [add_edge!(used_graph, e) for e in used_edge_lists[k]]

        graphplot!(ax, used_graph, 
            layout=fixed_layout,
            edge_plottype=:beziersegments,
            edge_color=used_edge_colors[k],
            edge_attr = (; linestyle = used_edge_styles[k]),
            edge_width= used_edge_widths[k],
            arrow_show = false,
            node_size = 0,
        )
    end


    empty_graph = SimpleDiGraph(length(layout))

    graphplot!(ax, empty_graph, 
        layout=fixed_layout,
        ilabels = repeat(p.id_subset, p.deployment_nums[selection_index] + p.collection_nums[selection_index]),
        ilabels_fontsize = 9,
        node_size = 30,
        node_marker=Circle;
    )

    time_positions = vcat(
        [stage_Δx*(stage-1) for stage in 1:p.deployment_nums[selection_index]],
        [stage_Δx*(stage-1+1) for stage in (p.deployment_nums[selection_index] + 1):(p.deployment_nums[selection_index] + p.collection_nums[selection_index])]
    )

    time_strings = [@sprintf("MJD %.0f", v) for v in convert_time_to_mjd.(mip_problem.times_journey[selection_index][2:end-1])]

    text!(ax, 
        time_positions,
        fill(0.7*node_Δy, length(time_positions)),
        text = time_strings,
        align = (:center, :center),
        # fontsize = 12,
    )

    hidedecorations!(ax)
    hidespines!(ax)

    ylims!(ax, [-node_Δy*(length(p.id_subset) - 0.5), node_Δy])
    

    # resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end




function plot_trajectory_paper(
    p::SequentialConvexProblem;
    plot_3d = false,
    solution_indices = nothing,
    rotating = false,
    output_file = nothing
)
    f = Figure(size = (800, 650), backgroundcolor = :white, figure_padding = 2)

    if isnothing(solution_indices)
        solution_indices = collect(1:p.mixing_number)
    end

    ax1 = if plot_3d
        Axis3(
            f[1:10, 1]; 
            xlabel = "x [AU]", 
            ylabel = "y [AU]", 
            zlabel = "z [AU]", 
            # limits = (2.2, 3.2, -0.5, 0.5, -0.5, 0.5),
            limits = (-3.0, 3.0, -3.0, 3.0, -3.0, 3.0),
            # aspect = 1,
            # perspectiveness = 1.0,
        )
    else
        Axis(
            f[1:10, 1]; 
            # xlabel = "x [AU]", 
            ylabel = "y [AU]", 
            # limits = (2.2, 3.2, -0.5, 0.5),
            # limits = (-0.25, 3.2, -1.1, 1.1),
            # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
            # yticks = [-1.0, -0.5, 0.0, 0.5, 1.0],
            # xticklabelsvisible = false,
            # limits = (0.5, 1.5, -3.0, -2.0),
            # aspect = 3.7/3
            # aspect = 1.568
        )
    end

    # ax2 = Axis(
    #     f[1:2, 3]; 
    #     xlabel = "x [AU]", 
    #     ylabel = "y [AU]", 
    #     # limits = (2.2, 3.2, -0.5, 0.5),
    #     limits = (2.55, 3.05, -0.45, 0.45),
    #     yticks = [-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4],
    #     # limits = (0.5, 1.5, -3.0, -2.0),
    #     aspect = 0.5/0.9
    # )

    ax3 = Axis(
        f[11:14, 1]; 
        xlabel = "x [AU]", 
        ylabel = "z [AU]", 
        # limits = (2.2, 3.2, -0.5, 0.5),
        # limits = (-0.25, 3.2, -0.35, 0.35),
        # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        # yticks = [-0.2, 0.0, 0.2],
        # limits = (0.5, 1.5, -3.0, -2.0),
        # aspect = 3.7/3
        # aspect = 1.32
    )

    # ax4 = Axis(
    #     f[3, 3]; 
    #     xlabel = "x [AU]", 
    #     ylabel = "z [AU]", 
    #     # limits = (2.2, 3.2, -0.5, 0.5),
    #     limits = (2.55, 3.05, -0.3, 0.3),
    #     # yticks = [-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4],
    #     # limits = (0.5, 1.5, -3.0, -2.0),
    #     # aspect = 0.5/0.9
    # )


    rotation_rate = if rotating
        visited_ids = setdiff(unique(vcat(p.id_journey[solution_indices]...)), [0, -3])

        a_visited = asteroids_classical[1, visited_ids]

        period_visited = 2π * sqrt.(a_visited.^3/μ)

        mean(period_visited)
    else
        0.0
    end

    visited_angle = -1.25

    # scatter!(ax1,
    #     [0.0],
    #     [0.0],
    #     [0.0],
    #     color = :black
    # )

    # ν_plot = collect(LinRange(0.0, 2π, 400)) 

    for (ax, order) in zip([ax1, ax3], [[1, 2, 3], [1, 3, 2]])

        t_planet = convert_mjd_to_time.(collect(mjd_start:10:mjd_end))

        for (i, planet_classical) in enumerate(eachcol(planets_classical))
            # temp = repeat(planet_classical, 1, 400)
            # temp[6, :] .= ν_plot
            # temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            # temp = 

            xp = integrate_trajectory(
                planets_cartesian[:, i],
                t_planet
            )

            if rotating
                xp = stack(get_state_rotation_matrix.(2π*t_planet/rotation_rate .+ visited_angle).*eachcol(xp[1:6, :]))
            end

            lines!(ax,
                xp[order[1], :],
                xp[order[2], :],
                xp[order[3], :],
                linestyle = :dashdot,
                color = Makie.wong_colors()[[2, 1, 6][i]],
                alpha = 0.35,
                label = ["Venus Orbit", "Earth Orbit", "Mars Orbit"][i]
            )
        end

        for n in solution_indices
            for k in 1:p.segment_number[n]
                x0 = p.x0[n][k]
                xf = p.xf[n][k]

                if rotating
                    x0 = get_state_rotation_matrix(2π*(p.times_journey[n][k])/rotation_rate + visited_angle)*x0[1:6]
                    xf = get_state_rotation_matrix(2π*(p.times_journey[n][k+1])/rotation_rate + visited_angle)*xf[1:6]
                end

                scatter!(ax,
                    x0[order[1]],
                    x0[order[2]],
                    x0[order[3]],
                )

                if k == p.segment_number[n]
                    scatter!(ax,
                        xf[order[1]],
                        xf[order[2]],
                        xf[order[3]],
                    )
                end

                t_fine = collect(p.t_nodes[n][k][1]:1.0*day_scale:p.t_nodes[n][k][end])

                if !(t_fine[end] ≈ p.t_nodes[n][k][end])
                    push!(t_fine, p.t_nodes[n][k][end])
                end

                x_fine = integrate_trajectory(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    t_fine;
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                x_target_fine = if p.id_journey[n][k + 1] != -3
                    integrate_trajectory(
                        ephermeris_cartesian_from_id(p.id_journey[n][k + 1], p.times_journey[n][k])[:],
                        t_fine
                    )
                else
                    zeros(6, 0)
                end

                if rotating
                    rotation_matrix = get_state_rotation_matrix.(2π*(t_fine .+ p.times_journey[n][k])/rotation_rate .+ visited_angle)

                    x_fine[1:6, :] = stack(rotation_matrix.*eachcol(x_fine[1:6, :]))

                    if p.id_journey[n][k + 1] != -3
                        x_target_fine[1:6, :] = stack(rotation_matrix.*eachcol(x_target_fine[1:6, :]))
                    end
                end

                lines!(ax,
                    x_fine[order[1], :],
                    x_fine[order[2], :],
                    x_fine[order[3], :],
                    color = :black
                )

                lines!(ax,
                    x_target_fine[order[1], :],
                    x_target_fine[order[2], :],
                    x_target_fine[order[3], :],
                    alpha = 0.2,
                    color = :black
                )


                # thrust_force = if p.objective_config == LoggedMassConfig()
                #     temp = exp.(p.x_nodes[n][k][7, 1:end-1])
                #     p.u_nodes[n][k][4, :] .* temp * thrust * m_scale * a_scale * 1e3
                # else
                #     p.u_nodes[n][k][4, :] * thrust * m_scale * a_scale * 1e3
                # end

                # sel = thrust_force .>= 0.01

                # thrust_vectors = stack(normalize.(eachcol(p.u_nodes[n][k][1:3, :]))).*thrust_force'

                # length_scale = 0.2

                # for i in collect(1:length(thrust_force))[sel]
                #     lines!(ax1,
                #         [p.x_nodes[n][k][1, i], p.x_nodes[n][k][1, i] + length_scale*thrust_vectors[1, i]],
                #         [p.x_nodes[n][k][2, i], p.x_nodes[n][k][2, i] + length_scale*thrust_vectors[2, i]],
                #         [p.x_nodes[n][k][3, i], p.x_nodes[n][k][3, i] + length_scale*thrust_vectors[3, i]],
                #         color = :black,
                #         alpha = 0.5,
                #         linewidth = 1,
                #     )
                # end
            end
        end
    end


    axislegend(ax1, merge = true, unique = true, orientation = :horizontal)
    

    # hidedecorations!(
    #     ax1
    # )

    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end





function plot_thrust_profile_paper(
    p::SequentialConvexProblem;
    output_file = nothing
)
    f = Figure(size = (800, 200), backgroundcolor = :white, figure_padding = 5)
    axs = []

    cs = ColorSchemes.tableau_10

    mixing = length(p.x0)
    
    n = 1
    i = 1

    deployments = sum(p.Δm0[n] .≈ -40/m_scale)

    m_fuel = if typeof(p.objective_config) == LoggedMassConfig
        m_scale*exp(p.x_nodes[n][1][7, 1]) - 500.0 - 40.0*deployments - (m_scale*exp(p.x_nodes[n][end][7, end]) - 500.0 + m_scale*p.Δm0[n][end])
    else
        m_scale*p.x_nodes[n][1][7, 1] - 500.0 - 40.0*deployments - (m_scale*p.x_nodes[n][end][7, end] - 500.0 + m_scale*p.Δm0[n][end])
    end

    m_returned = -m_scale*p.Δm0[n][end]

    t_nodes_combined = vcat(
        [p.t_nodes[n][k][1:end-1] .+ p.times_journey[n][k] for k in 1:length(p.x0[n])]...
    )
    u_nodes_combined = hcat(p.u_nodes[n]...)
    x_nodes_combined = hcat([val[:, 1:end-1] for val in p.x_nodes[n]]...)

    push!(t_nodes_combined, p.times_journey[end][end])
    u_nodes_combined = hcat(u_nodes_combined, p.u_nodes[end][end][:, end])
    x_nodes_combined = hcat(x_nodes_combined, p.x_nodes[end][end][:, end])

    thrust_force = if p.objective_config == LoggedMassConfig()
        temp = exp.(x_nodes_combined[7, :])
        u_nodes_combined[4, :] .* temp * thrust * m_scale * a_scale * 1e3
    else
        u_nodes_combined[4, :] * thrust * m_scale * a_scale * 1e3
    end
    
    ax = Axis(
        f[i, 1]; 
        xlabel = "time [MJD]", 
        # xticksvisible = false,
        # xticklabelsvisible = false,
        ylabel = "thrust", 
        xticks = [65000, 66000, 67000, 68000, 69000],
        yticks = ([0.0, 0.6], ["0.0", "max"]),
        xgridvisible = false,
        limits = (mjd_start, mjd_end, -0.05, 0.65)
    )

    push!(axs, ax)

    # vlines!(
    #     ax,
    #     convert_time_to_mjd(t_nodes_combined),
    #     color = :black,
    #     alpha = 0.1,
    # )

    Label(
        f[i, 1, Right()], 
        @sprintf("%.2f kg fuel\n%.2f kg mined", m_fuel, m_returned), 
        # font = :italic, 
        width = 120,
        justification = :left,
        padding = (10, 0, 0, 0)
    )


# text!(ax, 
#     time_positions,
#     fill(0.7*node_Δy, length(time_positions)),
#     text = time_strings,
#     align = (:center, :center),
#     fontsize = 20,
# )

    # plot_location_vertical_markers(
    #     ax, 
    #     convert_time_to_mjd(p.times_journey[n]), 
    #     p.id_journey[n], 
    #     p.Δm0[n]; 
    #     vertical_position = 1.05*maximum(thrust_force)
    # )

    stairs!(
        ax, 
        convert_time_to_mjd(t_nodes_combined), 
        thrust_force;
        step = :post,
        color = :black,
        linewidth = 1.5,
    )

    
    vlines!(
        ax,
        convert_time_to_mjd(p.times_journey[n][1]),
        color = cs[1],
        linewidth = 2,
        alpha = 0.8
    )

    
    vlines!(
        ax,
        convert_time_to_mjd(p.times_journey[n][end]),
        color = cs[1],
        linewidth = 2,
        alpha = 0.8
    )

    vlines!(
        ax,
        convert_time_to_mjd(p.times_journey[n][2:(end-1)]),
        color = cs[3],
        linewidth = 2,
        alpha = 0.8
    )




    # linkxaxes!(axs...)
    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end





function plot_team_improvements(
    submitted_files,
    reoptimized_files,
    team_names;
    output_file = nothing
)
    f = Figure(size = (800, 400), backgroundcolor = :white)
    axs = []

    cs = ColorSchemes.tableau_10

    ax = Axis(
        f[1, 1]; 
        xlabel = "improved per ship mass [kg]", 
        # yticklabelsize = 10,
        # xticksvisible = false,
        # xticklabelsvisible = false,
        # yticksvisible = false,
        # yticklabelsvisible = false,
        # ylabel = "thrust", 
        # xticks = [65000, 66000, 67000, 68000, 69000],
        # yticks = ([0.0, 0.6], ["0.0", "max"]),
        # xgridvisible = false,
        # ygridvisible = false,
        yaxisposition = :left,
        yticks = (collect(1:length(submitted_files)), team_names),
        # limits = (mjd_start, 65420, -0.05, 0.65)
    )



    for (i, (submitted_file, reoptimized_file)) in enumerate(zip(submitted_files, reoptimized_files))

        id_journey_submitted, times_journey_submitted, mined_mass_submitted, penalised_mass_submitted, groups_submitted, result_files_submitted = load_result_folders_grouping([submitted_file])

        id_journey_reoptimized, times_journey_reoptimized, mined_mass_reoptimized, penalised_mass_reoptimized, groups_reoptimized, result_files_reoptimized = load_result_folders_grouping([reoptimized_file])

        total_improved = 0.0

        y_position = i

        improved_mass = zeros(length(id_journey_submitted))
        colors = zeros(Int64, length(id_journey_submitted))

        for (i, id_journey) in enumerate(id_journey_submitted)
            k = findfirst(==(id_journey), id_journey_reoptimized)

            improved_mass[i] = mined_mass_reoptimized[k] - mined_mass_submitted[i]
            
            colors[i] = if length(findall(==(groups_submitted[i]), groups_submitted)) == 1
                1
            else
                2
            end
        end

        rangebars!(ax, 
            [i], 
            minimum(improved_mass), 
            maximum(improved_mass), 
            color = :black,
            whiskerwidth = 15, 
            direction = :x
        )

        total_improved = mean(improved_mass)

        lines!(ax,
            [total_improved, total_improved], 
            [y_position-0.35, y_position+0.35],
            color=:black,
            alpha=0.75,
            linewidth=2,
            # label=""
        )

        for (mass, color) in zip(improved_mass, colors)
            lines!(ax,
                [mass, mass], 
                [y_position-0.18, y_position+0.18],
                color=cs[color],
                label= ["self-cleaning ship", "mixed ship"][color],
                alpha=0.75,
                linewidth=2,
                # label=false
            )
        end

        push!(axs, ax)
    end

    axislegend(ax, merge = true, unique = true, orientation = :vertical, position = :rc)

    linkxaxes!(axs...)
    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return


end


function plot_graph_structure(
    p::MixedIntegerProblem;
    plot_optimal_path = false,
    plot_pruned = true,
    plot_labels = true,
    output_file = nothing,
    selection_index = 1,
    figure_size = (700, 200)
)

    stage_size = length(p.id_subset)

    stage_Δx = 0.1
    node_Δy = 0.05

    used_edge_expansion = 2.0
    # used_edge_expansion = 0.2

    unused_edge_list = []

    used_edge_lists = [[] for _ in 1:p.solutions]
    used_edge_colors = [[] for _ in 1:p.solutions]
    used_edge_styles = [[] for _ in 1:p.solutions]
    used_edge_widths = [[] for _ in 1:p.solutions]

    node_layout = []
    arc_layout = []

    # color_range = ColorSchemes.tab10
    # color_range = ColorSchemes.tableau_10
    color_range = ColorSchemes.tableau_10
    # color_range = get(colorschemes[:inferno], LinRange(1.0, 0.0, 50))
    # color_range = ColorSchemes.tableau_jewel_bright
    # color_range = ColorSchemes.tableau_miller_stone

    arc_node_offset = 0.01

    for stage in 1:p.deployment_nums[selection_index]
        for i in 1:stage_size
            push!(node_layout, ((stage - 1)*stage_Δx, -(i-1)*node_Δy))
            push!(arc_layout, ((stage - 1)*stage_Δx - arc_node_offset, -(i-1)*node_Δy))
            push!(arc_layout, ((stage - 1)*stage_Δx + arc_node_offset, -(i-1)*node_Δy))
        end
    end

    for stage in 1:p.deployment_arcs[selection_index]
        for i in 1:stage_size, j in 1:stage_size
            if !plot_pruned || ((p.deployment_cost[selection_index][i, j, stage] <= p.cost_limit) && (i != j))
                used_check = false

                temp = []

                if plot_optimal_path
                    for k in p.solutions:-1:1
                    # for k in 1:p.solutions
                        if p.id_journey_solutions[k][selection_index][stage + 1] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][stage + 2] == p.id_subset[j]
                            used_check = true

                            for l in temp
                                used_edge_widths[l][end] += used_edge_expansion
                            end

                            push!(used_edge_lists[k], (2*(i + stage_size*(stage-1)), 2*(j + stage_size*stage) - 1))
                            push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                            push!(used_edge_styles[k], :solid)
                            push!(used_edge_widths[k], 4)
                            push!(temp, k)
                        end
                    end
                end

                if !used_check
                    push!(unused_edge_list, (2*(i + stage_size*(stage-1)), 2*(j + stage_size*stage) - 1))
                end
            end
        end
    end

    collection_offset = length(node_layout)

    for i in 1:stage_size, j in 1:stage_size
        if !plot_pruned || (p.intermediate_cost[selection_index][i, j] <= p.cost_limit)
            used_check = false

            temp = []

            if plot_optimal_path
                for k in p.solutions:-1:1
                # for k in 1:p.solutions
                    if p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 1] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 2] == p.id_subset[j]
                        used_check = true

                        for l in temp
                            used_edge_widths[l][end] += used_edge_expansion
                        end
                    
                        push!(used_edge_lists[k], (2*(collection_offset - stage_size + i), 2*(collection_offset + j) - 1))
                        push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                        push!(used_edge_styles[k], :solid)
                        push!(used_edge_widths[k], 4)
                        push!(temp, k)
                    end
                end
            end

            if !used_check
                push!(unused_edge_list, (2*(collection_offset - stage_size + i), 2*(collection_offset + j) - 1))
            end
        end
    end

    for stage in (p.deployment_nums[selection_index] + 1):(p.deployment_nums[selection_index] + p.collection_nums[selection_index])
        for i in 1:stage_size
            push!(node_layout, ((stage - 1 + 1)*stage_Δx, -(i-1)*node_Δy))
            push!(arc_layout, ((stage - 1 + 1)*stage_Δx - arc_node_offset, -(i-1)*node_Δy))
            push!(arc_layout, ((stage - 1 + 1)*stage_Δx + arc_node_offset, -(i-1)*node_Δy))
        end
    end

    for stage in 1:p.collection_arcs[selection_index]
        for i in 1:stage_size, j in 1:stage_size
            if !plot_pruned || ((p.collection_cost[selection_index][i, j, stage] <= p.cost_limit) && (i != j))
                used_check = false

                temp = []

                if plot_optimal_path
                    for k in p.solutions:-1:1
                    # for k in 1:p.solutions
                        if p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 1 + stage] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 2 + stage] == p.id_subset[j]
                            used_check = true

                            for l in temp
                                used_edge_widths[l][end] += used_edge_expansion
                            end
                        
                            push!(used_edge_lists[k], (2*(collection_offset + i + stage_size*(stage-1)), 2*(collection_offset + j + stage_size*stage) - 1))
                            push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                            push!(used_edge_styles[k], :solid)
                            push!(used_edge_widths[k], 4)
                            push!(temp, k)
                        end
                    end
                end

                if !used_check
                    push!(unused_edge_list, (2*(collection_offset + i + stage_size*(stage-1)), 2*(collection_offset + j + stage_size*stage) - 1))
                end
            end
        end
    end

    unused_arc_graph = SimpleDiGraph(length(arc_layout))

    [add_edge!(unused_arc_graph, e) for e in unused_edge_list]

    fixed_node_layout(_) = node_layout
    fixed_arc_layout(_) = arc_layout

    f = Figure(size = figure_size, backgroundcolor = :white, figure_padding = 0)

    ax = Axis(
        f[1, 1];
        width = figure_size[1]
        # size = (800, 800),
    )

    graphplot!(ax, unused_arc_graph, 
        layout=fixed_arc_layout,
        # size = (800, 800),
        edge_plottype=:beziersegments,
        edge_color= (:black, 0.1),
        edge_width=2,
        arrow_show = false,
        node_size = 0,
    )

    for k in p.solutions:-1:1
    # for k in 1:p.solutions
        used_graph = SimpleDiGraph(length(arc_layout))

        [add_edge!(used_graph, e) for e in used_edge_lists[k]]

        graphplot!(ax, used_graph, 
            layout=fixed_arc_layout,
            edge_plottype=:beziersegments,
            edge_color=used_edge_colors[k],
            edge_attr = (; linestyle = used_edge_styles[k]),
            edge_width= used_edge_widths[k],
            arrow_show = false,
            node_size = 0,
        )
    end

    empty_graph = SimpleDiGraph(length(node_layout))

    graphplot!(ax, empty_graph, 
        layout=fixed_node_layout,
        ilabels = repeat(mip_problem.id_subset, p.deployment_nums[selection_index] + p.collection_nums[selection_index]),
        ilabels_fontsize = 12,
        node_size = (40, 20),
        # ilabels_fontsize = 8,
        # node_size = 15,
        node_marker=Rect;
        # node_marker=Circle;
    )

    if plot_labels
        time_positions = vcat(
            [stage_Δx*(stage-1+0.5) for stage in 1:p.deployment_arcs[selection_index]],
            [stage_Δx*(p.deployment_arcs[selection_index] + 1)],
            [stage_Δx*(stage-1+1+0.5) for stage in (p.deployment_nums[selection_index] + 1):(p.deployment_nums[selection_index] + p.collection_arcs[selection_index])]
        )

        time_strings = if p.mixing_number == 1 
            vcat(
                [latexstring("d_{i, j, $(v)}") for v in 1:p.deployment_arcs[selection_index]],
                [latexstring("m_{i, j}")],
                [latexstring("c_{i, j, $(v)}") for v in 1:p.collection_arcs[selection_index]],
            )
        else
            vcat(
                [latexstring("d_{i, j, $(v)}^{$(selection_index)}") for v in 1:p.deployment_arcs[selection_index]],
                [latexstring("m_{i, j}^{$(selection_index)}")],
                [latexstring("c_{i, j, $(v)}^{$(selection_index)}") for v in 1:p.collection_arcs[selection_index]],
            )
        end

        # display(time_positions)
        # display(time_strings)

        text!(ax, 
            time_positions,
            fill(0.7*node_Δy, length(time_positions)),
            text = time_strings,
            align = (:center, :center),
            fontsize = 20,
        )
    end

    hidedecorations!(ax)
    hidespines!(ax)

    temp = stage_Δx*(p.deployment_nums[selection_index] + p.collection_arcs[selection_index] + 1)

    xlims!(ax, [-0.04*temp, 1.04*temp])
    ylims!(ax, [-node_Δy*(length(p.id_subset) - 0.5), node_Δy*1.1])
    # ylims!(ax, [-node_Δy*(length(p.id_subset) - 0.5), node_Δy*0.5])
    

    # Colorbar(f[2, 1], limits = (0, 50), labelsize=16, colormap = color_range, vertical = false, flipaxis = false, width = 0.9*figure_size[1], label = "BIP solution rank")

    # , limits = (-1, 1), colormap = :heat,
    # label = "Temperature", vertical = false, flipaxis = false,
    # highclip = :cyan, lowclip = :red)

    # resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end


function plot_graph_structure_big(
    p::MixedIntegerProblem;
    plot_optimal_path = false,
    plot_pruned = true,
    plot_labels = true,
    output_file = nothing,
    selection_index = 1,
    figure_size = (900, 200)
)

    stage_size = length(p.id_subset)

    stage_Δx = 0.1
    node_Δy = 0.05

    # used_edge_expansion = 2.0
    used_edge_expansion = 0.2

    unused_edge_list = []

    used_edge_lists = [[] for _ in 1:p.solutions]
    used_edge_colors = [[] for _ in 1:p.solutions]
    used_edge_styles = [[] for _ in 1:p.solutions]
    used_edge_widths = [[] for _ in 1:p.solutions]

    node_layout = []
    arc_layout = []

    # color_range = ColorSchemes.tab10
    # color_range = ColorSchemes.tableau_10
    # color_range = ColorSchemes.tableau_10
    color_range = get(colorschemes[:inferno], LinRange(1.0, 0.0, 50))
    # color_range = ColorSchemes.tableau_jewel_bright
    # color_range = ColorSchemes.tableau_miller_stone
    arc_node_offset = 0.02

    for stage in 1:p.deployment_nums[selection_index]
        for i in 1:stage_size
            push!(node_layout, ((stage - 1)*stage_Δx, -(i-1)*node_Δy))
            push!(arc_layout, ((stage - 1)*stage_Δx - arc_node_offset, -(i-1)*node_Δy))
            push!(arc_layout, ((stage - 1)*stage_Δx + arc_node_offset, -(i-1)*node_Δy))
        end
    end

    for stage in 1:p.deployment_arcs[selection_index]
        for i in 1:stage_size, j in 1:stage_size
            if !plot_pruned || ((p.deployment_cost[selection_index][i, j, stage] <= p.cost_limit) && (i != j))
                used_check = false

                temp = []

                if plot_optimal_path
                    for k in p.solutions:-1:1
                    # for k in 1:p.solutions
                        if p.id_journey_solutions[k][selection_index][stage + 1] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][stage + 2] == p.id_subset[j]
                            used_check = true

                            for l in temp
                                used_edge_widths[l][end] += used_edge_expansion
                            end

                            push!(used_edge_lists[k], (2*(i + stage_size*(stage-1)), 2*(j + stage_size*stage) - 1))
                            push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                            push!(used_edge_styles[k], :solid)
                            push!(used_edge_widths[k], 4)
                            push!(temp, k)
                        end
                    end
                end

                if !used_check
                    push!(unused_edge_list, (2*(i + stage_size*(stage-1)), 2*(j + stage_size*stage) - 1))
                end
            end
        end
    end

    collection_offset = length(node_layout)

    for i in 1:stage_size, j in 1:stage_size
        if !plot_pruned || (p.intermediate_cost[selection_index][i, j] <= p.cost_limit)
            used_check = false

            temp = []

            if plot_optimal_path
                for k in p.solutions:-1:1
                # for k in 1:p.solutions
                    if p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 1] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 2] == p.id_subset[j]
                        used_check = true

                        for l in temp
                            used_edge_widths[l][end] += used_edge_expansion
                        end
                    
                        push!(used_edge_lists[k], (2*(collection_offset - stage_size + i), 2*(collection_offset + j) - 1))
                        push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                        push!(used_edge_styles[k], :solid)
                        push!(used_edge_widths[k], 4)
                        push!(temp, k)
                    end
                end
            end

            if !used_check
                push!(unused_edge_list, (2*(collection_offset - stage_size + i), 2*(collection_offset + j) - 1))
            end
        end
    end

    for stage in (p.deployment_nums[selection_index] + 1):(p.deployment_nums[selection_index] + p.collection_nums[selection_index])
        for i in 1:stage_size
            push!(node_layout, ((stage - 1 + 1)*stage_Δx, -(i-1)*node_Δy))
            push!(arc_layout, ((stage - 1 + 1)*stage_Δx - arc_node_offset, -(i-1)*node_Δy))
            push!(arc_layout, ((stage - 1 + 1)*stage_Δx + arc_node_offset, -(i-1)*node_Δy))
        end
    end

    for stage in 1:p.collection_arcs[selection_index]
        for i in 1:stage_size, j in 1:stage_size
            if !plot_pruned || ((p.collection_cost[selection_index][i, j, stage] <= p.cost_limit) && (i != j))
                used_check = false

                temp = []

                if plot_optimal_path
                    for k in p.solutions:-1:1
                    # for k in 1:p.solutions
                        if p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 1 + stage] == p.id_subset[i] && p.id_journey_solutions[k][selection_index][p.deployment_nums[selection_index] + 2 + stage] == p.id_subset[j]
                            used_check = true

                            for l in temp
                                used_edge_widths[l][end] += used_edge_expansion
                            end
                        
                            push!(used_edge_lists[k], (2*(collection_offset + i + stage_size*(stage-1)), 2*(collection_offset + j + stage_size*stage) - 1))
                            push!(used_edge_colors[k], (color_range[(k-1) % length(color_range) + 1], 1.0))
                            push!(used_edge_styles[k], :solid)
                            push!(used_edge_widths[k], 4)
                            push!(temp, k)
                        end
                    end
                end

                if !used_check
                    push!(unused_edge_list, (2*(collection_offset + i + stage_size*(stage-1)), 2*(collection_offset + j + stage_size*stage) - 1))
                end
            end
        end
    end

    unused_arc_graph = SimpleDiGraph(length(arc_layout))

    [add_edge!(unused_arc_graph, e) for e in unused_edge_list]

    fixed_node_layout(_) = node_layout
    fixed_arc_layout(_) = arc_layout
    
    f = Figure(size = figure_size, backgroundcolor = :white, figure_padding = 0)

    ax = Axis(
        f[1, 1];
        width = figure_size[1]
        # size = (800, 800),
    )


    graphplot!(ax, unused_arc_graph, 
        layout=fixed_arc_layout,
        # size = (800, 800),
        edge_plottype=:beziersegments,
        edge_color= (:black, 0.08),
        edge_width=2,
        arrow_show = false,
        node_size = 0,
    )

    for k in p.solutions:-1:1
    # for k in 1:p.solutions
        used_graph = SimpleDiGraph(length(arc_layout))

        [add_edge!(used_graph, e) for e in used_edge_lists[k]]

        graphplot!(ax, used_graph, 
            layout=fixed_arc_layout,
            edge_plottype=:beziersegments,
            edge_color=used_edge_colors[k],
            edge_attr = (; linestyle = used_edge_styles[k]),
            edge_width= used_edge_widths[k],
            arrow_show = false,
            node_size = 0,
        )
    end


    empty_graph = SimpleDiGraph(length(node_layout))

    graphplot!(ax, empty_graph, 
        layout=fixed_node_layout,
        ilabels = repeat(mip_problem.id_subset, p.deployment_nums[selection_index] + p.collection_nums[selection_index]),
        # ilabels_fontsize = 12,
        # node_size = 20,
        ilabels_fontsize = 7,
        node_size = (22, 10),
        node_marker=Rect;
    )

    # if plot_labels
    #     time_positions = vcat(
    #         [stage_Δx*(stage-1+0.5) for stage in 1:p.deployment_arcs[selection_index]],
    #         [stage_Δx*(p.deployment_arcs[selection_index] + 1)],
    #         [stage_Δx*(stage-1+1+0.5) for stage in (p.deployment_nums[selection_index] + 1):(p.deployment_nums[selection_index] + p.collection_arcs[selection_index])]
    #     )

    #     time_strings = if p.mixing_number == 1 
    #         vcat(
    #             [latexstring("d_{i, j, $(v)}") for v in 1:p.deployment_arcs[selection_index]],
    #             [latexstring("m_{i, j}")],
    #             [latexstring("c_{i, j, $(v)}") for v in 1:p.collection_arcs[selection_index]],
    #         )
    #     else
    #         vcat(
    #             [latexstring("d_{i, j, $(v)}^{$(selection_index)}") for v in 1:p.deployment_arcs[selection_index]],
    #             [latexstring("m_{i, j}^{$(selection_index)}")],
    #             [latexstring("c_{i, j, $(v)}^{$(selection_index)}") for v in 1:p.collection_arcs[selection_index]],
    #         )
    #     end

    #     # display(time_positions)
    #     # display(time_strings)

    #     text!(ax, 
    #         time_positions,
    #         fill(0.7*node_Δy, length(time_positions)),
    #         text = time_strings,
    #         align = (:center, :center),
    #         fontsize = 20,
    #     )
    # end

    hidedecorations!(ax)
    hidespines!(ax)

    temp = stage_Δx*(p.deployment_nums[selection_index] + p.collection_arcs[selection_index] + 1)

    xlims!(ax, [-0.025*temp, 1.025*temp])
    ylims!(ax, [-node_Δy*(length(p.id_subset) - 0.5), node_Δy*1.1])
    # ylims!(ax, [-node_Δy*(length(p.id_subset) - 0.5), node_Δy*0.5])
    

    Colorbar(f[2, 1], limits = (0, 50), labelsize=16, colormap = color_range, vertical = false, flipaxis = false, width = 0.98*figure_size[1], label = "BIP solution rank")

    # , limits = (-1, 1), colormap = :heat,
    # label = "Temperature", vertical = false, flipaxis = false,
    # highclip = :cyan, lowclip = :red)

    # resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end


function plot_discretization_comparison(
    p1::SequentialConvexProblem,
    p2::SequentialConvexProblem,
    p3::SequentialConvexProblem;
    output_file = nothing
)
    f = Figure(size = (800, 360), backgroundcolor = :white)
    axs = []

    cs = ColorSchemes.tableau_10

    for (i, p) in enumerate([p1, p2, p3])
        mixing = length(p.x0)
        
        n = 1

        deployments = sum(p.Δm0[n] .≈ -40/m_scale)

        m_fuel = if typeof(p.objective_config) == LoggedMassConfig
            m_scale*exp(p.x_nodes[n][1][7, 1]) - 500.0 - 40.0*deployments - (m_scale*exp(p.x_nodes[n][end][7, end]) - 500.0 + m_scale*p.Δm0[n][end])
        else
            m_scale*p.x_nodes[n][1][7, 1] - 500.0 - 40.0*deployments - (m_scale*p.x_nodes[n][end][7, end] - 500.0 + m_scale*p.Δm0[n][end])
        end

        m_returned = -m_scale*p.Δm0[n][end]

        t_nodes_combined = vcat(
            [p.t_nodes[n][k][1:end-1] .+ p.times_journey[n][k] for k in 1:length(p.x0[n])]...
        )
        u_nodes_combined = hcat(p.u_nodes[n]...)
        x_nodes_combined = hcat([val[:, 1:end-1] for val in p.x_nodes[n]]...)

        push!(t_nodes_combined, p.times_journey[end][end])
        u_nodes_combined = hcat(u_nodes_combined, p.u_nodes[end][end][:, end])
        x_nodes_combined = hcat(x_nodes_combined, p.x_nodes[end][end][:, end])

        thrust_force = if p.objective_config == LoggedMassConfig()
            temp = exp.(x_nodes_combined[7, :])
            u_nodes_combined[4, :] .* temp * thrust * m_scale * a_scale * 1e3
        else
            u_nodes_combined[4, :] * thrust * m_scale * a_scale * 1e3
        end
        
        ax = Axis(
            f[i, 1]; 
            xlabel = i == 3 ? "time [MJD]" : "", 
            xticksvisible = i == 3,
            xticklabelsvisible =  i == 3,
            ylabel = "thrust", 
            xticks = [64300, 64400, 64500, 64600, 64700, 64800, 64900, 65000, 65100, 65200, 65300, 65400, 65500],
            yticks = ([0.0, 0.6], ["0.0", "max"]),
            xgridvisible = false,
            limits = (64300, 65420, -0.05, 0.65)
        )

        push!(axs, ax)

        vlines!(
            ax,
            convert_time_to_mjd(t_nodes_combined),
            color = :black,
            alpha = 0.1,
        )

        text1 = [
            "Fixed Time",
            "Fixed Time",
            "Dynamic Time",
        ]

        text2 = [
            "Fixed Discretization",
            "Dynamic Discretization",
            "Dynamic Discretization",
        ]

        text = "$(text1[i])\n$(text2[i])\n"

        Label(
            f[i, 1, Right()], 
            @sprintf("%s\n%s\n%.2f kg fuel\n%.2f kg mined", text1[i], text2[i], m_fuel, m_returned), 
            # font = :italic, 
            width = 150,
            justification = :left,
            padding = (10, 0, 0, 0)
        )


    # text!(ax, 
    #     time_positions,
    #     fill(0.7*node_Δy, length(time_positions)),
    #     text = time_strings,
    #     align = (:center, :center),
    #     fontsize = 20,
    # )

        # plot_location_vertical_markers(
        #     ax, 
        #     convert_time_to_mjd(p.times_journey[n]), 
        #     p.id_journey[n], 
        #     p.Δm0[n]; 
        #     vertical_position = 1.05*maximum(thrust_force)
        # )

        stairs!(
            ax, 
            convert_time_to_mjd(t_nodes_combined), 
            thrust_force;
            step = :post,
            color = :black,
            linewidth = 1.5,
        )

        
        vlines!(
            ax,
            convert_time_to_mjd(p.times_journey[n][1:4]),
            color = cs[[1, 3, 3, 3]],
            linewidth = 2,
            alpha = 0.8
        )
    end

    linkxaxes!(axs...)
    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end




function plot_convergence_comparison(
    mip_problem_objectives, 
    scp_problem_objectives;
    output_file = nothing
)

    f = Figure(size = (800, 300), backgroundcolor = :white)

    cs = ColorSchemes.tableau_10
    
    ax = Axis(
        f[1, 1]; 
        xlabel = "BIP solution rank", 
        ylabel = "SCP per ship mass [kg]", 
        xticks = 0:5:50,
        limits = (0, 51, 720, 790)
        # xticks = 0:10:50
    )

    for i in 1:length(mip_problem_objectives)
        scatter!(ax,
            mip_problem_objectives[i],
            scp_problem_objectives[i],
            label = "Iteration $(i)"
        )
    end

    axislegend(ax, unique = true, merge = true, orientation = :horizontal)

    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end



function plot_gtoc12_problem(
    p::SequentialConvexProblem
)

    cs = ColorSchemes.tableau_10

    f = Figure(size = (1000, 480), backgroundcolor = :white, figure_padding = 0)

    ax1 = Axis3(
        f[1, 1]; 
        # xlabel = "x [AU]", 
        # ylabel = "y [AU]", 
        limits = 0.84.*(-3.0, 3.0, -3.0, 3.0, -0.75, 0.75),
        azimuth = -π/2 + 0.2,
        elevation = π/10 + 0.03,
        aspect = (3, 3, 1),
    )

    scatter!(ax1,
        [0.0],
        [0.0],
        [0.0],
        color = cs[6],
        markersize = 10,
        label = "Sun"
    )

    for ax in [ax1]
        ν_plot = collect(LinRange(0.0, 2π, 400)) 

        for (i, planet_classical) in enumerate(eachcol(planets_classical))
            temp = repeat(planet_classical, 1, 400)
            temp[6, :] .= ν_plot
            temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            lines!(ax,
                temp[1, :],
                temp[2, :],
                temp[3, :],
                # linestyle = :dashdot,
                linestyle = :solid,
                linewidth = 1,
                color = cs[[2, 1, 3][i]],
                label = ["Venus Orbit", "Earth Orbit", "Mars Orbit"][i]
            )
        end

        for (i, asteroid_classical) in enumerate(eachcol(asteroids_classical[:, [15184, 12286]]))
            temp = repeat(asteroid_classical, 1, 400)
            temp[6, :] .= ν_plot
            temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            lines!(ax,
                temp[1, :],
                temp[2, :],
                temp[3, :],
                linestyle = :solid,
                # linestyle = :dashdot,
                linewidth = 1,
                alpha = 0.2,
                color = :black,
                label = "Asteroid Orbits"
            )
        end


        for n in 1:p.mixing_number
            for k in 1:p.segment_number[n]
                t_fine = collect(p.t_nodes[n][k][1]:1.0*day_scale:p.t_nodes[n][k][end])

                if !(t_fine[end] ≈ p.t_nodes[n][k][end])
                    push!(t_fine, p.t_nodes[n][k][end])
                end

                x_fine = integrate_trajectory(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    t_fine;
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                lines!(ax,
                    x_fine[1, :],
                    x_fine[2, :],
                    x_fine[3, :],
                    color = :black,
                    alpha = 0.5,
                    linewidth = 2,
                    linestyle = :solid,
                    label = "Mining Ship\nTrajectory"
                )

                scatter!(ax,
                    p.x0[n][k][1],
                    p.x0[n][k][2],
                    p.x0[n][k][3],
                    color = :black,
                    label = "Rendezvous\nEvent"
                )

                if k == p.segment_number[n]
                    scatter!(ax,
                        p.xf[n][k][1],
                        p.xf[n][k][2],
                        p.xf[n][k][3],
                        color = :black
                    )
                end

                text = [
                    "Departure",
                    "Asteroid 1\nDeployment",
                    "Asteroid 2\nDeployment",
                    "Asteroid 1\nCollection",
                    "Asteroid 2\nCollection",
                    "Arrival"
                ]

                
                offsets = [
                    (-7.5, 0),
                    (5, 0),
                    (5, 0),
                    (5, 0),
                    (-5, 0),
                    (-5, 0),
                ]

                alignment = [
                    (:right, :bottom),
                    (:left, :bottom),
                    (:left, :bottom),
                    (:left, :bottom),
                    (:right, :bottom),
                    (:right, :bottom)
                ]
                
                text!(ax, 
                    p.x0[n][k][1],
                    p.x0[n][k][2],
                    p.x0[n][k][3],
                    text = text[k],
                    align = alignment[k],
                    # fontsize = 12,
                    offset = offsets[k]
                )

                if k == p.segment_number[n]
                    text!(ax, 
                        p.xf[n][k][1],
                        p.xf[n][k][2],
                        p.xf[n][k][3],
                        text = text[end],
                        align = alignment[end],
                        # fontsize = 12,
                        offset = offsets[end],
                    )
                end

            end
        end


        hidedecorations!(ax)
        hidespines!(ax)

        # axislegend(ax, unique = true, merge = true)

    end

    
    f[1, 2] = Legend(f, ax1, framevisible = false, unique = true, merge = true, labelsize=16)


    resize_to_layout!(f)
    display(f)

    save("output/plots/gtoc12_explainer.png", f)

    return
end




function plot_scp_details(
    p::SequentialConvexProblem
)

    cs = ColorSchemes.tableau_10
    cs2 = ColorSchemes.tableau_10[[4, 5, 6, 7, 9]]
    # cs2 = [:black, :black, :black, :black, :black]
    # cs2 = ColorSchemes.Archambault

    f = Figure(size = (1000, 480), backgroundcolor = :white)
    # f = Figure(size = (1000, 680), backgroundcolor = :white)

    # ax1 = Axis(
    #     f[1, 1]; 
    #     # xlabel = "x [AU]", 
    #     # ylabel = "y [AU]", 
    #     limits = 1.05.*(-3.0, 3.0, -3.0, 3.0),
    #     aspect = 1
    # )

    ax1 = Axis3(
        f[1, 1]; 
        # xlabel = "x [AU]", 
        # ylabel = "y [AU]", 
        limits = 0.84.*(-3.0, 3.0, -3.0, 3.0, -0.75, 0.75),
        azimuth = -π/2,
        elevation = π/10,
        aspect = (3, 3, 1)
    )

    scatter!(ax1,
        [0.0],
        [0.0],
        color = cs[6],
        markersize = 10,
        label = "Sun"
    )

    for ax in [ax1]
        ν_plot = collect(LinRange(0.0, 2π, 400)) 

        for (i, planet_classical) in enumerate(eachcol(planets_classical))
            temp = repeat(planet_classical, 1, 400)
            temp[6, :] .= ν_plot
            temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            lines!(ax,
                temp[1, :],
                temp[2, :],
                temp[3, :],
                linestyle = :solid,
                linewidth = 1,
                # linestyle = :dashdot,
                color = cs[[2, 1, 3][i]],
                label = ["Venus Orbit", "Earth Orbit", "Mars Orbit"][i]
            )
        end

        for (i, asteroid_classical) in enumerate(eachcol(asteroids_classical[:, [15184, 12286]]))
            temp = repeat(asteroid_classical, 1, 400)
            temp[6, :] .= ν_plot
            temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            lines!(ax,
                temp[1, :],
                temp[2, :],
                temp[3, :],
                linestyle = :solid,
                # linestyle = :dashdot,
                linewidth = 1,
                alpha = 0.2,
                color = :black,
                label = "Asteroid Orbits"
            )
        end


        for n in 1:p.mixing_number
            for k in 1:p.segment_number[n]
                t_fine = collect(p.t_nodes[n][k][1]:1.0*day_scale:p.t_nodes[n][k][end])

                if !(t_fine[end] ≈ p.t_nodes[n][k][end])
                    push!(t_fine, p.t_nodes[n][k][end])
                end

                x_fine = integrate_trajectory(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    t_fine;
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                scatter!(
                    ax,
                    p.x_nodes[n][k][1, :],
                    p.x_nodes[n][k][2, :],
                    p.x_nodes[n][k][3, :],
                    color = cs2[k]
                )

                text = [
                    "Leg 1",
                    "Leg 2",
                    "Leg 3",
                    "Leg 4",
                    "Leg 5",
                ]
                
                lines!(ax,
                    x_fine[1, :],
                    x_fine[2, :],
                    x_fine[3, :],
                    # color = :black,
                    color = cs2[k],
                    alpha = 1.0,
                    linestyle = :solid,
                    linewidth = 2.0,
                    label = text[k]
                )


                # text = [
                #     "Earth",
                #     "Asteroid 1",
                #     "Asteroid 2",
                #     "Asteroid 1",
                #     "Asteroid 2",
                #     "Earth"
                # ]

                # # alignment = [
                # #     (:left, :bottom),
                # #     (:left, :bottom),
                # #     (:right, :top),
                # #     (:left, :bottom),
                # #     (:right, :bottom),
                # #     (:right, :bottom)
                # # ]

                
                # alignment = [
                #     (:left, :bottom),
                #     (:left, :bottom),
                #     (:left, :bottom),
                #     (:left, :bottom),
                #     (:right, :bottom),
                #     (:right, :bottom)
                # ]

                # # offsets = [
                # #     (5, 0),
                # #     (5, 0),
                # #     (-5, 0),
                # #     (5, 0),
                # #     (-5, 0),
                # #     (-5, 0),
                # # ]

                # offsets = [
                #     (5, 0),
                #     (5, 0),
                #     (5, 0),
                #     (5, 0),
                #     (-5, 0),
                #     (-5, 0),
                # ]
                
                # text!(ax, 
                #     p.x0[n][k][1],
                #     p.x0[n][k][2],
                #     p.x0[n][k][3],
                #     text = text[k],
                #     align = alignment[k],
                #     # fontsize = 12,
                #     offset = offsets[k]
                # )

                # if k == p.segment_number[n]
                #     text!(ax, 
                #         p.xf[n][k][1],
                #         p.xf[n][k][2],
                #         p.xf[n][k][3],
                #         text = text[end],
                #         align = alignment[end],
                #         # fontsize = 12,
                #         offset = offsets[end],
                #     )
                # end

                # temp = Int64(round(size(p.x_nodes[n][k], 2) *0.7))

                # text = [
                #     "Leg 1\nEarth->Asteroid 1",
                #     "Leg 2\nAsteroid 1->Asteroid 2",
                #     "Leg 3\nAsteroid 2->Asteroid 1",
                #     "Leg 4\nAsteroid 1->Asteroid 2",
                #     "Leg 5\nAsteroid 2->Earth",
                # ]


                # text!(ax,
                #     p.x_nodes[n][k][1, temp],
                #     p.x_nodes[n][k][2, temp],
                #     text = text[k]
                # )



            end

            for k in 1:p.segment_number[n]
                scatter!(ax,
                    p.x0[n][k][1],
                    p.x0[n][k][2],
                    p.x0[n][k][3],
                    color = :black,
                    label = "Rendezvous\nEvent"
                )

                if k == p.segment_number[n]
                    scatter!(ax,
                        p.xf[n][k][1],
                        p.xf[n][k][2],
                        p.xf[n][k][3],
                        color = :black,
                        label = "Rendezvous\nEvent"
                    )
                end
            end
        end


        hidedecorations!(ax)
        hidespines!(ax)

        # axislegend(ax, unique = true, merge = true)

    end

    f[1, 2] = Legend(f, ax1, framevisible = false, unique = true, merge = true, labelsize=16)

    resize_to_layout!(f)
    display(f)

    save("output/plots/gtoc12_scp_explainer.png", f)

    return
end



function plot_trajectory_final(
    p1::SequentialConvexProblem,
    p2::SequentialConvexProblem;
    output_file = nothing
)
    f = Figure(size = (800, 550), backgroundcolor = :white, figure_padding = 2)

    ax1 = Axis(
        f[1:4, 1]; 
        xlabel = "x [AU]", 
        ylabel = "y [AU]", 
        limits = (-3.2, 3.2, -3.2, 3.2),
        # limits = (-0.25, 3.2, -1.1, 1.1),
        # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        # yticks = [-1.0, -0.5, 0.0, 0.5, 1.0],
        # xticklabelsvisible = false,
        # limits = (0.5, 1.5, -3.0, -2.0),
        # aspect = 3.7/3
        # aspect = 1.568
    )

    ax2 = Axis(
        f[1:4, 2]; 
        xlabel = "x [AU]", 
        # ylabel = "y [AU]", 
        limits = (-3.2, 3.2, -3.2, 3.2),
        # limits = (2.2, 3.2, -0.5, 0.5),
        # limits = (-0.25, 3.2, -1.1, 1.1),
        # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        # yticks = [-1.0, -0.5, 0.0, 0.5, 1.0],
        # xticklabelsvisible = false,
        yticklabelsvisible = false,
        # limits = (0.5, 1.5, -3.0, -2.0),
        # aspect = 3.7/3
        # aspect = 1.568
    )

    ax3 = Axis(
        f[5, 1]; 
        xlabel = "ship ID", 
        ylabel = "rendezvous number", 
        # limits = (-3.2, 3.2, -0.5, 0.5),
        # limits = (2.2, 3.2, -0.5, 0.5),
        # limits = (-0.25, 3.2, -0.35, 0.35),
        # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        # yticks = [-0.2, 0.0, 0.2],
        # yticklabelsvisible = false,
        # limits = (0.5, 1.5, -3.0, -2.0),
        # aspect = 3.7/3
        # aspect = 1.32
    )

    ax4 = Axis(
        f[5, 2]; 
        xlabel = "ship ID", 
        yticklabelsvisible = false,
        # limits = (-3.2, 3.2, -0.5, 0.5),
        # ylabel = "z [AU]", 
        # limits = (2.2, 3.2, -0.5, 0.5),
        # limits = (-0.25, 3.2, -0.35, 0.35),
        # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        # yticks = [-0.2, 0.0, 0.2],
        # limits = (0.5, 1.5, -3.0, -2.0),
        # aspect = 3.7/3
        # aspect = 1.32
    )


    # rotation_rate = if rotating
    #     visited_ids = setdiff(unique(vcat(p.id_journey[solution_indices]...)), [0, -3])

    #     a_visited = asteroids_classical[1, visited_ids]

    #     period_visited = 2π * sqrt.(a_visited.^3/μ)

    #     mean(period_visited)
    # else
    #     0.0
    # end

    # visited_angle = -1.25

    # scatter!(ax1,
    #     [0.0],
    #     [0.0],
    #     [0.0],
    #     color = :black
    # )

    # ν_plot = collect(LinRange(0.0, 2π, 400)) 

    for (ax, ax_aux, order, p) in zip([ax1, ax2], [ax3, ax4], [[1, 2, 3], [1, 2, 3]], [p1, p2])
        
        solution_indices = collect(1:p.mixing_number)

        # t_planet = convert_mjd_to_time.(collect(mjd_start:10:mjd_end))
        
        ν_plot = collect(LinRange(0.0, 2π, 400)) 

        for (i, planet_classical) in enumerate(eachcol(planets_classical))
            temp = repeat(planet_classical, 1, 400)
            temp[6, :] .= ν_plot
            temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            lines!(ax,
                temp[order[1], :],
                temp[order[2], :],
                temp[order[3], :],
                # linestyle = :dashdot,
                color = Makie.wong_colors()[[2, 1, 6][i]]
            )
        end

        for n in solution_indices
            for k in 1:p.segment_number[n]
                x0 = p.x0[n][k]
                xf = p.xf[n][k]

                # if rotating
                #     x0 = get_state_rotation_matrix(2π*(p.times_journey[n][k])/rotation_rate + visited_angle)*x0[1:6]
                #     xf = get_state_rotation_matrix(2π*(p.times_journey[n][k+1])/rotation_rate + visited_angle)*xf[1:6]
                # end

                scatter!(ax,
                    x0[order[1]],
                    x0[order[2]],
                    x0[order[3]],
                )

                if k == p.segment_number[n]
                    scatter!(ax,
                        xf[order[1]],
                        xf[order[2]],
                        xf[order[3]],
                    )
                end

                t_fine = collect(p.t_nodes[n][k][1]:1.0*day_scale:p.t_nodes[n][k][end])

                if !(t_fine[end] ≈ p.t_nodes[n][k][end])
                    push!(t_fine, p.t_nodes[n][k][end])
                end

                x_fine = integrate_trajectory(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    t_fine;
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                # x_target_fine = if p.id_journey[n][k + 1] != -3
                #     integrate_trajectory(
                #         ephermeris_cartesian_from_id(p.id_journey[n][k + 1], p.times_journey[n][k])[:],
                #         t_fine
                #     )
                # else
                #     zeros(6, 0)
                # end

                # if rotating
                #     rotation_matrix = get_state_rotation_matrix.(2π*(t_fine .+ p.times_journey[n][k])/rotation_rate .+ visited_angle)

                #     x_fine[1:6, :] = stack(rotation_matrix.*eachcol(x_fine[1:6, :]))

                #     # if p.id_journey[n][k + 1] != -3
                #     #     x_target_fine[1:6, :] = stack(rotation_matrix.*eachcol(x_target_fine[1:6, :]))
                #     # end
                # end

                lines!(ax,
                    x_fine[order[1], :],
                    x_fine[order[2], :],
                    x_fine[order[3], :],
                    color = :black,
                    alpha = 0.25
                )

                # lines!(ax,
                #     x_target_fine[order[1], :],
                #     x_target_fine[order[2], :],
                #     x_target_fine[order[3], :],
                #     alpha = 0.2,
                #     color = :black
                # )


                # thrust_force = if p.objective_config == LoggedMassConfig()
                #     temp = exp.(p.x_nodes[n][k][7, 1:end-1])
                #     p.u_nodes[n][k][4, :] .* temp * thrust * m_scale * a_scale * 1e3
                # else
                #     p.u_nodes[n][k][4, :] * thrust * m_scale * a_scale * 1e3
                # end

                # sel = thrust_force .>= 0.01

                # thrust_vectors = stack(normalize.(eachcol(p.u_nodes[n][k][1:3, :]))).*thrust_force'

                # length_scale = 0.2

                # for i in collect(1:length(thrust_force))[sel]
                #     lines!(ax1,
                #         [p.x_nodes[n][k][1, i], p.x_nodes[n][k][1, i] + length_scale*thrust_vectors[1, i]],
                #         [p.x_nodes[n][k][2, i], p.x_nodes[n][k][2, i] + length_scale*thrust_vectors[2, i]],
                #         [p.x_nodes[n][k][3, i], p.x_nodes[n][k][3, i] + length_scale*thrust_vectors[3, i]],
                #         color = :black,
                #         alpha = 0.5,
                #         linewidth = 1,
                #     )
                # end
            end
        end

        # display(p.id_journey)

        # temp = vcat([fill(i, length(k) - 2) for (i, k) in enumerate(p.id_journey)]...)

        # display(temp)

        # colors = vcat([
        #     [count(==(j), k[2:end-1]) for j in k[2:end-1]]
        #     for k in p.id_journey
        # ]...)

        # display(colors)


        # barplot!(ax_aux,
        #     temp,
        #     fill(1, length(temp)),
        #     stack = colors,
        #     color = Makie.wong_colors()[[2, 1]][colors]
        # )

        
        temp = vcat([fill(i, 2) for (i, k) in enumerate(p.id_journey)]...)

        display(temp)

        temp2 = [[count(==(j), k[2:end-1]) for j in k[2:end-1]] for k in p.id_journey]

        heights = vcat(
            [vcat(count(==(1), k), count(==(2), k)) for k in temp2]...
        )

        colors = vcat([[1, 2] for (i, k) in enumerate(p.id_journey)]...)


        barplot!(ax_aux,
            temp,
            heights,
            stack = colors,
            color = ColorSchemes.tab10[[2, 1]][colors]
        )





    end


    # axislegend(ax1, merge = true, unique = true, orientation = :horizontal)
    

    # hidedecorations!(
    #     ax1
    # )

    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end


function plot_trajectory_and_thrust_profile_paper(
    p::SequentialConvexProblem;
    output_file = nothing,
    solution_indices = nothing,
    label_text = ""
)

    f = Figure(size = (550, 700), backgroundcolor = :white, figure_padding = 2)

    if isnothing(solution_indices)
        solution_indices = collect(1:p.mixing_number)
    end
    
    cs = ColorSchemes.tab10
    # cs = ColorSchemes.tab20
    # cs = ColorSchemes.tableau_miller_stone

    for i in 1:p.mixing_number

        ax1 = Axis(
            f[1:4, i]; 
            xlabel = "x [AU]", 
            ylabel = i==1 ? "y [AU]" : "", 
            limits = (-3.2, 3.2, -3.2, 3.2),
            # limits = (-0.25, 3.2, -1.1, 1.1),
            # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
            # yticks = [-1.0, -0.5, 0.0, 0.5, 1.0],
            yticklabelsvisible = i==1,
            # limits = (0.5, 1.5, -3.0, -2.0),
            # aspect = 3.7/3
            # aspect = 1.0
        )

        ax2 = Axis(
            f[5, i]; 
            xlabel = "time [MJD]", 
            # xticksvisible = false,
            # xticklabelsvisible = false,
            yticklabelsvisible = i==1,
            ylabel = i==1 ? "thrust" : "", 
            xticks = [65000, 66000, 67000, 68000, 69000],
            yticks = ([0.0, 0.6], ["0.0", "max"]),
            xgridvisible = false,
            limits = (64300, 69900, -0.05, 0.65)
        )

        # scatter!(ax1,
        #     [0.0],
        #     [0.0],
        #     [0.0],
        #     color = :black
        # )

        ν_plot = collect(LinRange(0.0, 2π, 400)) 

        for (i, planet_classical) in enumerate(eachcol(planets_classical))
            temp = repeat(planet_classical, 1, 400)
            temp[6, :] .= ν_plot
            temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            lines!(ax1,
                temp[1, :],
                temp[2, :],
                temp[3, :],
                # linestyle = :dashdot,
                color = Makie.wong_colors()[[2, 1, 6][i]],
                # alpha = 0.35,
                label = ["Venus Orbit", "Earth Orbit", "Mars Orbit"][i]
            )
        end

        for n in [i]
            temp = deepcopy(p.id_journey[n])
            temp[1] = -3

            chosen_colors = cs[[findfirst(==(i), temp) for i in temp]]

            for k in 1:p.segment_number[n]
                x0 = p.x0[n][k]
                xf = p.xf[n][k]

                t_fine = collect(p.t_nodes[n][k][1]:1.0*day_scale:p.t_nodes[n][k][end])

                if !(t_fine[end] ≈ p.t_nodes[n][k][end])
                    push!(t_fine, p.t_nodes[n][k][end])
                end

                x_fine = integrate_trajectory(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    t_fine;
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                lines!(ax1,
                    x_fine[1, :],
                    x_fine[2, :],
                    x_fine[3, :],
                    color = :black,
                    alpha = 0.8,
                )
                
                scatter!(ax1,
                    x0[1],
                    x0[2],
                    x0[3],
                    color = chosen_colors[k],
                    markersize = 10,
                )

                if k == p.segment_number[n]
                    scatter!(ax1,
                        xf[1],
                        xf[2],
                        xf[3],
                        color = chosen_colors[k+1],
                        markersize = 10,
                    )
                end

                # lines!(ax1,
                #     x_target_fine[1, :],
                #     x_target_fine[2, :],
                #     x_target_fine[3, :],
                #     alpha = 0.2,
                #     color = :black
                # )

                thrust_force = if p.objective_config == LoggedMassConfig()
                    temp = exp.(p.x_nodes[n][k][7, 1:end-1])
                    p.u_nodes[n][k][4, :] .* temp * thrust * m_scale * a_scale * 1e3
                else
                    p.u_nodes[n][k][4, :] * thrust * m_scale * a_scale * 1e3
                end

                sel = thrust_force .>= 0.01

                thrust_vectors = stack(normalize.(eachcol(p.u_nodes[n][k][1:3, :]))).*thrust_force'

                length_scale = 0.2

                for i in collect(1:length(thrust_force))[sel]
                    lines!(ax1,
                        [p.x_nodes[n][k][1, i], p.x_nodes[n][k][1, i] + length_scale*thrust_vectors[1, i]],
                        [p.x_nodes[n][k][2, i], p.x_nodes[n][k][2, i] + length_scale*thrust_vectors[2, i]],
                        [p.x_nodes[n][k][3, i], p.x_nodes[n][k][3, i] + length_scale*thrust_vectors[3, i]],
                        color = :black,
                        alpha = 0.5,
                        linewidth = 1,
                    )
                end
            end

                
            deployments = sum(p.Δm0[n] .≈ -40/m_scale)

            m_fuel = if typeof(p.objective_config) == LoggedMassConfig
                m_scale*exp(p.x_nodes[n][1][7, 1]) - 500.0 - 40.0*deployments - (m_scale*exp(p.x_nodes[n][end][7, end]) - 500.0 + m_scale*p.Δm0[n][end])
            else
                m_scale*p.x_nodes[n][1][7, 1] - 500.0 - 40.0*deployments - (m_scale*p.x_nodes[n][end][7, end] - 500.0 + m_scale*p.Δm0[n][end])
            end

            m_returned = -m_scale*p.Δm0[n][end]

            t_nodes_combined = vcat(
                [p.t_nodes[n][k][1:end-1] .+ p.times_journey[n][k] for k in 1:length(p.x0[n])]...
            )
            u_nodes_combined = hcat(p.u_nodes[n]...)
            x_nodes_combined = hcat([val[:, 1:end-1] for val in p.x_nodes[n]]...)

            push!(t_nodes_combined, p.times_journey[end][end])
            u_nodes_combined = hcat(u_nodes_combined, p.u_nodes[end][end][:, end])
            x_nodes_combined = hcat(x_nodes_combined, p.x_nodes[end][end][:, end])

            thrust_force = if p.objective_config == LoggedMassConfig()
                temp = exp.(x_nodes_combined[7, :])
                u_nodes_combined[4, :] .* temp * thrust * m_scale * a_scale * 1e3
            else
                u_nodes_combined[4, :] * thrust * m_scale * a_scale * 1e3
            end

            stairs!(
                ax2, 
                convert_time_to_mjd(t_nodes_combined), 
                thrust_force;
                step = :post,
                color = :black,
                linewidth = 1.5,
            )
        
            
            vlines!(
                ax2,
                convert_time_to_mjd(p.times_journey[n][1]),
                color = cs[1],
                linewidth = 2,
                alpha = 0.9
            )
        
            
            vlines!(
                ax2,
                convert_time_to_mjd(p.times_journey[n][end]),
                color = cs[1],
                linewidth = 2,
                alpha = 0.9
            )
        
            vlines!(
                ax2,
                convert_time_to_mjd(p.times_journey[n][2:(end-1)]),
                color = chosen_colors[2:end-1],
                linewidth = 2,
                alpha = 0.9
            )
        end
    end




    # Label(
    #     f[1:4, 3], 
    #     label_text, 
    #     # font = :italic, 
    #     width = 150,
    #     justification = :left,
    #     padding = (-50, 0, 0, 0)
    # )


    # axislegend(ax1, merge = true, unique = true, orientation = :horizontal)

    # linkxaxes!(axs...)
    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end



function plot_trajectory_paper_2(
    p::SequentialConvexProblem;
    output_file = nothing,
    solution_indices = nothing,
    label_text = ""
)

    f = Figure(size = (500, 500), backgroundcolor = :white, figure_padding = 2)

    if isnothing(solution_indices)
        solution_indices = collect(1:p.mixing_number)
    end
    


    # cs = ColorSchemes.tab10
    cs = ColorSchemes.tab20
    # cs = ColorSchemes.tableau_miller_stone

    for i in 1:p.mixing_number

        ax1 = Axis(
            f[1:4, i]; 
            xlabel = "x [AU]", 
            ylabel = i==1 ? "y [AU]" : "", 
            limits = (-3.2, 3.2, -3.2, 3.2),
            # limits = (-0.25, 3.2, -1.1, 1.1),
            # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
            # yticks = [-1.0, -0.5, 0.0, 0.5, 1.0],
            yticklabelsvisible = i==1,
            # limits = (0.5, 1.5, -3.0, -2.0),
            # aspect = 3.7/3
            # aspect = 1.0
        )

        # scatter!(ax1,
        #     [0.0],
        #     [0.0],
        #     [0.0],
        #     color = :black
        # )

        ν_plot = collect(LinRange(0.0, 2π, 400)) 

        for (i, planet_classical) in enumerate(eachcol(planets_classical))
            temp = repeat(planet_classical, 1, 400)
            temp[6, :] .= ν_plot
            temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            lines!(ax1,
                temp[1, :],
                temp[2, :],
                temp[3, :],
                # linestyle = :dashdot,
                color = Makie.wong_colors()[[2, 1, 6][i]],
                # alpha = 0.35,
                label = ["Venus Orbit", "Earth Orbit", "Mars Orbit"][i]
            )
        end

        for n in [i]
            temp = deepcopy(p.id_journey[n])
            temp[1] = -3

            chosen_colors = cs[[findfirst(==(i), temp) for i in temp]]

            for k in 1:p.segment_number[n]
                x0 = p.x0[n][k]
                xf = p.xf[n][k]

                t_fine = collect(p.t_nodes[n][k][1]:1.0*day_scale:p.t_nodes[n][k][end])

                if !(t_fine[end] ≈ p.t_nodes[n][k][end])
                    push!(t_fine, p.t_nodes[n][k][end])
                end

                x_fine = integrate_trajectory(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    t_fine;
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                lines!(ax1,
                    x_fine[1, :],
                    x_fine[2, :],
                    x_fine[3, :],
                    color = :black,
                    alpha = 0.8,
                )
                
                scatter!(ax1,
                    x0[1],
                    x0[2],
                    x0[3],
                    color = chosen_colors[k],
                    markersize = 10,
                )

                if k == p.segment_number[n]
                    scatter!(ax1,
                        xf[1],
                        xf[2],
                        xf[3],
                        color = chosen_colors[k+1],
                        markersize = 10,
                    )
                end

                # lines!(ax1,
                #     x_target_fine[1, :],
                #     x_target_fine[2, :],
                #     x_target_fine[3, :],
                #     alpha = 0.2,
                #     color = :black
                # )

                thrust_force = if p.objective_config == LoggedMassConfig()
                    temp = exp.(p.x_nodes[n][k][7, 1:end-1])
                    p.u_nodes[n][k][4, :] .* temp * thrust * m_scale * a_scale * 1e3
                else
                    p.u_nodes[n][k][4, :] * thrust * m_scale * a_scale * 1e3
                end

                sel = thrust_force .>= 0.01

                thrust_vectors = stack(normalize.(eachcol(p.u_nodes[n][k][1:3, :]))).*thrust_force'

                length_scale = 0.2

                for i in collect(1:length(thrust_force))[sel]
                    lines!(ax1,
                        [p.x_nodes[n][k][1, i], p.x_nodes[n][k][1, i] + length_scale*thrust_vectors[1, i]],
                        [p.x_nodes[n][k][2, i], p.x_nodes[n][k][2, i] + length_scale*thrust_vectors[2, i]],
                        [p.x_nodes[n][k][3, i], p.x_nodes[n][k][3, i] + length_scale*thrust_vectors[3, i]],
                        color = :black,
                        alpha = 0.5,
                        linewidth = 1,
                    )
                end
            end

                
            deployments = sum(p.Δm0[n] .≈ -40/m_scale)

            m_fuel = if typeof(p.objective_config) == LoggedMassConfig
                m_scale*exp(p.x_nodes[n][1][7, 1]) - 500.0 - 40.0*deployments - (m_scale*exp(p.x_nodes[n][end][7, end]) - 500.0 + m_scale*p.Δm0[n][end])
            else
                m_scale*p.x_nodes[n][1][7, 1] - 500.0 - 40.0*deployments - (m_scale*p.x_nodes[n][end][7, end] - 500.0 + m_scale*p.Δm0[n][end])
            end

            m_returned = -m_scale*p.Δm0[n][end]

            t_nodes_combined = vcat(
                [p.t_nodes[n][k][1:end-1] .+ p.times_journey[n][k] for k in 1:length(p.x0[n])]...
            )
            u_nodes_combined = hcat(p.u_nodes[n]...)
            x_nodes_combined = hcat([val[:, 1:end-1] for val in p.x_nodes[n]]...)

            push!(t_nodes_combined, p.times_journey[end][end])
            u_nodes_combined = hcat(u_nodes_combined, p.u_nodes[end][end][:, end])
            x_nodes_combined = hcat(x_nodes_combined, p.x_nodes[end][end][:, end])

            thrust_force = if p.objective_config == LoggedMassConfig()
                temp = exp.(x_nodes_combined[7, :])
                u_nodes_combined[4, :] .* temp * thrust * m_scale * a_scale * 1e3
            else
                u_nodes_combined[4, :] * thrust * m_scale * a_scale * 1e3
            end
        end
    end


    # axislegend(ax1, merge = true, unique = true, orientation = :horizontal)

    # linkxaxes!(axs...)
    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end


function plot_trajectory_showcase(
    p::SequentialConvexProblem;
    output_file = nothing,
    solution_indices = nothing,
)
    f = Figure(size = (800, 400), backgroundcolor = :white, figure_padding = 2)

    if isnothing(solution_indices)
        solution_indices = collect(1:p.mixing_number)
    end
    
    cs = ColorSchemes.tab20
    # cs = ColorSchemes.tableau_miller_stone

    for i in [1]
    # for i in 1:p.mixing_number
        ax1 = Axis3(
            f[1, i]; 
            # xlabel = "x [AU]", 
            # ylabel = i==1 ? "y [AU]" : "", 
            limits = (-3.2, 3.2, -3.2, 3.2, -0.5, 0.5),
            elevation = 0.4,
            aspect = (4, 4, 1),
            # viewmode = :stretch,
            xspinesvisible = false,
            yspinesvisible = false,
            zspinesvisible = false,
            # limits = (-0.25, 3.2, -1.1, 1.1),
            # xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
            # yticks = [-1.0, -0.5, 0.0, 0.5, 1.0],
            # yticklabelsvisible = i==1,
            # limits = (0.5, 1.5, -3.0, -2.0),
            # aspect = 3.7/3
            # aspect = 1.0
        )

        ν_plot = collect(LinRange(0.0, 2π, 400)) 

        for (i, planet_classical) in enumerate(eachcol(planets_classical))
            temp = repeat(planet_classical, 1, 400)
            temp[6, :] .= ν_plot
            temp = hcat(classical_to_cartesian.(eachcol(temp))...)

            lines!(ax1,
                temp[1, :],
                temp[2, :],
                temp[3, :],
                # linestyle = :dashdot,
                color = Makie.wong_colors()[[2, 1, 6][i]],
                # alpha = 0.35,
                # label = ["Venus Orbit", "Earth Orbit", "Mars Orbit"][i]
            )
        end

        for n in 1:p.mixing_number
            temp = deepcopy(p.id_journey[n])
            temp[1] = -3

            # display(temp)

            chosen_colors = cs[([findfirst(==(i), temp) for i in temp] .- 1) .% 20 .+ 1]

            for k in 1:p.segment_number[n]
                x0 = p.x0[n][k]
                xf = p.xf[n][k]

                t_fine = collect(p.t_nodes[n][k][1]:1.0*day_scale:p.t_nodes[n][k][end])

                if !(t_fine[end] ≈ p.t_nodes[n][k][end])
                    push!(t_fine, p.t_nodes[n][k][end])
                end

                x_fine = integrate_trajectory(
                    p.x0[n][k] .+ vcat([0.0, 0.0, 0.0], p.Δv0[n][k], [0.0]),
                    t_fine;
                    t_nodes = p.t_nodes[n][k],
                    u_nodes = p.u_nodes[n][k],
                    p.objective_config,
                )

                lines!(ax1,
                    x_fine[1, :],
                    x_fine[2, :],
                    x_fine[3, :],
                    color = :black,
                    alpha = 0.2,
                )
                
                scatter!(ax1,
                    x0[1],
                    x0[2],
                    x0[3],
                    color = chosen_colors[k],
                    markersize = 10,
                )

                if k == p.segment_number[n]
                    scatter!(ax1,
                        xf[1],
                        xf[2],
                        xf[3],
                        color = chosen_colors[k+1],
                        markersize = 10,
                    )
                end

                # lines!(ax1,
                #     x_target_fine[1, :],
                #     x_target_fine[2, :],
                #     x_target_fine[3, :],
                #     alpha = 0.2,
                #     color = :black
                # )

                # thrust_force = if p.objective_config == LoggedMassConfig()
                #     temp = exp.(p.x_nodes[n][k][7, 1:end-1])
                #     p.u_nodes[n][k][4, :] .* temp * thrust * m_scale * a_scale * 1e3
                # else
                #     p.u_nodes[n][k][4, :] * thrust * m_scale * a_scale * 1e3
                # end

                # sel = thrust_force .>= 0.01

                # thrust_vectors = stack(normalize.(eachcol(p.u_nodes[n][k][1:3, :]))).*thrust_force'

                # length_scale = 0.2

                # for i in collect(1:length(thrust_force))[sel]
                #     lines!(ax1,
                #         [p.x_nodes[n][k][1, i], p.x_nodes[n][k][1, i] + length_scale*thrust_vectors[1, i]],
                #         [p.x_nodes[n][k][2, i], p.x_nodes[n][k][2, i] + length_scale*thrust_vectors[2, i]],
                #         [p.x_nodes[n][k][3, i], p.x_nodes[n][k][3, i] + length_scale*thrust_vectors[3, i]],
                #         color = :black,
                #         alpha = 0.5,
                #         linewidth = 1,
                #     )
                # end
            end
        end

        hidedecorations!(ax1)
    end




    # Label(
    #     f[1:4, 3], 
    #     label_text, 
    #     # font = :italic, 
    #     width = 150,
    #     justification = :left,
    #     padding = (-50, 0, 0, 0)
    # )


    # axislegend(ax1, merge = true, unique = true, orientation = :horizontal)

    

    # linkxaxes!(axs...)
    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end


function plot_solution_ships(
    average_mined_mass_all,
    allowed_mined_mass_all;
    max_ships = 30,
    max_mass = 800.0,
    output_file = nothing,
)

    f = Figure(size = (800, 350), backgroundcolor = :white)

    ax = Axis(f[1, 1],
        # title = "Solution Ships",
        xlabel = "number of ships in campaign",
        ylabel = "average ship mined mass [kg]",
        limits = ((0.5, max_ships + 0.5), (600.0, max_mass)),
        xticks = 0:1:max_ships,
        xticklabelsize = 10,
        # yticklabelsize = 10,
        # leftpadding = 20,
        # bottompadding = 20,
        # xgridvisible = false
    )

    cs = ColorSchemes.tab10

    val1 = 0.3
    val2 = 0.2

    for i in 1:max_ships
        lines!(ax, [i-val1, i+val1], [average_mined_mass_all[i], average_mined_mass_all[i]], color=:black, linewidth=2, alpha=0.8, label="average mined mass")

        # rangebars!(ax, 
        #     [i], 
        #     minimum(mined_mass_all[chosen_ships_all[i]]), 
        #     maximum(mined_mass_all[chosen_ships_all[i]]), 
        #     color = :black,
        #     whiskerwidth = 15, 
        #     direction = :y
        # )

        lines!(ax, [i-val1, i+val1], [allowed_mined_mass_all[i], allowed_mined_mass_all[i]], color=cs[4], linewidth=2, alpha=0.8, label="permitted mined mass")

        for j in chosen_ships_all[i]
            temp = if length(findall(==(groups_all[j]), groups_all)) == 1
                lines!(ax, [i-val2, i+val2], [mined_mass_all[j], mined_mass_all[j]], color=cs[1], alpha=0.4, linewidth=2, label="self-cleaning ship")
            else
                lines!(ax, [i-val2, i+val2], [mined_mass_all[j], mined_mass_all[j]], color=cs[2], alpha=0.4, linewidth=2, label="mixed ship")
            end

        end
    end

    axislegend(ax, unique = true, merge = true, orientation = :vertical, position = :lb)


    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end


function plot_solution_ships(
    average_mined_mass_all,
    allowed_mined_mass_all;
    max_ships = 30,
    max_mass = 800.0,
    output_file = nothing,
)

    f = Figure(size = (800, 350), backgroundcolor = :white)

    ax = Axis(f[1, 1],
        # title = "Solution Ships",
        xlabel = "number of ships in campaign",
        ylabel = "average ship mined mass [kg]",
        limits = ((0.5, max_ships + 0.5), (600.0, max_mass)),
        xticks = 0:1:max_ships,
        xticklabelsize = 10,
        # yticklabelsize = 10,
        # leftpadding = 20,
        # bottompadding = 20,
        # xgridvisible = false
    )

    cs = ColorSchemes.tab10

    val1 = 0.3
    val2 = 0.2

    for i in 1:max_ships
        lines!(ax, [i-val1, i+val1], [average_mined_mass_all[i], average_mined_mass_all[i]], color=:black, linewidth=2, alpha=0.8, label="average mined mass")

        # rangebars!(ax, 
        #     [i], 
        #     minimum(mined_mass_all[chosen_ships_all[i]]), 
        #     maximum(mined_mass_all[chosen_ships_all[i]]), 
        #     color = :black,
        #     whiskerwidth = 15, 
        #     direction = :y
        # )

        lines!(ax, [i-val1, i+val1], [allowed_mined_mass_all[i], allowed_mined_mass_all[i]], color=cs[4], linewidth=2, alpha=0.8, label="permitted mined mass")

        for j in chosen_ships_all[i]
            temp = if length(findall(==(groups_all[j]), groups_all)) == 1
                lines!(ax, [i-val2, i+val2], [mined_mass_all[j], mined_mass_all[j]], color=cs[1], alpha=0.4, linewidth=2, label="self-cleaning ship")
            else
                lines!(ax, [i-val2, i+val2], [mined_mass_all[j], mined_mass_all[j]], color=cs[2], alpha=0.4, linewidth=2, label="mixed ship")
            end

        end
    end

    axislegend(ax, unique = true, merge = true, orientation = :vertical, position = :lb)


    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end



function plot_bip_solution_values(
    mip_problem1, 
    mip_problem2;
    output_file = nothing
)

    f = Figure(size = (800, 240), backgroundcolor = :white)

    cs = ColorSchemes.tableau_10
    
    ax = Axis(
        f[1, 1]; 
        xlabel = "total transfer Δv [km/s]", 
        ylabel = "count", 
        # xticks = 0:5:50,
        # limits = (0, 51, 720, 790)
        # xticks = 0:10:50
    )

    thresholded_solutions = [val in mip_problem1.id_journey_solutions for val in mip_problem2.id_journey_solutions]

    solution_indices = collect(1:mip_problem2.solutions)

    bins = collect(25.0:0.5:60.0)
    # bins = collect(10.0:0.5:40.0)


    bin_centers = 0.5*(bins[1:end-1] + bins[2:end])

    bin_counts_thresholded = [sum((mip_problem1.objective_solutions.*v_scale .>= bins[i]) .& (mip_problem1.objective_solutions.*v_scale .< bins[i+1])) for i in 1:length(bins)-1]
    
    bin_counts_not_thresholded = [sum((mip_problem2.objective_solutions[.!thresholded_solutions].*v_scale .>= bins[i]) .& (mip_problem2.objective_solutions[.!thresholded_solutions].*v_scale .< bins[i+1])) for i in 1:length(bins)-1]

    for i in 1:length(bin_centers)
        display(bin_counts_thresholded[i])

        scatter!(
            ax,
            fill(bin_centers[i], bin_counts_thresholded[i]),
            collect(1:bin_counts_thresholded[i]),
            color = cs[1],
            label = "remaining",
            markersize = 8,
        )

        scatter!(
            ax,
            fill(bin_centers[i], bin_counts_not_thresholded[i]),
            bin_counts_thresholded[i] .+ collect(1:bin_counts_not_thresholded[i]),
            color = cs[2],
            label = "pruned",
            markersize = 8,
        )
    end



    # scatter!(
    #     ax,
    #     solution_indices[.!thresholded_solutions],
    #     mip_problem2.objective_solutions[.!thresholded_solutions].*v_scale,
    #     color = cs[2],
    #     label = "pruned",
    #     markersize = 8,
    # )

    # scatter!(
    #     ax,
    #     solution_indices[thresholded_solutions],
    #     mip_problem2.objective_solutions[thresholded_solutions].*v_scale,
    #     color = cs[1],
    #     label = "remaining",
    #     markersize = 8,
    # )



    axislegend(ax, unique = true, merge = true, orientation = :vertical, position = :rt)

    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return
end