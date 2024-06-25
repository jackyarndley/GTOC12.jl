


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


function plot_thrust_information(p::SequentialConvexProblem)
    mixing = length(p.x0)

    f = Figure(size = (800, 250*mixing), backgroundcolor = :white)

    axs = []

    for n in 1:mixing
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
            f[n, 1]; 
            xlabel = n == mixing ? "t [MJD]" : "", 
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
)
    f = Figure(size = (900, 800), backgroundcolor = :white)

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
            limits = (-3.0, 3.0, -3.0, 3.0),
            aspect = 1
        )
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

    for n in 1:p.mixing_number
        for k in 1:p.segment_number[n]

            scatter!(ax1,
                p.x0[n][k][1],
                p.x0[n][k][2],
                p.x0[n][k][3],
            )

            if k == p.segment_number[n]
                scatter!(ax1,
                    p.xf[n][k][1],
                    p.xf[n][k][2],
                    p.xf[n][k][3],
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
    color_range = ColorSchemes.tableau_jewel_bright
    # color_range = ColorSchemes.tableau_miller_stone

    for stage in 1:p.deployment_nums[1]
        for i in 1:stage_size
            push!(layout, ((stage - 1)*stage_Δx, -(i-1)*node_Δy))
        end
    end

    for stage in 1:p.deployment_arcs[1]
        for i in 1:stage_size, j in 1:stage_size
            if (p.deployment_cost[1][i, j, stage] <= p.cost_limit) && (i != j)
                used_check = false

                temp = []

                for k in p.solutions:-1:1
                # for k in 1:p.solutions
                    if p.id_journey_solutions[k][1][stage + 1] == p.id_subset[i] && p.id_journey_solutions[k][1][stage + 2] == p.id_subset[j]
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

                if !used_check
                    push!(unused_edge_list, (i + stage_size*(stage-1), j + stage_size*stage))
                end
            end
        end
    end

    collection_offset = length(layout)

    for i in 1:stage_size, j in 1:stage_size
        if (p.intermediate_cost[1][i, j] <= p.cost_limit)
            used_check = false

            temp = []

            for k in p.solutions:-1:1
            # for k in 1:p.solutions
                if p.id_journey_solutions[k][1][p.deployment_nums[1] + 1] == p.id_subset[i] && p.id_journey_solutions[k][1][p.deployment_nums[1] + 2] == p.id_subset[j]
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

            if !used_check
                push!(unused_edge_list, (collection_offset - stage_size + i, collection_offset + j))
            end
        end
    end

    for stage in (p.deployment_nums[1] + 1):(p.deployment_nums[1] + p.collection_nums[1])
        for i in 1:stage_size
            push!(layout, ((stage - 1 + 1)*stage_Δx, -(i-1)*node_Δy))
        end
    end

    for stage in 1:p.collection_arcs[1]
        for i in 1:stage_size, j in 1:stage_size
            if (p.collection_cost[1][i, j, stage] <= p.cost_limit) && (i != j)
                used_check = false

                temp = []

                for k in p.solutions:-1:1
                # for k in 1:p.solutions
                    if p.id_journey_solutions[k][1][p.deployment_nums[1] + 1 + stage] == p.id_subset[i] && p.id_journey_solutions[k][1][p.deployment_nums[1] + 2 + stage] == p.id_subset[j]
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

                if !used_check
                    push!(unused_edge_list, (collection_offset + i + stage_size*(stage-1), collection_offset + j + stage_size*stage))
                end
            end
        end
    end
    
    unused_graph = SimpleDiGraph(length(layout))

    [add_edge!(unused_graph, e) for e in unused_edge_list]

    fixed_layout(_) = layout

    f = Figure(size = (1000, 600), backgroundcolor = :white)

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
        ilabels = repeat(p.id_subset, p.deployment_nums[1] + p.collection_nums[1]),
        ilabels_fontsize = 9,
        node_size = 30,
        node_marker=Circle;
    )

    hidedecorations!(ax)
    hidespines!(ax)
    




    resize_to_layout!(f)
    display(f)

    # save("test.png", f)

    return
end