include("header.jl")

GLMakie.activate!()
CairoMakie.activate!()

using LaTeXStrings











id_journey = [
    [0, 15184, 12286, 15184, 12286, -3]
]

times_journey = [
    convert_mjd_to_time([64428.0 65013.0 65258.0 65778.0 66013.0 66583.0][:]),
]




scp_iterations = 80

node_time_spacing = 20.0*day_scale





p = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.05,
    mass_overhead = 1.0/m_scale
);

solve!(p,
    MixedTimeAdaptive(); 
    adaptive_time = true
)




plot_gtoc12_problem(p)





function plot_gtoc12_problem(
    p::SequentialConvexProblem
)

    cs = ColorSchemes.tableau_10

    f = Figure(size = (1200, 600), backgroundcolor = :white)

    ax1 = Axis3(
        f[1, 1]; 
        # xlabel = "x [AU]", 
        # ylabel = "y [AU]", 
        limits = 0.84.*(-3.0, 3.0, -3.0, 3.0, -1.0, 1.0),
        azimuth = -π/2,
        elevation = π/10,
        aspect = (3, 3, 1)
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
                linestyle = :dashdot,
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
                linestyle = :dashdot,
                linewidth = 2,
                alpha = 0.5,
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
                    alpha = 1.0,
                    linestyle = :solid,
                    label = "Mining Ship\nTrajectory"
                )

                scatter!(ax,
                    p.x0[n][k][1],
                    p.x0[n][k][2],
                    p.x0[n][k][3],
                    color = :black
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
                    "Earth Departure",
                    "Asteroid 1\nDeployment",
                    "Asteroid 2\nDeployment",
                    "Asteroid 1\nCollection",
                    "Asteroid 2\nCollection",
                    "Earth Arrival"
                ]

                alignment = [
                    (:left, :bottom),
                    (:left, :bottom),
                    (:left, :bottom),
                    (:left, :bottom),
                    (:left, :top),
                    (:right, :bottom)
                ]
                
                text!(ax, 
                    p.x0[n][k][1],
                    p.x0[n][k][2],
                    p.x0[n][k][3],
                    text = text[k],
                    align = alignment[k],
                    # fontsize = 12,
                    offset = (5, 0)
                )

                if k == p.segment_number[n]
                    text!(ax, 
                        p.xf[n][k][1],
                        p.xf[n][k][2],
                        p.xf[n][k][3],
                        text = text[end],
                        align = alignment[end],
                        # fontsize = 12,
                        offset = (-5, 0),
                    )
                end

            end
        end


        hidedecorations!(ax)
        hidespines!(ax)

        axislegend(ax, unique = true, merge = true)

    end

    resize_to_layout!(f)
    display(f)

    save("output/plots/gtoc12_explainer.png", f)

    return
end



function plot_graph_structure(
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
        ilabels = repeat(1:5, p.deployment_nums[selection_index] + p.collection_nums[selection_index]),
        ilabels_fontsize = 9,
        node_size = 30,
        node_marker=Circle;
    )

    time_positions = vcat(
        [stage_Δx*(stage-1+0.5) for stage in 1:p.deployment_arcs[selection_index]],
        [stage_Δx*(p.deployment_arcs[selection_index] + 1)],
        [stage_Δx*(stage-1+1+0.5) for stage in (p.deployment_nums[selection_index] + 1):(p.deployment_nums[selection_index] + p.collection_arcs[selection_index])]
    )

    time_strings = vcat(
        [latexstring("d_{i, j, $(v)}") for v in 1:p.deployment_arcs[selection_index]],
        [latexstring("m_{i, j}")],
        [latexstring("c_{i, j, $(v)}") for v in 1:p.collection_arcs[selection_index]],
    )

    display(time_positions)
    display(time_strings)
    

    text!(ax, 
        time_positions,
        fill(0.7*node_Δy, length(time_positions)),
        text = time_strings,
        align = (:center, :center),
        fontsize = 20,
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




id_subset = [3241, 15184, 19702, 46418, 53592]


mip_problem = MixedIntegerProblem(id_subset, [3], [3])
mip_problem.cost_limit = 6/v_scale



solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 3
)

plot_graph_structure(
    mip_problem;
    plot_pruned = false,
    plot_optimal_path = false,
    output_file = "output/plots/bip_connected.png"
)


plot_graph_structure(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = false,
    output_file = "output/plots/bip_pruned.png"
)

plot_graph_structure(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_solutions.png"
)





id_subset = [3241, 15184, 19702, 46418, 53592]


mip_problem = MixedIntegerProblem(id_subset, [3, 2], [2, 3])
mip_problem.cost_limit = 6/v_scale



solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 3
)


plot_graph_structure(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_mixed_1.png",
    selection_index = 1
)

plot_graph_structure(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_mixed_2.png",
    selection_index = 2
)








id_subset = sort([15184, 3241, 2032, 53592, 46418, 19702, 23056, 46751, 32088, 23987])






mip_problem = MixedIntegerProblem(id_subset, [10], [10])
mip_problem.cost_limit = 8/v_scale




solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 40,
    time_limit_seconds = 300
)




scp_problem = SequentialConvexProblem(
    [mip_problem.id_journey_solutions[k][1] for k in 1:20], 
    [mip_problem.times_journey[1] for k in 1:20];
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.025,
    mass_overhead = 1.0/m_scale
);


solve!(scp_problem,
    MixedTimeAdaptive(); 
    adaptive_time = true
)







mip_problem = MixedIntegerProblem(id_subset, [10], [10];
    times_journey = [scp_problem.times_journey[20]]
)

mip_problem.cost_limit = 8/v_scale




solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 80
)













times_5day = collect(0:5*day_scale:maximum_time)
asteroids_cartesian_5day = ephemeris_cartesian_at(asteroids_classical, times_5day)



check_times_deploy = collect(convert_mjd_to_5day_time(64900):20:convert_mjd_to_5day_time(66500))
check_times_collect = collect(convert_mjd_to_5day_time(68200):20:convert_mjd_to_5day_time(69000))


times_join_check = collect(convert_mjd_to_time(65000):100*day_scale:convert_mjd_to_time(69000))
asteroids_join_check = ephemeris_cartesian_at(asteroids_classical, times_join_check)









while length(id_subset) < 50
    print("\n$(length(id_subset)) ID $(id_subset[end])")

    temp = []

    for id in id_subset
        temp10 = all_ids[[minimum(norm.(eachcol(a .- asteroids_join_check[1:3, id, :]))) for a in eachslice(asteroids_join_check[1:3, :, :], dims=2)] .<= 0.1]

        temp = unique(vcat(temp, temp10))
    end

    temp = setdiff(temp, id_subset)
    # temp = setdiff(temp, used_ids)

    print(", near asteroids = $(length(temp))")

    average_minimum_dv = fill(0.0, length(temp))
    # average_minimum_dv = fill(1000.0, length(temp))

    for id in id_subset
        transfer_dvs_deploy = [v_scale.*get_lambert_phase(id, i, check_times_deploy) for i in temp]
        transfer_dvs_collect = [v_scale.*get_lambert_phase(id, i, check_times_collect) for i in temp]

        average_minimum_dv .+= sum.(transfer_dvs_deploy)./length(check_times_deploy) .+ sum.(transfer_dvs_collect)./length(check_times_collect)
        # average_minimum_dv .+= minimum.(transfer_dvs_deploy) .+ minimum.(transfer_dvs_collect)
        # average_minimum_dv = min.(average_minimum_dv, minimum.(transfer_dvs_deploy) .+ minimum.(transfer_dvs_collect))
    end

    minimum_id = temp[argmin(average_minimum_dv)]

    print(", lowest Δv = $(minimum(average_minimum_dv))")

    push!(id_subset, minimum_id)
end









