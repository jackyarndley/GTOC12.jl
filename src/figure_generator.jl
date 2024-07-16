include("header.jl")

GLMakie.activate!()
CairoMakie.activate!()

using LaTeXStrings











id_journey = [
    [0, 15184, 12286, 15184, 12286, -3]
]

times_journey = [
    convert_mjd_to_time([64428.0 64913.0 65208.0 65648.0 66145.0 66550.0][:]),
]




scp_iterations = 40

node_time_spacing = 20.0*day_scale





p = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.02,
    mass_overhead = 1.0/m_scale
);

solve!(p)









plot_trajectory(p)



plot_gtoc12_problem(p)


plot_scp_details(p)


plot_thrust_information(p)


mjd_end = 69807



id_subset = [3241, 15184, 19702, 46418, 53592]






mip_problem = MixedIntegerProblem(id_subset, [3], [3])
mip_problem.cost_limit = 6/v_scale

solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 3
)

times_journey = mip_problem.times_journey
id_journey = mip_problem.id_journey_solutions[1]



p1 = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.05,
    mass_overhead = 1.0/m_scale
);

solve!(p1; 
    fixed_segments = true,
    fixed_rendezvous = true
)

p2 = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.025,
    mass_overhead = 1.0/m_scale
);

solve!(p2; 
    fixed_segments = false,
    fixed_rendezvous = true
)

p3 = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.05,
    mass_overhead = 1.0/m_scale
);

solve!(p3; 
    fixed_segments = false,
    fixed_rendezvous = false
)




plot_discretization_comparison(p1, p2, p3; output_file = "output/plots/discretization_comparison.png")











function plot_gtoc12_problem(
    p::SequentialConvexProblem
)

    cs = ColorSchemes.tableau_10

    f = Figure(size = (1000, 480), backgroundcolor = :white)

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
                    "",
                    "Asteroid 1\nDeployment",
                    "Asteroid 2\nDeployment",
                    "Asteroid 1\nCollection",
                    "Asteroid 2\nCollection",
                    ""
                ]

                
                offsets = [
                    (5, 0),
                    (5, 0),
                    (5, 0),
                    (5, 0),
                    (-5, 0),
                    (-5, 0),
                ]

                alignment = [
                    (:left, :bottom),
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





function plot_graph_structure(
    p::MixedIntegerProblem;
    plot_optimal_path = false,
    plot_pruned = true,
    output_file = nothing,
    selection_index = 1,
    figure_size = (900, 200)
)

    stage_size = length(p.id_subset)

    stage_Δx = 0.1
    node_Δy = 0.05

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

    f = Figure(size = figure_size, backgroundcolor = :white, figure_padding = 0)

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
        ilabels = repeat(1:length(mip_problem.id_subset), p.deployment_nums[selection_index] + p.collection_nums[selection_index]),
        # ilabels_fontsize = 12,
        # node_size = 20,
        ilabels_fontsize = 8,
        node_size = 15,
        node_marker=Circle;
    )

    # time_positions = vcat(
    #     [stage_Δx*(stage-1+0.5) for stage in 1:p.deployment_arcs[selection_index]],
    #     [stage_Δx*(p.deployment_arcs[selection_index] + 1)],
    #     [stage_Δx*(stage-1+1+0.5) for stage in (p.deployment_nums[selection_index] + 1):(p.deployment_nums[selection_index] + p.collection_arcs[selection_index])]
    # )

    # time_strings = vcat(
    #     [latexstring("d_{i, j, $(v)}") for v in 1:p.deployment_arcs[selection_index]],
    #     [latexstring("m_{i, j}")],
    #     [latexstring("c_{i, j, $(v)}") for v in 1:p.collection_arcs[selection_index]],
    # )

    # # display(time_positions)
    # # display(time_strings)

    # text!(ax, 
    #     time_positions,
    #     fill(0.7*node_Δy, length(time_positions)),
    #     text = time_strings,
    #     align = (:center, :center),
    #     fontsize = 20,
    # )

    hidedecorations!(ax)
    hidespines!(ax)

    xlims!(ax, [-stage_Δx*0.2, stage_Δx*(p.deployment_nums[selection_index] + p.collection_arcs[selection_index] + 1) + stage_Δx*0.2])
    ylims!(ax, [-node_Δy*(length(p.id_subset) - 0.5), node_Δy*0.5])
    

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
    f = Figure(size = (800, 400), backgroundcolor = :white)
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
            xlabel = i == 3 ? "time " : "", 
            xticksvisible = false,
            xticklabelsvisible = false,
            ylabel = "thrust", 
            # xticks = [65000, 66000, 67000, 68000, 69000],
            yticks = ([0.0, 0.6], ["0.0", "max"]),
            xgridvisible = false,
            limits = (mjd_start, 65420, -0.05, 0.65)
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
    scp_problem_objectives
)

    f = Figure(size = (800, 400), backgroundcolor = :white)

    cs = ColorSchemes.tableau_10
    
    ax = Axis(
        f[1, 1]; 
        xlabel = "BIP problem index", 
        ylabel = "SCP problem objective kg", 
    )

    for i in 1:length(mip_problem_objectives)
        scatter!(ax,
            mip_problem_objectives[i],
            scp_problem_objectives[i],
            label = "Iteration $(i)"
        )
    end

    f[1, 2] = Legend(f, ax, framevisible = false, unique = true, merge = true, labelsize=16)

    resize_to_layout!(f)
    display(f)

    # if !isnothing(output_file)
    #     save(output_file, f)
    # end

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
# id_subset = sort([2032, 3241, 15184, 19702, 23056, 23987, 32088, 46418, 46751, 53592, 3896, 37818, 15083, 5707, 19434, 981, 48748, 40804, 23483, 47817, 2174, 28289, 43836, 39557, 9260, 17983, 13655, 22108, 3302, 57913])
# id_subset = sort([2032, 3241, 15184, 19702, 23056, 23987, 32088, 46418, 46751, 53592, 3896, 37818, 15083, 5707, 19434, 981, 48748, 40804, 23483, 47817])


mip_problem = MixedIntegerProblem(id_subset, [10], [10])
mip_problem.cost_limit = 10/v_scale


join([@sprintf("%5s ", val) for val in id_subset])

join([@sprintf("%5i ", val) for val in convert_time_to_mjd.(mip_problem.times_journey[1])])

join([@sprintf("%5s ", val) for val in mip_problem.id_journey_solutions[5][1]])




mip_problem_objectives = Vector{Float64}[]
scp_problem_objectives = Vector{Float64}[]





solution_number = 20


solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.2,
    solutions_count_maximum = 3*solution_number,
    time_limit_seconds = 300
)


push!(mip_problem_objectives, mip_problem.objective_solutions.*v_scale)

scp_problem = SequentialConvexProblem(
    [mip_problem.id_journey_solutions[k][1] for k in 1:min(solution_number, mip_problem.solutions)], 
    [mip_problem.times_journey[1] for k in 1:min(solution_number, mip_problem.solutions)];
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.05,
    mass_overhead = 0.0/m_scale,
    dynamical_error = 1e-4
);

solve!(scp_problem)



push!(scp_problem_objectives, [-m_scale*scp_problem.Δm0[n][end] for n in 1:scp_problem.mixing_number])




plot_convergence_comparison(
    [collect(1:min(solution_number, length(mip_problem_objectives[i]))) for i in 1:length(mip_problem_objectives)],
    scp_problem_objectives
)

# plot_convergence_comparison(
#     mip_problem_objectives,
#     scp_problem_objectives
# )




mip_problem = MixedIntegerProblem(id_subset, [10], [10];
    times_journey = [scp_problem.times_journey[argmax(scp_problem_objectives[end])]]
)

mip_problem.cost_limit = 10/v_scale






plot_thrust_information(scp_problem; solution_indices = [1])

plot_trajectory(scp_problem; solution_indices = [1], plot_3d = false, rotating = true)






plot_graph_structure(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_big.png",
    figure_size = (900, 500)
)







# times_5day = collect(0:5*day_scale:maximum_time)
# asteroids_cartesian_5day = ephemeris_cartesian_at(asteroids_classical, times_5day)



# check_times_deploy = collect(convert_mjd_to_5day_time(64900):20:convert_mjd_to_5day_time(66500))
# check_times_collect = collect(convert_mjd_to_5day_time(68200):20:convert_mjd_to_5day_time(69000))


# times_join_check = collect(convert_mjd_to_time(65000):100*day_scale:convert_mjd_to_time(69000))
# asteroids_join_check = ephemeris_cartesian_at(asteroids_classical, times_join_check)









while length(id_subset) < 20
    print("Size: \n$(length(id_subset)), last ID: $(id_subset[end])")

    temp = []

    for id in id_subset
        temp10 = all_ids[[minimum(norm.(eachcol(a .- asteroids_join_check[1:3, id, :]))) for a in eachslice(asteroids_join_check[1:3, :, :], dims=2)] .<= 0.05]

        temp = unique(vcat(temp, temp10))
    end

    temp = setdiff(temp, id_subset)
    # temp = setdiff(temp, used_ids)

    print(", near asteroids: $(length(temp))")

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







# id_journey, times_journey, mined_mass_submitted, penalised_mass_submitted, groups_submitted, result_files_submitted = load_result_folders_grouping(["data/bundled/submitted/OptimiCS.txt"])




# p = SequentialConvexProblem(
#     id_journey, 
#     times_journey;
#     objective_config = LoggedMassConfig(),
#     trust_region_factor = 0.025,
#     mass_overhead = 0.0/m_scale,
#     dynamical_error = 1e-4
# );

# solve!(p)







plot_team_improvements(
    [
        "data/bundled/submitted/JPL.txt",
        "data/bundled/submitted/BIT.txt",
        "data/bundled/submitted/OptimiCS.txt",
        "data/bundled/submitted/ESA.txt",
        "data/bundled/submitted/TheAntipodes.txt",
        "data/bundled/submitted/NUDT.txt",
        "data/bundled/submitted/SIGMA.txt",
        "data/bundled/submitted/ATQ.txt",
        "data/bundled/submitted/ADL.txt",
        "data/bundled/submitted/NUAA.txt",
    ],
    [
        "data/bundled/reoptimized/JPL"
        "data/bundled/reoptimized/BIT"
        "data/bundled/reoptimized/OptimiCS"
        "data/bundled/reoptimized/ESA"
        "data/bundled/reoptimized/TheAntipodes"
        "data/bundled/reoptimized/NUDT"
        "data/bundled/reoptimized/SIGMA"
        "data/bundled/reoptimized/ATQ"
        "data/bundled/reoptimized/ADL"
        "data/bundled/reoptimized/NUAA"

    ],
    [
        "JPL",
        "BIT-CAS-DFH",
        "OptimiCS",
        "ESA ACT\n& Friends",
        "TheAntipodes",
        "NUDT-LIPSAM",
        "∑ TEAM",
        "ATQ",
        "ADL",
        "NUAA-ASTL",
    ];
    output_file = "output/plots/top10_convex.png"
)











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
        xgridvisible = false,
        ygridvisible = false,
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
            label=false
        )

        for (mass, color) in zip(improved_mass, colors)
            lines!(ax,
                [mass, mass], 
                [y_position-0.18, y_position+0.18],
                color=cs[color],
                alpha=0.75,
                linewidth=2,
                label=false
            )
        end

        push!(axs, ax)
    end

    linkxaxes!(axs...)
    resize_to_layout!(f)
    display(f)

    if !isnothing(output_file)
        save(output_file, f)
    end

    return


end


