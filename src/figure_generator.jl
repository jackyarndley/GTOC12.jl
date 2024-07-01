include("header.jl")

GLMakie.activate!()
CairoMakie.activate!()












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





id_subset = [3241, 15184, 19702, 46418, 53592]


mip_problem = MixedIntegerProblem(id_subset, [3], [3])
mip_problem.cost_limit = 8/v_scale



solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 3
)


plot_graph(
    mip_problem;
    plot_pruned = false,
    plot_optimal_path = false,
    output_file = "output/plots/bip_connected.png"
)

plot_graph(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = false,
    output_file = "output/plots/bip_pruned.png"
)

plot_graph(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_solutions.png"
)




id_subset = [3241, 15184, 19702, 46418, 53592]


mip_problem = MixedIntegerProblem(id_subset, [3, 2], [2, 3])
mip_problem.cost_limit = 8/v_scale



solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 1
)



plot_graph(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_mixed_1.png",
    selection_index = 1,
)


plot_graph(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_mixed_2.png",
    selection_index = 2,
)







id_subset = [298 3889 4445 7103 10290 10649 10916 12577 14062 14079 15291 15906 16070 16110 18835 18913 19496 20616 20651 22313 24024 25663 27293 27414 28050 30492 31751 32796 34679 36517 36565 36716 38270 38514 41013 41509 41629 43596 46392 46789 46877 47175 47291 47647 49082 49251 50836 50898 53789 54572 55377 56717 57669 57998 59572 59916][:][1:2:end]



mip_problem = MixedIntegerProblem(id_subset, [7], [7])
mip_problem.cost_limit = 6/v_scale



solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.2,
    solutions_count_maximum = 3
)



plot_graph(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_big.png",
    figure_size = (1200, 1200)
)