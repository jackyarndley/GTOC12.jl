


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
    solution_index = 1
)
    n = solution_index

    f = Figure(size = (900, 800), backgroundcolor = :white)

    x_nodes_combined = hcat([val[:, 1:end-1] for val in p.x_nodes[n]]...)

    ax1 = Axis(
        f[1:2, 1:2]; 
        xlabel = "x [AU]", 
        ylabel = "y [AU]", 
        aspect = 1
    )

    for k in 1:length(p.x0[n])
        scatter!(ax1,
            p.x0[n][k][1],
            p.x0[n][k][2],
        )

        if k == length(p.x0[n])
            scatter!(ax1,
                p.xf[n][k][1],
                p.xf[n][k][2],
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
        )



        thrust_force = if p.objective_config == LoggedMassConfig()
            temp = exp.(p.x_nodes[n][k][7, 1:end-1])
            p.u_nodes[n][k][4, :] .* temp * thrust * m_scale * a_scale * 1e3
        else
            p.u_nodes[n][k][4, :] * thrust * m_scale * a_scale * 1e3
        end

        sel = thrust_force .>= 0.01

        thrust_vectors = stack(normalize.(eachcol(p.u_nodes[n][k][1:3, :]))).*thrust_force'

        arrows!(ax1,
            p.x_nodes[n][k][1, 1:end-1][sel],
            p.x_nodes[n][k][2, 1:end-1][sel],
            thrust_vectors[1, sel],
            thrust_vectors[2, sel],
            lengthscale = 0.2,
            arrowhead = ' ',
            linecolor = (:black, 0.4),
        )

    end




    limits!(ax1, -3.0, 3.0, -3.0, 3.0)

    # hidedecorations!(
    #     ax1
    # )

    resize_to_layout!(f)
    display(f)

    return
end