


function get_event_label_and_color(
    id,
    mass_change
)
    label = if id == 0
        "Earth"
    elseif id < 0
        ["Venus", "Earth", "Mars"][-id - 1]
    else
        string(id)
    end

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
    id_journey,
    times_journey,
    mass_changes,
    t_nodes,
    u_nodes,
    x0,
    objective_config
)
    x0 = copy(x0)
    t_nodes_mjd = convert_time_to_mjd(t_nodes)
    times_journey_mjd = convert_time_to_mjd(times_journey)

    if objective_config == LoggedMassConfig()
        x0[7] = log(x0[7])
        extra_callback = get_log_mass_change_callback_dynamic(times_journey .- times_journey[1], mass_changes)
    else
        extra_callback = get_mass_change_callback_dynamic(times_journey .- times_journey[1], mass_changes)
    end

    x_nodes = propagate_spacecraft(
        x0, 
        t_nodes .+ 1e-10; 
        t_nodes = t_nodes, 
        u_nodes = u_nodes, 
        thrust_config = ZeroOrderHoldConfig(),
        objective_config,
        extra_callback
    )

    thrust_force = if objective_config == LoggedMassConfig()
        temp = exp.(x_nodes[7, :])
        u_nodes[4, :] .* temp * thrust * m_scale * a_scale * 1e3
    else
        u_nodes[4, :] * thrust * m_scale * a_scale * 1e3
    end

    f = Figure(size = (800, 300), backgroundcolor = :white)

    ax1 = Axis(
        f[1, 1]; 
        xlabel = "t [MJD]", 
        ylabel = "thrust [N]", 
    )

    ylims!(ax1, (-0.05 * maximum(thrust_force), 1.2 * maximum(thrust_force)))

    plot_location_vertical_markers(ax1, times_journey_mjd, id_journey, mass_changes; vertical_position = 1.05*maximum(thrust_force))

    stairs!(
        ax1, 
        t_nodes_mjd, 
        thrust_force;
        step = :post,
        color = :black,
        linewidth = 1.5,
    )

    resize_to_layout!(f)
    display(f)

    return
end
