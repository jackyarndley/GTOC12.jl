include("header.jl")

id_journey = [0, 56608, 56608, -3]
times_journey = [1.0321259370000504, 10.827001079130554, 85.01277301090416, 92.92573852790454]




mass_changes = get_mass_change_at_ids(id_journey, times_journey)
# mass_changes = get_mass_change_at_ids_mixed(
#     id_journey,
#     times_journey,
#     id_journey_mixed,
#     times_journey_mixed
# )

locations_journey, Δv_journey_combined, Δv_journey_departure, Δv_journey_arrival = get_lambert_trajectory(
    times_journey,
    id_journey
)


# Plot the test trajectory
# plot_ship_trajectory_lambert(
#     times_journey, 
#     id_journey, 
#     Δv_journey_combined,
#     planets_classical
# )


node_time_spacing = 5.0 * day_scale

mass_starting = 3000.0 / m_scale
reverse_mass = false
segments = length(id_journey) - 1

x_departure_cartesian, x_arrival_cartesian, t_nodes, u_nodes, Δv_departure_injection, Δv_departure_injection_limit, Δv_arrival_injection, Δv_arrival_injection_limit = get_lambert_guess_for_scp(
    Δv_journey_departure,
    Δv_journey_arrival,
    locations_journey,
    id_journey,
    times_journey,
    mass_starting,
    mass_changes;
    reverse_mass
);

x_nodes, u_nodes, Δv_departure_injection, Δv_arrival_injection, termination, maximum_error = solve_scp_full(
    x_departure_cartesian,
    x_arrival_cartesian,
    t_nodes,
    u_nodes,
    Δv_departure_injection,
    Δv_departure_injection_limit,
    Δv_arrival_injection,
    Δv_arrival_injection_limit,
    mass_changes,
    id_journey,
    LoggedMassConfig();
    display_output = true,
    linearization_error = 1e-4,
    reverse_mass,
);

# Refine as guess for normal mass problem
u_nodes, x_departure_cartesian = convert_log_mass_control_to_mass_control(locations_journey, mass_starting, t_nodes, u_nodes, Δv_departure_injection, mass_changes, x_departure_cartesian);

x_nodes, u_nodes, Δv_departure_injection, Δv_arrival_injection, termination, maximum_error = solve_scp_full(
    x_departure_cartesian,
    x_arrival_cartesian,
    t_nodes,
    u_nodes,
    Δv_departure_injection,
    Δv_departure_injection_limit,
    Δv_arrival_injection,
    Δv_arrival_injection_limit,
    mass_changes,
    id_journey,
    MassConfig();
    display_output = true,
    linearization_error = 1e-7,
    reverse_mass,
);


mass_starting = min(1.0, x_nodes[1][7, 1])

ga_times, ga_amount = get_gravity_assist_information(times_journey, id_journey, Δv_departure_injection, Δv_arrival_injection)

t_nodes_total, u_nodes_total = process_output_for_plotting(times_journey, t_nodes, u_nodes)

x_initial_departure_cartesian = locations_journey[:, 1]
x_initial_departure_cartesian[4:6] .+= Δv_departure_injection[1]
x_initial_departure_cartesian = vcat(x_initial_departure_cartesian, mass_starting)
# x_initial_departure_cartesian = vcat(x_initial_departure_cartesian, mass_starting - 40/m_scale)


plot_thrust_information(
    t_nodes_total, 
    u_nodes_total, 
    x_initial_departure_cartesian,
    times_journey, 
    id_journey, 
    mass_changes,
    MassConfig();
    return_plt = false
)

plot_ship_trajectory_low_thrust(
    t_nodes_total, 
    u_nodes_total, 
    x_initial_departure_cartesian, 
    times_journey, 
    id_journey, 
    mass_changes,
    planets_classical,
    MassConfig();
    ga_times, 
    ga_amount,
    zoom = 1.0,
)


for i in 1:length(u_nodes)
    display(maximum(norm.(norm.(eachcol(u_nodes[i][1:3, :])) .- u_nodes[i][4, :])))
end


print("\n\nMass mined: $(-m_scale*mass_changes[end])kg")
print("\n\nFinal mass remaining: $(m_scale*(x_nodes[end][7, end] + mass_changes[end]) - 500)kg")

mass_starting = x_nodes[1][7, 1]

# mass_starting = 1.0

write_solution(
    # "results/testing/Result.txt",
    result_files_all[chosen_ships_all[36]][10],
    mass_starting,
    mass_changes,
    locations_journey,
    times_journey,
    id_journey,
    Δv_departure_injection,
    Δv_arrival_injection,
    t_nodes,
    u_nodes
)


























