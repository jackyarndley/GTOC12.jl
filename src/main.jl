include("header.jl")






journey_length = 9
mixing_number = 2

Δv_cutoff = 50.0/journey_length


times_journey_mixed = generate_default_times(
    journey_length, 
    mixing_number, 
    subset_ids;
    algorithm = 0,
    time_parameter_days = 165
)


id_journey_mixed = mip_ordering_mixing(
    journey_length, 
    subset_ids, 
    times_journey_mixed;
    transfer_dv_limit = Δv_cutoff,
    time_limit_seconds = 80.0,
    include_middle_transfer_dv = true,
    allow_middle_transfer = true,
    self_cleaning = false,
)


for k in 1:length(id_journey_mixed)
    id_journey_mixed[k], times_journey_mixed[k] = pad_long_middle_transfer(
        id_journey_mixed[k], 
        times_journey_mixed[k], 
        id_journey_mixed,
        times_journey_mixed
    )
end

# id_journey = id_journey_mixed[1]
# times_journey = times_journey_mixed[1]





for i in 1:length(times_journey_mixed)
    mass_changes = get_mass_change_at_ids_mixed(
        id_journey_mixed[i],
        times_journey_mixed[i],
        id_journey_mixed,
        times_journey_mixed,
    )

    index_to_remove = findlast(==(0.0), mass_changes)

    if index_to_remove != 1
        times_journey_mixed[i] = vcat(times_journey_mixed[i][1:(index_to_remove - 1)], times_journey_mixed[i][(index_to_remove + 1):end])
    end
end

id_journey_mixed = mip_ordering_mixing(
    journey_length, 
    subset_ids, 
    times_journey_mixed;
    transfer_dv_limit = Δv_cutoff,
    time_limit_seconds = 80.0,
    include_middle_transfer_dv = true,
    # allow_middle_transfer = false
    # self_cleaning = true
)

for k in 1:length(id_journey_mixed)
    id_journey_mixed[k], times_journey_mixed[k] = pad_long_middle_transfer(
        id_journey_mixed[k], 
        times_journey_mixed[k], 
        id_journey_mixed,
        times_journey_mixed
    )
end
















times_journey_mixed = Vector{Float64}[]
id_journey_mixed = Vector{Int64}[]


for i in 1:2
    for algorithm in [2]
    # for algorithm in [0, 1, 2, 2]
        times_journey_mixed_temp = generate_default_times(
            journey_length, 
            mixing_number, 
            subset_ids;
            algorithm,
            # time_parameter_days = 100+i*20
            time_parameter_days = 125+i*20
            # time_parameter_days = 135+i*20
        )

        id_journey_mixed_temp = mip_ordering_mixing(
            journey_length, 
            subset_ids, 
            times_journey_mixed_temp;
            transfer_dv_limit = Δv_cutoff,
            time_limit_seconds = 40.0,
            include_middle_transfer_dv = true,
            # allow_middle_transfer = false
            self_cleaning = true
        )

        for k in 1:length(id_journey_mixed_temp)
            id_journey_mixed_temp[k], times_journey_mixed_temp[k] = pad_long_middle_transfer(
                id_journey_mixed_temp[k], 
                times_journey_mixed_temp[k], 
                id_journey_mixed_temp,
                times_journey_mixed_temp
            )
        end

        id_journey_mixed = vcat(id_journey_mixed, id_journey_mixed_temp)
        times_journey_mixed = vcat(times_journey_mixed, times_journey_mixed_temp)
    end
end




for i in 1:length(times_journey_mixed)
    mass_changes = get_mass_change_at_ids_mixed(
        id_journey_mixed[i],
        times_journey_mixed[i],
        id_journey_mixed,
        times_journey_mixed,
    )

    index_to_remove = findlast(==(0.0), mass_changes)

    if index_to_remove != 1
        times_journey_mixed[i] = vcat(times_journey_mixed[i][1:(index_to_remove - 1)], times_journey_mixed[i][(index_to_remove + 1):end])
    end

    id_journey_mixed[i] = mip_ordering_mixing(
        journey_length, 
        subset_ids, 
        [times_journey_mixed[i]];
        transfer_dv_limit = Δv_cutoff,
        time_limit_seconds = 40.0,
        include_middle_transfer_dv = true,
        self_cleaning = true
        # allow_middle_transfer = false
    )[1]

    id_journey_mixed[i], times_journey_mixed[i] = pad_long_middle_transfer(
        id_journey_mixed[i], 
        times_journey_mixed[i], 
        id_journey_mixed,
        times_journey_mixed
    )
end




plot_trajectory_lambert_information(
    id_journey_mixed[2],
    times_journey_mixed[2]
)






# id_journey_mixed = mip_ordering_mixing(
#     journey_length, 
#     subset_ids, 
#     times_journey_mixed;
#     transfer_dv_limit = Δv_cutoff,
#     time_limit_seconds = 40.0,
#     include_middle_transfer_dv = true,
#     self_cleaning = true
#     # allow_middle_transfer = false
# )

# for i in 1:length(times_journey_mixed)
#     id_journey_mixed[i], times_journey_mixed[i] = pad_long_middle_transfer(
#         id_journey_mixed[i], 
#         times_journey_mixed[i], 
#         id_journey_mixed,
#         times_journey_mixed
#     )
# end









for i in 1:length(times_journey_mixed)
    mass_changes = get_mass_change_at_ids_mixed(
        id_journey_mixed[i],
        times_journey_mixed[i],
        id_journey_mixed,
        times_journey_mixed,
    )

    index_to_remove = findlast(==(0.0), mass_changes)

    if index_to_remove != 1
        times_journey_mixed[i] = vcat(times_journey_mixed[i][1:(index_to_remove - 1)], times_journey_mixed[i][(index_to_remove + 1):end])
    end
end


temp_a = deepcopy(times_journey_mixed)
temp_b = deepcopy(id_journey_mixed)




for (k, subset_ids) in enumerate(unique.([id_journey[2:end-1] for id_journey in id_journey_all]))
    journey_length = min(9, length(subset_ids))
    Δv_cutoff = 150.0/journey_length

    id_journey_mixed[k] = mip_ordering_mixing(
        journey_length, 
        subset_ids, 
        [times_journey_mixed[k]];
        transfer_dv_limit = Δv_cutoff,
        time_limit_seconds = 40.0,
        include_middle_transfer_dv = false,
        self_cleaning = true
        # allow_middle_transfer = false
    )[1]

    id_journey_mixed[k], times_journey_mixed[k] = pad_long_middle_transfer(
        id_journey_mixed[k], 
        times_journey_mixed[k], 
        [id_journey_mixed[k]],
        [times_journey_mixed[k]]
    )
end








































































id_journey = [
    [0, 56608, 56608, -3]
]

times_journey = [
    [1.0321259370000504, 10.827001079130554, 85.01277301090416, 92.92573852790454]
]


mass_changes = get_mass_change_at_ids_mixed(
    id_journey_mixed[1],
    times_journey_mixed[1],
    id_journey_mixed,
    times_journey_mixed
)


locations_journey, Δv_journey_combined, Δv_journey_departure, Δv_journey_arrival = get_lambert_trajectory(
    id_journey_mixed[1],
    times_journey_mixed[1],
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







x_nodes, u_nodes, Δv_departure_injection, Δv_arrival_injection, termination, maximum_error = solve_scp_full_mixed_adaptive_time(
    x_departure_cartesian_mixed,
    x_arrival_cartesian_mixed,
    t_nodes_mixed,
    u_nodes_mixed,
    Δv_departure_injection_mixed,
    Δv_departure_injection_limit_mixed,
    Δv_arrival_injection_mixed,
    Δv_arrival_injection_limit_mixed,
    mass_changes_mixed,
    id_journey_mixed,
    LoggedMassConfig();
    linearization_error = 1e-4,
    display_output = true,
    reverse_mass,
    initial_scaling_factor = 0.025,
    # optimizer = Gurobi.Optimizer,
    optimizer = Mosek.Optimizer,
    time_change = true
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
    id_journey,
    times_journey,
    mass_changes,
    t_nodes_total, 
    u_nodes_total, 
    x_initial_departure_cartesian,
    MassConfig()
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


























