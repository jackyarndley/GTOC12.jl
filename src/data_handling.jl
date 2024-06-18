using DelimitedFiles, LinearAlgebra

function read_file_classical_data(filename)
    classical_data = readdlm(filename; skipstart=1)
    classical_data =  classical_data[:, 3:end]
    classical_data[:, 3:end] .*= π/180

    return classical_data'
end

function remove_gravity_assists_from_journeys(id_journey_all)
    id_journey_all_removed = deepcopy(id_journey_all)

    for i in 1:length(id_journey_all_removed)
        temp = filter(>=(0), id_journey_all_removed[i])

        push!(temp, -3)

        id_journey_all_removed[i] = temp
    end

    return id_journey_all_removed
end


# Format for the files
# write the mass before and after 
# shipid asteroidid time state mass_before_dropoff
# shipid asteroidid time state mass_after_dropoff

# shipid -1 time_first_node 0.0 0.0 0.0 
# shipid -1 time_first_node thrust1

# shipid -1 time_second_node thrust1
# shipid -1 time_second_node thrust2

# ...
# shipid -1 time_last_node thrustlast-1
# shipid -1 time_last_node thrustlast
# shipid -1 time_end 0.0 0.0 0.0

# shipid asteroidid time state mass_before_dropoff
# shipid asteroidid time state mass_after_dropoff

function write_solution(
    filename,
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
    # Create file if does not exist
    touch(filename)

    open(filename, "w") do f
        ship_mass = copy(mass_starting)

        x_current = locations_journey[:, 1]
        x_current = vcat(x_current, [ship_mass])

        times_mjd_journey = convert_time_to_mjd(times_journey)

        maximum_position_error = 0.0

        for ship_id in 1:1
            for i in 1:size(id_journey, 1)
                target_id = id_journey[i]
                
                if !(mass_changes[i] ≈ 0 && target_id > 0)
                    # Apply the arrival injection state
                    if i > 1
                        x_current[4:6] -= Δv_arrival_injection[i - 1]
                    end

                    # Write arrival state
                    write(f, "$ship_id $target_id $(times_mjd_journey[i]) $(convert_to_actual_state_string(x_current))\n")

                    # Apply departure injection state
                    if i < size(id_journey, 1)
                        x_current[4:6] .= locations_journey[4:6, i] + Δv_departure_injection[i]
                    end

                    # Apply mass dropoff/pickup
                    x_current[7] += mass_changes[i]

                    # If we are at the last state terminate
                    if i == size(id_journey, 1)
                        write(f, "$ship_id $target_id $(times_mjd_journey[i]) $(convert_to_actual_state_string(x_current))")
                        break
                    end

                    # Write departure state
                    write(f, "$ship_id $target_id $(times_mjd_journey[i]) $(convert_to_actual_state_string(x_current))\n")
                end

                # Normalise thrust to the correct value
                u_nodes[i][4, :] = norm.(eachcol(u_nodes[i][1:3, :]))

                # Propagate to get the final mass
                x_nodes = integrate_trajectory(
                    x_current, 
                    t_nodes[i]; 
                    t_nodes = t_nodes[i], 
                    u_nodes = u_nodes[i], 
                    thrust_config = ZeroOrderHoldConfig(),
                    objective_config = MassConfig()
                )

                maximum_position_error = max(maximum_position_error, r_scale*norm(x_nodes[1:3, end] - locations_journey[1:3, i + 1]))

                x_current = copy(x_nodes[:, end])

                x_current[1:6] = locations_journey[:, i + 1]

                # Write the thrusting lines
                thrust_nodes = u_nodes[i][1:3, :] * thrust * m_scale * a_scale * 1e3

                current_time = convert_time_to_mjd(times_journey[i] + t_nodes[i][1])

                for j in 1:(length(t_nodes[i]) - 1)
                    write(f, "$ship_id -1 $current_time 0.0, 0.0, 0.0\n")

                    current_thrust = thrust_nodes[:, j]

                    if norm(current_thrust) > thrust * m_scale * a_scale * 1e3
                        current_thrust = normalize(current_thrust) * thrust * m_scale * a_scale * 1e3
                    end

                    write(f, "$ship_id -1 $current_time $(current_thrust[1]) $(current_thrust[2]) $(current_thrust[3])\n")

                    next_time = convert_time_to_mjd(times_journey[i] + t_nodes[i][j + 1])

                    if j + 1 == length(t_nodes[i])
                        next_time = convert_time_to_mjd(times_journey[i + 1])
                    end

                    while !(current_time ≈ next_time)
                        current_time += min(1.0, next_time - current_time)
                        write(f, "$ship_id -1 $current_time $(current_thrust[1]) $(current_thrust[2]) $(current_thrust[3])\n")
                    end

                    # current_time = next_time

                    write(f, "$ship_id -1 $current_time 0.0, 0.0, 0.0\n")
                end
            end
        end

        print("\nMaximum position error: $(maximum_position_error)km")
    end
end


function combine_solutions(
    filenames,
    output_filename
)
    touch(output_filename)

    open(output_filename, "w") do f1
        for (i, file) in enumerate(filenames)
            data = open(file) do f2
                readlines(f2)
            end

            for line in data
                write(f1, "$(i) ")
                # write(f1, convert(String, split(line, " "; limit=2)[2]))
                write(f1, convert(String, Base.split(line, " "; limit=2)[2]))
                write(f1, "\n")
            end
        end
    end

    temp = open(output_filename, "r") do f1
        convert(String, chomp(read(f1, String)))
    end

    open(output_filename, "w") do f1
        write(f1, temp)
    end

    return
end
