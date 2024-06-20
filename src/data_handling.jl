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
    p::SequentialConvexProblem,
    filename
)
    touch(filename)

    open(filename, "w") do f
        for n in 1:p.mixing_number
            for k in 1:(p.segment_number[n] + 1)

                # Write current state
                if k == 1
                    write(f, "$n $(p.id_journey[n][k]) $(convert_time_to_mjd(p.times_journey[n][k])) $(convert_to_actual_state_string(p.x0[n][k]))\n")
                else
                    temp = deepcopy(p.xf[n][k - 1])
                    temp[4:6] .-= p.Δvf[n][k - 1]
                    # temp[4:6] = p.x_nodes[n][k - 1][4:6, end]

                    write(f, "$n $(p.id_journey[n][k]) $(convert_time_to_mjd(p.times_journey[n][k])) $(convert_to_actual_state_string(temp))\n")
                end

                if k == (p.segment_number[n] + 1)
                    temp = deepcopy(p.xf[n][end])
                    temp[4:6] .-= p.Δvf[n][end]
                    temp[7] += p.Δm0[n][k]
        
                    write(f, "$n $(p.id_journey[n][k]) $(convert_time_to_mjd(p.times_journey[n][k])) $(convert_to_actual_state_string(temp))")

                    if n != p.mixing_number
                        write(f, "\n")
                    end

                    break
                end

                temp = deepcopy(p.x0[n][k])
                temp[4:6] += p.Δv0[n][k]

                write(f, "$n $(p.id_journey[n][k]) $(convert_time_to_mjd(p.times_journey[n][k])) $(convert_to_actual_state_string(temp))\n")


                for j in 1:(length(p.t_nodes[n][k]) - 1)
                    current_time = convert_time_to_mjd(p.times_journey[n][k] + p.t_nodes[n][k][j])
                    current_thrust = deepcopy(p.u_nodes[n][k][1:3, j]) * thrust * m_scale * a_scale * 1e3
                    next_time = min(
                        convert_time_to_mjd(p.times_journey[n][k] + p.t_nodes[n][k][j + 1]),
                        convert_time_to_mjd(p.times_journey[n][k + 1])
                    )

                    # Limit thrust in the case of small overhead
                    if norm(current_thrust) > thrust * m_scale * a_scale * 1e3
                        current_thrust = normalize(current_thrust) * thrust * m_scale * a_scale * 1e3
                    end

                    write(f, "$n -1 $current_time 0.0, 0.0, 0.0\n")

                    for time in current_time:1.0:next_time
                        current_time = time
                        write(f, "$n -1 $current_time $(current_thrust[1]) $(current_thrust[2]) $(current_thrust[3])\n")
                    end

                    if !(current_time ≈ next_time)
                        write(f, "$n -1 $next_time $(current_thrust[1]) $(current_thrust[2]) $(current_thrust[3])\n")
                    end

                    write(f, "$n -1 $next_time 0.0, 0.0, 0.0\n")
                end
            end
        end

        # print("\nMaximum position error: $(maximum_position_error)km")
    end
end



