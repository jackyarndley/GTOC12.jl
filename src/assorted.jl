function create_acceleration_guess_from_lambert_spread(Δv_departure, Δv_arrival, Δv_injection, Δt, nodes)
    half_nodes = Int(round(nodes/2))

    u_departure = (Δv_departure./Δt - Δv_injection./Δt) / half_nodes
    u_arrival = (Δv_arrival./Δt) / half_nodes

    u_nodes = hcat(
        repeat(u_departure, outer = [1, half_nodes]),
        repeat(u_arrival, outer = [1, nodes - half_nodes])
    ) ./ thrust/maximum_mass

    u_nodes = vcat(
        u_nodes,
        stack(norm.(eachcol(u_nodes)))'
    )

    return u_nodes
end

function control_guess_from_lambert_impulsive(Δv_departure, Δv_arrival, Δv_departure_injection, Δv_arrival_injection, Δt, nodes)
    u_nodes = zeros(3, nodes)

    u_departure = Δv_departure - Δv_departure_injection
    u_arrival = Δv_arrival - Δv_arrival_injection

    # Get the acceleration needed
    u_nodes[:, 1] = u_departure ./ Δt[1]
    u_nodes[:, end] = u_arrival ./ Δt[end]

    u_nodes = vcat(
        u_nodes,
        stack(norm.(eachcol(u_nodes)))'
    )
    
    return u_nodes ./ thrust
end

function calculate_Δv_injection_from_lambert(Δv_injection, injection_condition; limit = 6.0/v_scale)
    Δv_injection = copy(Δv_injection)

    Δv_injection_limit = if injection_condition
        if norm(Δv_injection) > limit
            Δv_injection = normalize(Δv_injection) * limit
        end

        limit
    else
        Δv_injection .*= 0.0

        0.0
    end

    return Δv_injection, Δv_injection_limit
end

function convert_mjd_to_time(mjd)
    return (mjd .- mjd_start).*day_scale
end

function convert_time_to_mjd(time)
    return time./day_scale .+ mjd_start
end

function convert_to_state_string(x)
    x = copy(x)

    x[1:3, :] .*= r_scale
    x[4:6, :] .*= v_scale
    x[7, :] .= x[7, :] * m_scale

    return "$(x[1]) $(x[2]) $(x[3]) $(x[4]) $(x[5]) $(x[6]) $(x[7])"
end

function ephermeris_cartesian_from_id(id, times)
    if id < 0
        ephemeris_cartesian_at(planets_classical[:, -id - 1], times)
    elseif id == 0
        ephemeris_cartesian_at(planets_classical[:, 2], times)
    else
        ephemeris_cartesian_at(asteroids_classical[:, id], times)
    end
end

function find_best_lambert_transfer(
    x_start_journey,
    x_end_journey,
    tof
)
    n_revolutions = 0
    current_dv = 1e5
    best_dv = 1e6
    best_transfer = nothing

    while n_revolutions < 1
        best_dv = current_dv

        transfer = multi_revolution_lambert(
            x_start_journey,
            x_end_journey,
            tof,
            n_revolutions
        )

        current_dv = norm(transfer[1] - x_start_journey[4:6, :]) + norm(x_end_journey[4:6, :] - transfer[2])

        if current_dv > best_dv
            break
        end

        best_dv = current_dv
        best_transfer = transfer

        n_revolutions += 1
    end
    
    return best_transfer[1] - x_start_journey[4:6, :], x_end_journey[4:6, :] - best_transfer[2]
end

function get_lambert_trajectory(
    times_journey,
    id_journey
)
    Δv_journey_combined = zeros(Float64, 3, 0)
    Δv_journey_departure = zeros(Float64, 3, 0)
    Δv_journey_arrival = zeros(Float64, 3, 0)

    x_start_journey = ephermeris_cartesian_from_id(id_journey[1], times_journey[1])[:]

    locations_journey = hcat(x_start_journey)

    for i in 2:length(id_journey)
        x_end_journey = ephermeris_cartesian_from_id(id_journey[i], times_journey[i])[:]

        tof = times_journey[i] - times_journey[i-1]

        Δv_departure, Δv_arrival = if id_journey[i] == id_journey[i-1]
            zeros(Float64, 3), zeros(Float64, 3)
        else
            temp = deepcopy(x_end_journey)

            # Prevent near 180° transfers
            while norm(cross(x_start_journey[1:3], temp[1:3])) < 0.6
                tof += 1*day_scale
                temp = ephermeris_cartesian_from_id(id_journey[i], times_journey[i - 1] + tof)
            end

            find_best_lambert_transfer(x_start_journey, temp, tof)
        end

        if i == 2
            Δv_journey_combined = hcat(Δv_journey_combined, Δv_departure)
        else
            Δv_journey_combined[:, i-1] .+= Δv_departure
        end
        
        Δv_journey_combined = hcat(Δv_journey_combined, Δv_arrival)
        Δv_journey_departure = hcat(Δv_journey_departure, Δv_departure)
        Δv_journey_arrival = hcat(Δv_journey_arrival, Δv_arrival)

        locations_journey = hcat(locations_journey, x_end_journey)
        x_start_journey = x_end_journey
    end

    return locations_journey, Δv_journey_combined, Δv_journey_departure, Δv_journey_arrival
end

function get_mass_change_at_ids(
    id_journey,
    times_journey
)
    mass_change = zeros(Float64, length(id_journey))
    id_visited = []

    for i in eachindex(id_journey)
        if id_journey[i] in id_visited
            first_index = findfirst(==(id_journey[i]), id_journey)
            last_index = findlast(==(id_journey[i]), id_journey)

            mass_change[i] = if i != last_index
                0.0
            else
                (times_journey[i] - times_journey[first_index]) * mining_rate
            end
        elseif id_journey[i] > 0
            mass_change[i] = -miner_mass
            push!(id_visited, id_journey[i])
        elseif id_journey[i] == -3 && i == length(id_journey)
            mass_change[i] = -sum(mass_change[mass_change .> 0.0])
        end
    end

    return mass_change
end

function pad_long_middle_transfer(
    id_journey,
    times_journey,
    id_journey_all,
    times_journey_all,
)
    mass_changes = get_mass_change_at_ids_mixed(
        id_journey,
        times_journey,
        id_journey_all,
        times_journey_all,
    )

    i = findfirst(>(0.0), mass_changes) - 1

    if (times_journey[i + 1] - times_journey[i] > 1200*day_scale) && (id_journey[i] != id_journey[i+1])
        id_journey = vcat(id_journey[1:i], [id_journey[i]], id_journey[i+1:end])
        times_journey = vcat(times_journey[1:i], [times_journey[i+1] - 800*day_scale], times_journey[i+1:end])
    end

    return id_journey, times_journey
end

function is_self_cleaning(id_journey)
    for id in id_journey
        if id <= 0
            continue
        end

        if length(findall(==(id), id_journey)) == 1
            return false
        end
    end

    return true
end

function get_mass_change_at_ids_mixed(
    id_journey,
    times_journey,
    id_journey_all,
    times_journey_all,
)
    mass_change = zeros(Float64, length(id_journey))

    if is_self_cleaning(id_journey)
        return get_mass_change_at_ids(id_journey, times_journey)
    end

    for i in eachindex(id_journey)
        times = []

        if id_journey[i] == 0
            continue
        elseif id_journey[i] == -3
            mass_change[i] = -sum(mass_change[mass_change .> 0.0])
            continue
        end

        for j in eachindex(id_journey_all)
            if is_self_cleaning(id_journey_all[j])
                continue
            end

            for k in eachindex(id_journey_all[j])
                if id_journey_all[j][k] == id_journey[i]
                    push!(times, times_journey_all[j][k])
                end
            end
        end

        if minimum(times) ≈ times_journey[i]
            mass_change[i] = -miner_mass
        elseif maximum(times) ≈ times_journey[i]
            mass_change[i] = (maximum(times) - minimum(times)) * mining_rate
        end
    end

    return mass_change
end

function convert_log_mass_control_to_mass_control(
    locations_journey,
    mass_starting,
    t_nodes,
    u_nodes,
    Δv_departure_injection,
    mass_changes,
    x_departure_cartesian
)
    # Final refinement with non-logged formulation so that the thrust is constant over an arc
    x_current = locations_journey[:, 1]
    x_current = vcat(x_current, [log(mass_starting)])
    x_current[4:6] .+= Δv_departure_injection[1]
    x_current[7] = log(exp(x_current[7]) + mass_changes[1]) 

    for i in 1:length(u_nodes)
        x_departure_cartesian[i][7] = exp(x_current[7])

        x_nodes_segment = propagate_spacecraft(
            x_current, 
            t_nodes[i]; 
            t_nodes = t_nodes[i], 
            u_nodes = u_nodes[i], 
            thrust_config = ZeroOrderHoldConfig(),
            objective_config = LoggedMassConfig()
        )

        x_current = x_nodes_segment[:, end]

        x_current[7] = log(exp(x_current[7]) + mass_changes[i + 1])

        x_nodes_segment[7, :] = exp.(x_nodes_segment[7, :])

        u_nodes[i] = u_nodes[i][1:4, :]
        u_nodes[i][1, :] .*= x_nodes_segment[7, 1:(end-1)]
        u_nodes[i][2, :] .*= x_nodes_segment[7, 1:(end-1)]
        u_nodes[i][3, :] .*= x_nodes_segment[7, 1:(end-1)]
        u_nodes[i][4, :] .*= x_nodes_segment[7, 1:(end-1)]
    end

    return u_nodes, x_departure_cartesian
end

function get_gravity_assist_information(
    times_journey,
    id_journey,
    Δv_departure_injection,
    Δv_arrival_injection;
    display_output = true
)
    ga_times = []
    ga_amount = []
    
    if display_output
        @printf "\n%8s  %5s  %9s  %12s  %12s  %6s  %6s\n" "TIME" "PLACE" "V∞" "RP MIN" "RP ACT" "A MAX" "A ACT"
    end

    for i in 1:segments
        if id_journey[i] < 0
            push!(ga_times, times_journey[i])
            push!(ga_amount, Δv_arrival_injection[i - 1] + Δv_departure_injection[i])
    
            name = ["Venus", "Earth", "Mars"][-id_journey[i] - 1]
            μ_gravity_assist = μ_planets[-id_journey[i] .- 1]
            rp_min = rp_min_planets_actual[-id_journey[i] .- 1]
    
            vinf = norm(Δv_arrival_injection[i - 1])
            rp_actual = get_gravity_assist_radius(Δv_arrival_injection[i - 1], Δv_departure_injection[i], μ_gravity_assist)
    
            θ_max = rad2deg(get_gravity_assist_angle(Δv_arrival_injection[i - 1], μ_gravity_assist, rp_min))
            θ_actual = rad2deg(get_gravity_assist_angle(Δv_arrival_injection[i - 1], μ_gravity_assist, rp_actual))
    
            if display_output
                @printf "%7.2f  %5s  %5.2fkm/s  %10.2fkm  %10.2fkm  %3.2f°  %3.2f°\n" convert_time_to_mjd(times_journey[i]) name vinf*v_scale rp_min*r_scale rp_actual*r_scale θ_max θ_actual
            end
        end
    end

    return ga_times, ga_amount
end

function process_output_for_plotting(times_journey, t_nodes, u_nodes)
    t_nodes_total = vcat([t_nodes[i][1:(end-1)] .+ times_journey[i] for i in 1:segments]...)
    u_nodes_total = hcat(u_nodes...)
    
    t_nodes_total = vcat(t_nodes_total, [times_journey[segments + 1]])
    t_nodes_total = vcat(t_nodes_total, [times_journey[segments + 1]])
    u_nodes_total = hcat(u_nodes_total, zeros(size(u_nodes_total, 1)))
    u_nodes_total = hcat(u_nodes_total, zeros(size(u_nodes_total, 1)))
    
    return t_nodes_total, u_nodes_total
end

function get_lambert_guess_for_scp(
    Δv_journey_departure,
    Δv_journey_arrival,
    locations_journey,
    id_journey,
    times_journey,
    mass_current,
    mass_changes;
    reverse_mass = false,
    ignore_long_transfers = true,
)
    Δv_journey_departure = copy(Δv_journey_departure)
    Δv_journey_arrival = copy(Δv_journey_arrival)

    x_departure_cartesian = fill(Float64[], segments)
    x_arrival_cartesian = fill(Float64[], segments)
    t_nodes = fill(Float64[], segments)
    u_nodes = fill(Matrix{Float64}(undef, 0, 0), segments)
    Δv_departure_injection = fill(Float64[], segments)
    Δv_departure_injection_limit = fill(0.0, segments)
    Δv_arrival_injection = fill(Float64[], segments)
    Δv_arrival_injection_limit = fill(0.0, segments)
    
    for i in 1:segments
        # When a transfer is very long the Lambert is not a good guess so just use a zero one
        if ignore_long_transfers && ((mass_changes[i] < 0.0) && (mass_changes[i+1] > 0.0))
            Δv_journey_departure[:, i] *= 0.0
            Δv_journey_arrival[:, i] *= 0.0
        end

        # Only allow departure or arrival Δv at the beginning and end
        Δv_departure_injection[i], Δv_departure_injection_limit[i] = calculate_Δv_injection_from_lambert(Δv_journey_departure[:, i], id_journey[i] == 0)
        Δv_arrival_injection[i], Δv_arrival_injection_limit[i] = calculate_Δv_injection_from_lambert(Δv_journey_arrival[:, i], (i + 1) == length(id_journey) && id_journey[i+1] == -3)

        # Gravity assist departure
        if id_journey[i] < 0
            Δv_departure_injection[i], Δv_departure_injection_limit[i] = calculate_Δv_injection_from_lambert(Δv_journey_departure[:, i], true; limit = 50.0 / v_scale)
        end
    
        # Gravity assist arrival
        if id_journey[i + 1] < 0 && (i + 1) != length(id_journey)
            Δv_arrival_injection[i], Δv_arrival_injection_limit[i] = calculate_Δv_injection_from_lambert(Δv_journey_arrival[:, i], true; limit = 50.0 / v_scale)
        end
    end

    for i in 1:segments
        t_max = times_journey[i+1] - times_journey[i]
        t_nodes_segment = collect(0.0:node_time_spacing:(t_max - t_max % node_time_spacing))

        if !(t_nodes_segment[end] ≈ t_max)
            t_nodes_segment = vcat(t_nodes_segment, t_max)
        end
    
        if id_journey[i] == id_journey[i + 1]
            t_nodes_segment = [0.0, t_max/2, t_max]
        end
        
        nodes = length(t_nodes_segment)

        # Get mass estimate with lambert burn. In reverse mode this is more than before
        mass_current, mass_end = if reverse_mass
            mass_current*exp((norm(Δv_journey_departure[:, i] .- Δv_departure_injection[i]) + norm(Δv_journey_arrival[:, i] .- Δv_arrival_injection[i])) / g0_isp), mass_current
        else
            mass_current, mass_current/exp((norm(Δv_journey_departure[:, i] .- Δv_departure_injection[i]) + norm(Δv_journey_arrival[:, i] .- Δv_arrival_injection[i])) / g0_isp)
        end
    
        u_nodes_segment = control_guess_from_lambert_impulsive(
            Δv_journey_departure[:, i], 
            Δv_journey_arrival[:, i],
            Δv_departure_injection[i],
            Δv_arrival_injection[i],
            t_nodes_segment[2:end] .- t_nodes_segment[1:(end-1)],
            nodes - 1
        )

        x_departure_cartesian[i] = vcat(locations_journey[:, i], mass_current)
        x_arrival_cartesian[i] = vcat(locations_journey[:, i+1], mass_end)
    
        t_nodes[i] = t_nodes_segment
        u_nodes[i] = u_nodes_segment
    
        # Apply the mass change for the linking constraints
        if reverse_mass
            mass_current = mass_current - mass_changes[i + 1]
        else
            mass_current = mass_end + mass_changes[i + 1]
        end
    end
    
    return x_departure_cartesian, x_arrival_cartesian, t_nodes, u_nodes, Δv_departure_injection, Δv_departure_injection_limit, Δv_arrival_injection, Δv_arrival_injection_limit
end

function get_gravity_assist_radius(va, vd, μ)
    δ = acos(dot(-normalize(va), normalize(vd)))
    e = 1/sin(δ/2)

    vinf = min(norm(va), norm(vd))

    return μ*(e - 1)/vinf^2
end

function get_maximum_gravity_assist_angle(v_ref, μ, rp_min)
    v_inf_2 = sum([a^2 for a in v_ref])

    return 2*asin((μ/rp_min)/(v_inf_2 + (μ/rp_min)))
end

function get_maximum_gravity_assist_deviation(v_ref, μ, rp_min)
    v_inf_2 = sum([a^2 for a in v_ref])

    θ = 2*asin((μ/rp_min)/(v_inf_2 + (μ/rp_min)))

    return sqrt(v_inf_2*(2 - 2*cos(θ)))
end

function get_gravity_assist_angle(v, μ, rp)
    v_inf_2 = sum([a^2 for a in v])

    return 2*asin((μ/rp)/(v_inf_2 + (μ/rp)))
end

function get_norm_linearisation(v)
    stm_function(values) = begin
        return [norm(values)]
    end

    return ForwardDiff.jacobian(
        stm_function, 
        vcat(v)
    )
end

function get_angle_linearisation(va, vd, μ)
    stm_function(values) = begin
        va = values[1:3]
        vd = values[4:6]
        rp = get_gravity_assist_radius(va, vd, μ)
        temp = get_gravity_assist_angle(va, μ, rp)
        temp = min(temp, get_gravity_assist_angle(vd, μ, rp))
        return [temp]
    end

    return ForwardDiff.jacobian(
        stm_function, 
        vcat(va, vd)
    )
end

function convert_time_to_5day_time(
    time
)
    return Int(floor((time / day_scale) / 5)) + 1
end

function convert_mjd_to_5day_time(
    mjd_time
)
    time = convert_mjd_to_time(mjd_time)
    return convert_time_to_5day_time(time)
end

function find_first_repeating_index(
    id_journey
)
    for i in 1:(length(id_journey) - 1)
        if id_journey[i] == id_journey[i + 1]
            return i
        end
    end
end

function load_result_files(result_files; low_thrust_dv = false)
    bonus_optimal = true

    bonus = if bonus_optimal
        open("results/bonuses/bonus_coefficients.txt") do f
            [parse(Float64, split(line)[1]) for line in readlines(f)]
        end
    else
        ones(60000)
    end

    ship_number = 0
    id_journey_all = Vector{Vector{Int64}}()
    times_journey_all = Vector{Vector{Float64}}()
    mined_mass_all = Vector{Float64}()
    penalised_mass_all = Vector{Float64}()
    low_thrust_dv_all = Vector{Vector{Float64}}()
    fuel_remaining_all = Vector{Float64}()

    for file_name in result_files
        temp = open(file_name) do f
            readlines(f)
        end

        filter!(x->split(x)[2]!="-1", temp)

        local_number = 0
        previous_mass = 0.0

        for i in 1:2:length(temp)
            temp2 = split(temp[i])

            if parse(Int64, temp2[1]) != local_number
                local_number = parse(Int64, temp2[1])
                ship_number += 1

                push!(id_journey_all, Int64[])
                push!(times_journey_all, Float64[])
                push!(mined_mass_all, 0.0)
                push!(penalised_mass_all, 0.0)
                push!(low_thrust_dv_all, Float64[])
                push!(fuel_remaining_all, 0.0)
            end

            temp3 = split(temp[i + 1])

            m1 = parse(Float64, temp2[end])
            m2 = parse(Float64, temp3[end])

            if parse(Int64, temp2[2]) != 0
                push!(low_thrust_dv_all[ship_number], g0_isp*log(previous_mass/m1))
            end

            previous_mass = m2

            # Pickup
            if m1 < m2
                mined_mass_all[end] += m2 - m1
                penalised_mass_all[end] += bonus[parse(Int64, temp2[2])] * (m2 - m1)
            end

            fuel_remaining_all[end] = m2 - 500

            push!(id_journey_all[ship_number], parse(Int64, temp2[2]))
            push!(times_journey_all[ship_number], convert_mjd_to_time(parse(Float64, temp2[3])))
        end
    end 

    return if low_thrust_dv
        id_journey_all, times_journey_all, mined_mass_all, penalised_mass_all, low_thrust_dv_all, fuel_remaining_all
    else
        id_journey_all, times_journey_all, mined_mass_all, penalised_mass_all, fuel_remaining_all
    end
    
end

function load_result_folders_grouping(result_folders)
    id_journey_all = Vector{Vector{Int64}}()
    times_journey_all = Vector{Vector{Float64}}()
    mined_mass_all = Vector{Float64}()
    penalised_mass_all = Vector{Float64}()
    fuel_remaining_all = Vector{Float64}()
    groups_all = Vector{Int64}()
    result_files_all = []

    for folder in result_folders
        if folder[(end-3):end] == ".txt"
            result_files = [folder]
        else
            result_files = readdir(folder, join=true)
        end

        id_journey_all_folder, times_journey_all_folder, mined_mass_all_folder, penalised_mass_all_folder, fuel_remaining_all_folder = load_result_files(result_files)

        if folder[(end-3):end] == ".txt"
            result_files_all = vcat(result_files_all, repeat(result_files, length(id_journey_all_folder)))
        else
            result_files_all = vcat(result_files_all, result_files)
        end

        id_journey_all = vcat(id_journey_all, id_journey_all_folder)
        times_journey_all = vcat(times_journey_all, times_journey_all_folder)
        mined_mass_all = vcat(mined_mass_all, mined_mass_all_folder)
        penalised_mass_all = vcat(penalised_mass_all, penalised_mass_all_folder)
        fuel_remaining_all = vcat(fuel_remaining_all, fuel_remaining_all_folder)

        group_index = 1 

        groups = zeros(Int64, length(id_journey_all_folder))

        for i in 1:length(id_journey_all_folder)
            # Not assigned to a group
            if groups[i] == 0
                groups[i] = group_index

                group_ids = unique(id_journey_all_folder[i][2:(end-1)])

                filter!(x->count(==(x), id_journey_all_folder[i][2:(end-1)])==1, group_ids)

                previous_length = 0

                while length(group_ids) != previous_length
                    previous_length = length(group_ids)

                    for j in (i+1):length(id_journey_all_folder)
                        ids = unique(id_journey_all_folder[j][2:(end-1)])

                        filter!(x->count(==(x), id_journey_all_folder[j][2:(end-1)])==1, ids)

                        if length(intersect(group_ids, ids)) > 0
                            groups[j] = group_index
                            group_ids = unique(vcat(group_ids, ids))
                        end
                    end
                end

                group_index += 1
            end
        end

        if length(groups_all) > 0
            groups_all = vcat(groups_all, groups .+ maximum(groups_all))
        else
            groups_all = groups
        end
    end

    return id_journey_all, times_journey_all, mined_mass_all, penalised_mass_all, groups_all, result_files_all, fuel_remaining_all
end


function find_best_transfer_time_forward(
    id_start,
    id_end,
    search_start_time,
)
    minimum_Δv_norm = 1000.0
    minimum_Δv_time = 0.0
    minimum_Δv_duration = 0.0

    for start_time in search_start_time:10*day_scale:(search_start_time + 300*day_scale)
        starting_location = ephermeris_cartesian_from_id(id_start, start_time)
        times = collect(500:50:700).*day_scale

        for time in times
            ending_location = ephermeris_cartesian_from_id(id_end, start_time + time)

            lambert_departure, lambert_arrival = zero_revolution_lambert(
                starting_location,
                ending_location,
                time
            )

            departure_dv = lambert_departure - starting_location[4:6]
            arrival_dv = ending_location[4:6] - lambert_arrival
        
            departure_dv_norm = norm.(eachcol(departure_dv))
            arrival_dv_norm = norm.(eachcol(arrival_dv))

            if id_start == 0
                departure_dv_norm = max.(departure_dv_norm .- 6.0/v_scale, 0.0)
            end

            if id_end == -3
                arrival_dv_norm = max.(arrival_dv_norm .- 6.0/v_scale, 0.0)
            end
        
            total_dv_norm = departure_dv_norm + arrival_dv_norm

            if total_dv_norm[1] < minimum_Δv_norm
                minimum_Δv_time = start_time
                minimum_Δv_duration = time
                minimum_Δv_norm = total_dv_norm[1]
            end
        end
    end

    # Add a bit to the duration because the first transfer is usually very low thrust
    return minimum_Δv_time, minimum_Δv_duration + 100*day_scale, minimum_Δv_norm
end

function find_best_transfer_time_reverse(
    id_start,
    id_end,
    search_start_time,
)
    # 500-700 day transfers

    minimum_Δv_norm = 1000.0
    minimum_Δv_time = 0.0
    minimum_Δv_duration = 0.0

    for start_time in (search_start_time - 365*day_scale):10*day_scale:search_start_time
        times = collect(500:50:700).*day_scale

        for time in times
            starting_location = ephermeris_cartesian_from_id(id_start, start_time - time)
            ending_location = ephermeris_cartesian_from_id(id_end, start_time)

            lambert_departure, lambert_arrival = zero_revolution_lambert(
                starting_location,
                ending_location,
                time
            )

            departure_dv = lambert_departure - starting_location[4:6]
            arrival_dv = ending_location[4:6] - lambert_arrival
        
            departure_dv_norm = norm.(eachcol(departure_dv))
            arrival_dv_norm = norm.(eachcol(arrival_dv))

            if id_start == 0
                departure_dv_norm = max.(departure_dv_norm .- 6.0/v_scale, 0.0)
            end

            if id_end == -3
                arrival_dv_norm = max.(arrival_dv_norm .- 6.0/v_scale, 0.0)
            end
        
            total_dv_norm = departure_dv_norm + arrival_dv_norm

            if total_dv_norm[1] < minimum_Δv_norm
                minimum_Δv_time = start_time
                minimum_Δv_duration = time
                minimum_Δv_norm = total_dv_norm[1]
            end
        end
    end

    return minimum_Δv_time, minimum_Δv_duration, minimum_Δv_norm
end

function get_transfer_dv(id_start, time_start, id_set, Δt)
    starting_location = ephermeris_cartesian_from_id(id_start, [time_start])

    Δv = Float64[]

    for id in id_set
        if id == id_start
            push!(Δv, 0.0)
            continue
        end

        ending_location = ephermeris_cartesian_from_id(id, time_start + Δt)

        lambert_departure, lambert_arrival = zero_revolution_lambert(
            starting_location,
            ending_location,
            Δt
        )

        if isnan(lambert_departure[1])
            lambert_departure, lambert_arrival = multi_revolution_lambert(
                starting_location,
                ending_location,
                Δt,
                1
            )
        end

        departure_dv = lambert_departure - starting_location[4:6]
        arrival_dv = ending_location[4:6] - lambert_arrival
    
        departure_dv_norm = norm.(eachcol(departure_dv))
        arrival_dv_norm = norm.(eachcol(arrival_dv))

        if id_start == 0
            departure_dv_norm = max.(departure_dv_norm .- 6.0/v_scale, 0.0)
        end

        if id == -3
            arrival_dv_norm = max.(arrival_dv_norm .- 6.0/v_scale, 0.0)
        end
    
        total_dv_norm = departure_dv_norm + arrival_dv_norm

        push!(Δv, total_dv_norm[1])
    end

    return Δv
end

function get_transfer_dv_low_thrust(id_start, time_start, id_set, Δt; start_mass = 1.0)
    Δv = Float64[]

    times_journey = [time_start, time_start + Δt]
    mass_changes = [0.0, 0.0]

    global segments = 1
    global scp_iterations = 10
    global node_time_spacing = 25.0 * day_scale

    for id in id_set
        if id == id_start
            push!(Δv, 0.0)
            continue
        end

        id_journey = [id_start, id]

        locations_journey, Δv_journey_combined, Δv_journey_departure, Δv_journey_arrival = get_lambert_trajectory(
            times_journey,
            id_journey
        )

        x_departure_cartesian, x_arrival_cartesian, t_nodes, u_nodes, Δv_departure_injection, Δv_departure_injection_limit, Δv_arrival_injection, Δv_arrival_injection_limit = get_lambert_guess_for_scp(
            Δv_journey_departure,
            Δv_journey_arrival,
            locations_journey,
            id_journey,
            times_journey,
            start_mass,
            mass_changes
        );

        try
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
                display_output = false,
                linearization_error = 1e-4,
                fix_start_mass = true
            );

            if (termination == OPTIMAL) && (maximum_error < 1e-2)
                push!(Δv, g0_isp*log(start_mass/x_nodes[end][7, end]))
            else
                push!(Δv, 100.0)
            end
        catch
            push!(Δv, 100.0)
        end
            
        # push!(Δv, g0_isp*log(x_nodes[1][7, 1]/x_nodes[end][7, end]))
    end

    return Δv
end

function get_lambert_phase(
    id_start,
    id_targets,
    times,
)
    starting_locations = repeat(asteroids_cartesian_5day[:, id_start, times], 1, size(id_targets, 1))
    total_dv_norm_all = fill(1000.0, length(times))

    # for tof in 20:10:80
    for tof in 20:10:50
        ending_locations = asteroids_cartesian_5day[:, id_targets, times .+ tof]
    
        lambert_departure, lambert_arrival = zero_revolution_lambert(
            starting_locations[1:3, :],
            ending_locations[1:3, :],
            fill(5*tof*day_scale, length(times))
        )
    
        departure_dv = lambert_departure - starting_locations[4:6, :]
        arrival_dv = ending_locations[4:6, :] - lambert_arrival
    
        departure_dv_norm = norm.(eachcol(departure_dv))
        arrival_dv_norm = norm.(eachcol(arrival_dv))
    
        total_dv_norm = departure_dv_norm + arrival_dv_norm
    
        total_dv_norm_all = min.(total_dv_norm, total_dv_norm_all)
    end

    return total_dv_norm_all
end


function get_lambert_cost(
    id_start,
    id_targets,
    time,
)
    starting_locations = repeat(asteroids_cartesian_5day[:, id_start, time], 1, size(id_targets, 1))
    total_dv_norm_all = fill(1000.0, size(starting_locations, 2))

    # for tof in 20:10:80
    for tof in 20:10:50
        ending_locations = asteroids_cartesian_5day[:, id_targets, time + tof]

        lambert_departure, lambert_arrival = zero_revolution_lambert(
            starting_locations[1:3, :],
            ending_locations[1:3, :],
            fill(5*tof*day_scale, length(ending_locations))
        )

        departure_dv = lambert_departure - starting_locations[4:6, :]
        arrival_dv = ending_locations[4:6, :] - lambert_arrival
    
        departure_dv_norm = norm.(eachcol(departure_dv))
        arrival_dv_norm = norm.(eachcol(arrival_dv))
    
        total_dv_norm = departure_dv_norm + arrival_dv_norm
    
        total_dv_norm_all = min.(total_dv_norm, total_dv_norm_all)
    end

    return total_dv_norm_all
end

function get_periodic_transfer_cost(
    start_id
)
    times = vcat(
        collect(convert_mjd_to_5day_time(65000):10:convert_mjd_to_5day_time(66500)),
        collect(convert_mjd_to_5day_time(67800):10:convert_mjd_to_5day_time(69300))
    )

    transfer_dvs = [v_scale.*get_lambert_phase(start_id, i, times) for i in pruned_ids]


    # temp = sum.(good_transfer_dvs)./length.(good_transfer_dvs) .<= 3.0

    return sum.(transfer_dvs)./length.(transfer_dvs)
end

function get_periodic_transfer_cost2(
    start_id
)
    times = vcat(
        collect(convert_mjd_to_5day_time(65000):10:convert_mjd_to_5day_time(66500)),
        collect(convert_mjd_to_5day_time(67800):10:convert_mjd_to_5day_time(69300))
    )

    transfer_dvs = [v_scale.*get_lambert_phase(start_id, i, times) for i in pruned_ids]


    # temp = sum.(good_transfer_dvs)./length.(good_transfer_dvs) .<= 3.0

    # [max(minimum(dvs[1:(length(dvs) ÷ 2)]), minimum(dvs[(length(dvs) ÷ 2 + 1):end])) for dvs in transfer_dvs]

    return [max(minimum(dvs[1:(length(dvs) ÷ 2)]), minimum(dvs[(length(dvs) ÷ 2 + 1):end])) for dvs in transfer_dvs]
end

function get_close_periodic_transfers(
    start_id,
)
    removed_ids = Int64[]

    for i in 1:5
        temp = get_lambert_cost(
            start_id,
            pruned_ids,
            convert_mjd_to_5day_time(mjd_start + i*800),
        )

        removed_ids = vcat(removed_ids, pruned_ids[temp .>= 10.0 / v_scale])
    end

    transfer_ids = setdiff(pruned_ids, removed_ids)

    times = collect(convert_mjd_to_5day_time(mjd_start):convert_mjd_to_5day_time(mjd_end - 400))

    times = vcat(
        collect(convert_mjd_to_5day_time(65000):convert_mjd_to_5day_time(66500)),
        collect(convert_mjd_to_5day_time(67800):convert_mjd_to_5day_time(69300))
    )

    # start_period = 2π * asteroids_classical[1, start_id]^1.5
    # transfer_periods = 2π * asteroids_classical[1, transfer_ids].^1.5
    # period_ratio = start_period./transfer_periods

    # good_transfer_ids = transfer_ids[abs.(1.0 .- period_ratio) .< 0.1]
    # good_transfer_ids = transfer_ids[abs.(1.0 .- period_ratio) .< 0.025]
    good_transfer_ids = transfer_ids
    good_transfer_dvs = [v_scale.*get_lambert_phase(start_id, i, times) for i in good_transfer_ids]


    # temp = sum.(good_transfer_dvs[1:length(times[1])])

    temp1 = [minimum(i[1:(length(times) ÷ 2)]) for i in good_transfer_dvs]
    temp2 = [minimum(i[(length(times) ÷ 2 + 1):end]) for i in good_transfer_dvs]


    # temp = (temp1 .<= 1.25) .& (temp2 .<= 1.25)
    # temp = (temp1 .<= 1.0) .& (temp2 .<= 1.0)
    temp = (temp1 .<= 1.5) .& (temp2 .<= 1.5)
    # temp .|= sum.(good_transfer_dvs)./length.(good_transfer_dvs) .<= 2.5
    temp .|= sum.(good_transfer_dvs)./length.(good_transfer_dvs) .<= 2.5
    # temp .&= sum.(good_transfer_dvs)./length.(good_transfer_dvs) .<= 4.0

    temp = sum.(good_transfer_dvs)./length.(good_transfer_dvs) .<= 3.0

    good_transfer_ids = good_transfer_ids[temp]
    good_transfer_dvs = good_transfer_dvs[temp]

    return good_transfer_ids, good_transfer_dvs
end


function get_lambert_cost_over_time(
    start_id,
    good_transfer_ids
)
    times = collect(convert_mjd_to_5day_time(mjd_start):convert_mjd_to_5day_time(mjd_end - 400))

    good_transfer_dvs = [v_scale.*get_lambert_phase(start_id, i, times) for i in good_transfer_ids]

    return good_transfer_dvs
end



function get_periodic_transfers(
    start_id,
    chosen_transfer_ids,
)
    times = collect(convert_mjd_to_5day_time(mjd_start):convert_mjd_to_5day_time(mjd_end - 400))

    good_transfer_ids = chosen_transfer_ids
    good_transfer_dvs = [v_scale.*get_lambert_phase(start_id, i, times) for i in good_transfer_ids]

    return good_transfer_ids, good_transfer_dvs
end

function get_best_earth_transfer_times(
    subset_ids;
    print_information = true
)
    best_start = 100.0
    best_end = 0.0

    best_start_id = 0
    best_end_id = 0

    for start_id in subset_ids
        temp_start, start_delta = find_best_transfer_time_forward(0, start_id, 0.0)
        temp_end, end_delta = find_best_transfer_time_reverse(start_id, -3, convert_mjd_to_time(mjd_end))
            
        if temp_start < best_start
            best_start_id = start_id
        end

        if temp_end > best_end
            best_end_id = start_id
        end

        best_start = min(best_start, temp_start)
        best_end = max(best_start, temp_end)
    end

    start_time, start_delta = find_best_transfer_time_forward(0, best_start_id, 0.0)
    end_time, end_delta = find_best_transfer_time_reverse(best_end_id, -3, convert_mjd_to_time(mjd_end))

    if print_information
        println("\nStart time: $(start_time)\nEnd time: $(end_time)")
    end

    return [start_time, start_time + start_delta, end_time - end_delta, end_time]
end

function get_control_for_all_plotting(
    id_journey_all,
    times_journey_all
)
    global node_time_spacing = 5.0 * day_scale
    global scp_iterations = 10
    global mass_starting = 3000.0 / m_scale
    global reverse_mass = false

    t_nodes_all = []
    u_nodes_all = []
    x_initial_departure_cartesian_all = []

    for (id_journey, times_journey) in zip(id_journey_all, times_journey_all)
        global segments = length(id_journey) - 1

        mass_changes = get_mass_change_at_ids_mixed(
            id_journey,
            times_journey,
            id_journey_all,
            times_journey_all
        )
    
        locations_journey, Δv_journey_combined, Δv_journey_departure, Δv_journey_arrival = get_lambert_trajectory(
            times_journey,
            id_journey
        )
    
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

        t_nodes, u_nodes = process_output_for_plotting(times_journey, t_nodes, u_nodes)

        x_initial_departure_cartesian = locations_journey[:, 1]
        x_initial_departure_cartesian[4:6] .+= Δv_departure_injection[1]
        x_initial_departure_cartesian = vcat(x_initial_departure_cartesian, mass_starting)
        
        push!(t_nodes_all, t_nodes)
        push!(u_nodes_all, u_nodes)
        push!(x_initial_departure_cartesian_all, x_initial_departure_cartesian)
    end

    return t_nodes_all, u_nodes_all, x_initial_departure_cartesian_all
end

function get_ids_and_times_from_mat(filename)
    file = matopen(filename)

    id_journey = round.(Int64, vcat(read(file, "MSD")["astID"]', read(file, "MSC")["astID"]'))[:]
    times_journey = convert_mjd_to_time(vcat(read(file, "MSD")["tt"]', read(file, "MSC")["tt"]')[:])

    id_journey[end] = -3

    return id_journey, times_journey
end

function choose_best_ships(
    id_journey_all,
    penalised_mass_all,
    mined_mass_all,
    groups_all;
    max_ships = 40
)
    ship_number = length(id_journey_all)

    # model = Model(HiGHS.Optimizer)
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    # @variable(model, 0 <= id[1:60000] <= 2, Int)
    @variable(model, ship[1:ship_number], Bin)

    for group in unique(groups_all)
        group_indices = findall(==(group), groups_all)

        @constraint(model, 
            [i=1:length(group_indices)],
            sum(ship[group_indices]) == length(group_indices)*ship[group_indices[i]]
        )
    end


    for asteroid_id in 1:60000
        count_for_id = [count(==(asteroid_id), id_journey_all[ship_id]) for ship_id in 1:ship_number]

        if sum(count_for_id) != 0
            @constraint(model,
                sum(count_for_id .* ship) <= 2
            )
        end
    end

    # Resolve problem for each number of ships

    chosen_ships_all = []
    average_mined_mass_all = []
    average_penalised_mass_all = []
    allowed_mined_mass_all = []

    for allowed_number in 1:max_ships
        maximum_constraint = @constraint(model,
            sum(ship) <= allowed_number
        )

        push!(allowed_mined_mass_all, max(0.0, log(allowed_number/2)/0.004))

        mass_constraint = @constraint(model,
            sum(mined_mass_all .* ship) / allowed_number >= allowed_mined_mass_all[end]
        )

        @objective(model, Max, sum(penalised_mass_all .* ship))
        JuMP.optimize!(model)

        status = termination_status(model)

        # println(status)

        # The solution is not feasible, so find the best we can do
        if status != OPTIMAL
            delete(model, mass_constraint)
            @objective(model, Max, sum(mined_mass_all .* ship))
            JuMP.optimize!(model)
        end

        chosen_ships = (1:ship_number)[Bool.(round.(value.(ship)))]

        push!(chosen_ships_all, chosen_ships)
        push!(average_mined_mass_all, sum(mined_mass_all[chosen_ships]) / allowed_number)
        push!(average_penalised_mass_all, sum(penalised_mass_all[chosen_ships]) / allowed_number)

        delete(model, maximum_constraint)

        if status == OPTIMAL
            delete(model, mass_constraint)
        end
    end

    return chosen_ships_all, average_mined_mass_all, average_penalised_mass_all, allowed_mined_mass_all
end

function generate_default_times(
    journey_length, 
    mixing_number,
    subset_ids;
    algorithm = 0,
    time_parameter_days = 165
)
    best_transfer_times = get_best_earth_transfer_times(subset_ids)

    if algorithm == 0
        # 'speading'

        times_deploy = best_transfer_times[2] .+ 1.2*cumsum([0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]*time_parameter_days*day_scale)[1:journey_length]

        times_collect = reverse(best_transfer_times[end-1] .- cumsum([0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]*time_parameter_days*day_scale)[1:journey_length])
    elseif algorithm == 1
        # 'fixed'

        times_deploy = best_transfer_times[2] .+ 1.2*cumsum([0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]*time_parameter_days*day_scale)[1:journey_length]

        times_collect = reverse(best_transfer_times[end-1] .- cumsum([0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]*time_parameter_days*day_scale)[1:journey_length])
    elseif algorithm == 2
        # 'random'

        times_deploy = best_transfer_times[2] .+ 1.2*cumsum(vcat([0.0], 0.8 .+ 0.5*rand(10))*time_parameter_days*day_scale)[1:journey_length]

        times_collect = reverse(best_transfer_times[end-1] .- cumsum(vcat([0.0], 0.8 .+ 0.5*rand(10))*time_parameter_days*day_scale)[1:journey_length])
    elseif algorithm == 3
        # Good testing

        times_deploy = best_transfer_times[2] .+ time_parameter_days*day_scale*cumsum([0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0])[1:journey_length]

        times_collect = reverse(best_transfer_times[end-1] .- time_parameter_days*day_scale*cumsum([0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0])[1:journey_length])
    end

    return [
        vcat(
            [best_transfer_times[1]],
            times_deploy,
            times_collect, 
            [best_transfer_times[end]]
        ) for _ in 1:mixing_number
    ]
end

function get_deploy_collect_lengths(
    id_journey_mixed,
    times_journey_mixed
)
    deploy_lengths = []
    collect_lengths = []
    
    for i in 1:length(times_journey_mixed)
        int_mass_changes = get_mass_change_at_ids_mixed(
            id_journey_mixed[i],
            times_journey_mixed[i],
            id_journey_mixed,
            times_journey_mixed,
        )[2:end-1]
    
        push!(deploy_lengths, length(int_mass_changes[int_mass_changes .< 0]))
        push!(collect_lengths, length(int_mass_changes[int_mass_changes .> 0]))
    end
    
    return deploy_lengths, collect_lengths
end

function extract_result_statistics(result_file)
    id_journey_all, times_journey_all, mined_mass_all, penalised_mass_all, low_thrust_dv_all = load_result_files([result_file]; low_thrust_dv = true);

    initial_transfer_time_all = []
    final_transfer_time_all = []
    intermediate_transfer_time_all = []
    middle_transfer_time_all = []

    initial_transfer_lambert_dv_all = []
    final_transfer_lambert_dv_all = []

    initial_transfer_low_thrust_dv_all = []
    final_transfer_low_thrust_dv_all = []
    intermediate_transfer_low_thrust_dv_all = []
    middle_transfer_low_thrust_dv_all = []

    for (n, (id_journey, times_journey)) in enumerate(zip(id_journey_all, times_journey_all))
        mass_changes = get_mass_change_at_ids_mixed(
            id_journey,
            times_journey,
            id_journey_all,
            times_journey_all
        )

        push!(initial_transfer_time_all, times_journey[2] - times_journey[1])
        push!(final_transfer_time_all, times_journey[end] - times_journey[end - 1])

        middle_transfer_index = findfirst(>(0.0), mass_changes) - 1

        temp_transfer_time = []

        for i in 2:(length(times_journey) - 2)
            # push!(intermediate_transfer_time_all, times_journey[i + 1] - times_journey[i])
            push!(temp_transfer_time, times_journey[i + 1] - times_journey[i])

            if i == middle_transfer_index
                push!(middle_transfer_time_all, times_journey[i + 1] - times_journey[i])
            else
                # push!(intermediate_transfer_time_all, times_journey[i + 1] - times_journey[i])
                # if (times_journey[i + 1] - times_journey[i])/day_scale < 10000
                    
                #     push!(intermediate_transfer_time_all, times_journey[i + 1] - times_journey[i])
                # else
                #     println("asd")
                # end

                # push!(temp_transfer_time, times_journey[i + 1] - times_journey[i])
            end
        end

        new_sort = sortperm(temp_transfer_time)

        temp_transfer_time = temp_transfer_time[new_sort][1:(end-1)]
        intermediate_transfer_time_all = vcat(intermediate_transfer_time_all, temp_transfer_time)


        push!(initial_transfer_lambert_dv_all, 
            get_transfer_dv(id_journey[1], times_journey[1], id_journey[2], times_journey[2] - times_journey[1])[1]
        )

        push!(final_transfer_lambert_dv_all, 
            get_transfer_dv(id_journey[end-1], times_journey[end-1], id_journey[end], times_journey[end] - times_journey[end-1])[1]
        )

        push!(initial_transfer_low_thrust_dv_all, low_thrust_dv_all[n][1])
        push!(final_transfer_low_thrust_dv_all, low_thrust_dv_all[n][end])

        temp_low_thrust_dv = []

        for i in 2:(length(low_thrust_dv_all[n]) - 1)
            # push!(intermediate_transfer_low_thrust_dv_all, low_thrust_dv_all[n][i])
            push!(temp_low_thrust_dv, low_thrust_dv_all[n][i])

            if i == middle_transfer_index
                push!(middle_transfer_low_thrust_dv_all, low_thrust_dv_all[n][i])
            else
                # push!(temp_low_thrust_dv, low_thrust_dv_all[n][i])

                # if (times_journey[i + 1] - times_journey[i])/day_scale < 800
                #     push!(intermediate_transfer_low_thrust_dv_all, low_thrust_dv_all[n][i])
                # end


                # push!(intermediate_transfer_low_thrust_dv_all, low_thrust_dv_all[n][i])
            end
        end

        temp_low_thrust_dv = temp_low_thrust_dv[new_sort][1:(end-1)]
        intermediate_transfer_low_thrust_dv_all = vcat(intermediate_transfer_low_thrust_dv_all, temp_low_thrust_dv)
    end



    return [
        split(split(result_file, "/")[end], ".")[1],
        mean(initial_transfer_time_all)/day_scale,
        mean(final_transfer_time_all)/day_scale,
        mean(middle_transfer_time_all)/day_scale,
        mean(intermediate_transfer_time_all)/day_scale,
        v_scale*mean(initial_transfer_low_thrust_dv_all),
        v_scale*mean(final_transfer_low_thrust_dv_all),
        v_scale*mean(middle_transfer_low_thrust_dv_all),
        v_scale*mean(intermediate_transfer_low_thrust_dv_all),
        sum(mined_mass_all)/length(mined_mass_all)
    ]
end

function ephemeris_classical_at(classical, times; μ=1.0)
    n = sqrt.(μ./classical[1,:].^3)
    
    out = repeat(classical, 1, 1, size(times, 1))

    out[6,:,:] .+= n.*times'
    out[6,:,:] .%= 2π

    return out
end

function ephemeris_cartesian_at(classical, times; μ=1.0)
    out_classical = ephemeris_classical_at(classical, times; μ)

    out = stack(classical_to_cartesian.(eachslice(out_classical; dims=3); μ))

    return out
end