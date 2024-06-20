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


































































# Self cleaning
# id_journey_mixed = [
#     [0 35457 22813 45607 25996 21927 45750 19821 19821 22813 35457 25996 45750 45607 21927 -3][:],
#     [0 14875 40259 665 34568 59701 37913 39279 39279 37913 40259 665 34568 59701 14875 -3][:],
#     [0 595 30058 59481 10976 21375 34927 38331 38331 10976 21375 30058 595 59481 34927 -3][:],
# ]

# # times_journey_mixed = [
#     convert_mjd_to_time([64428.0 65013.0 65258.0 65458.0 65718.0 65948.0 66118.0 66313.0 67843.0 68083.0 68378.0 68593.0 68793.0 69078.0 69313.0 69783.0][:]),
#     convert_mjd_to_time([64448.0 65008.0 65248.0 65468.0 65778.0 65933.0 66118.0 66248.0 67813.0 68053.0 68278.0 68708.0 68893.0 69038.0 69343.0 69783.0][:]),
#     convert_mjd_to_time([64448.0 65023.0 65128.0 65253.0 65463.0 65643.0 65848.0 66093.0 67983.0 68233.0 68483.0 68653.0 68863.0 69138.0 69348.0 69763.0][:]),
# ]

# # Mixed 
# id_journey_mixed = [
#     [0 35457 22813 45607 25996 21927 45750 19821 19821 21375 59481 22813 45607 665 14875 -3][:],
#     [0 14875 40259 665 34568 59701 37913 39279 39279 595 10976 35457 34568 59701 34927 -3][:],
#     [0 595 30058 59481 10976 21375 34927 38331 38331 45750 40259 30058 37913 21927 25996 -3][:],
# ]

# times_journey_mixed = [
#     convert_mjd_to_time([64428.0 65013.0 65258.0 65458.0 65718.0 65948.0 66118.0 66313.0 68168.0 68408.0 68538.0 68748.0 68888.0 69168.0 69343.0 69803.0][:]),
#     convert_mjd_to_time([64448.0 65008.0 65248.0 65468.0 65778.0 65933.0 66118.0 66248.0 67973.0 68213.0 68558.0 68823.0 69033.0 69148.0 69353.0 69803.0][:]),
#     convert_mjd_to_time([64448.0 65023.0 65128.0 65253.0 65463.0 65643.0 65848.0 66093.0 67478.0 67813.0 68358.0 68713.0 68963.0 69153.0 69328.0 69803.0][:]),
# ]







id_journey = [
    [0, 56608, 56608, -3]
]

times_journey = [
    [1.0321259370000504, 10.827001079130554, 85.01277301090416, 92.92573852790454]
]




id_journey = [
    # [0 35457 22813 22813 35457][:],
    # [0 35457 22813 22813 35457 -3][:],
    # [0 35457 22813 45607 25996 21927 45750 19821 19821 22813 35457 25996 45750 45607 21927 -3][:],
    # [0 14875 40259 665 34568 59701 37913 39279 39279 37913 40259 665 34568 59701 14875 -3][:],

    [0 35457 22813 45607 25996 21927 45750 19821 19821 21375 59481 22813 45607 665 14875 -3][:],
    [0 14875 40259 665 34568 59701 37913 39279 39279 595 10976 35457 34568 59701 34927 -3][:],
    [0 595 30058 59481 10976 21375 34927 38331 38331 45750 40259 30058 37913 21927 25996 -3][:],
]

times_journey = [
    # convert_mjd_to_time([64428.0 65013.0 65258.0 69078.0 69313.0][:]),
    # convert_mjd_to_time([64448.0 65008.0 65248.0 65468.0 65778.0 65933.0 66118.0 66248.0 67813.0 68053.0 68278.0 68708.0 68893.0 69038.0 69343.0 69783.0][:]),
    # convert_mjd_to_time([64448.0 65008.0 65248.0 65468.0 65778.0 65933.0 66118.0 66248.0 67813.0 68053.0 68278.0 68708.0 68893.0 69038.0 69343.0 69783.0][:]),

    convert_mjd_to_time([64428.0 65013.0 65258.0 65458.0 65718.0 65948.0 66118.0 66313.0 68168.0 68408.0 68538.0 68748.0 68888.0 69168.0 69343.0 69803.0][:]),
    convert_mjd_to_time([64448.0 65008.0 65248.0 65468.0 65778.0 65933.0 66118.0 66248.0 67973.0 68213.0 68558.0 68823.0 69033.0 69148.0 69353.0 69803.0][:]),
    convert_mjd_to_time([64448.0 65023.0 65128.0 65253.0 65463.0 65643.0 65848.0 66093.0 67478.0 67813.0 68358.0 68713.0 68963.0 69153.0 69328.0 69803.0][:]),
]


scp_iterations = 80

node_time_spacing = 20.0*day_scale





id_journey, times_journey = load_result_files(
    ["data/PART_1.txt"]
)





p = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.1,
    mass_overhead = 1.0/m_scale
);

solve!(p,
    MixedTimeAdaptive(); 
    adaptive_time = true
)

convert_logged_mass_to_mass!(p)

solve!(p,
    MixedTimeAdaptive(); 
    adaptive_time = true
)

solve!(p,
    MixedTimeAdaptive(); 
    adaptive_time = false
)



plot_thrust_information(p)

plot_trajectory(p)










write_solution(p, "output/Result.txt")


for i in 1:length(p.u_nodes[1])
    display(maximum(norm.(norm.(eachcol(p.u_nodes[1][i][1:3, :])) .- p.u_nodes[1][i][4, :])))
end

