include("header.jl")

GLMakie.activate!()
CairoMakie.activate!()

using LaTeXStrings
using Latexify
using DataFrames
using Format

Latexify.set_default(; env=:tabular)

scp_iterations = 80

node_time_spacing = 20.0*day_scale


id_journey = [
    [0, 15184, 12286, 15184, 12286, -3]
]

times_journey = [
    convert_mjd_to_time([64438.0 64913.0 65208.0 65648.0 66145.0 66440.0][:]),
]

p = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.01,
    mass_overhead = 1.0/m_scale
);

solve!(p)

convert_logged_mass_to_mass!(p)

solve!(p)

solve!(p; fixed_rendezvous = true)


plot_trajectory(p)

plot_gtoc12_problem(p)

plot_scp_details(p)

plot_thrust_information(p)


id_subset = [3241, 15184, 19702, 46418, 53592]
id_subset = sort([2032, 3241, 15184, 19702, 23056, 23987, 32088, 46418, 46751, 53592, 3896, 37818, 15083, 5707, 19434, 981, 48748, 40804, 23483, 47817])


# Create DataFrame
df = DataFrame(
    [string.(id_subset),
    asteroids_classical[1, id_subset],
    asteroids_classical[2, id_subset],
    asteroids_classical[3, id_subset],
    asteroids_classical[4, id_subset],
    asteroids_classical[5, id_subset],
    asteroids_classical[6, id_subset]],
    [:ID, :a, :e, :i, :Ω, :ω, :M]
)


df[:, :i] = rad2deg.(Vector(df[:, :i]))
df[:, :Ω] = rad2deg.(Vector(df[:, :Ω]))
df[:, :ω] = rad2deg.(Vector(df[:, :ω]))
df[:, :M] = rad2deg.(Vector(df[:, :M]))


print(latexify(df; fmt = "%.3f"))





mip_problem = MixedIntegerProblem(
    id_subset, 
    [3], 
    [3],
    time_parameter_days = 175
)
mip_problem.cost_limit = 6/v_scale

solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 3,
)

times_journey = mip_problem.times_journey
id_journey = mip_problem.id_journey_solutions[1]


df = DataFrame(
    convert_time_to_mjd.(times_journey[1])',
    ["Earth Departure", "Deployment 1", "Deployment 2", "Deployment 3", "Collection 1", "Collection 2", "Collection 3", "Earth Arrival"]
)

print(latexify(df; fmt = "%.2f", side = [L"Time [MJD]"]))




p1 = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.02,
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

p2.trust_region_factor = 0.01

solve!(p2; 
    fixed_segments = false,
    fixed_rendezvous = true
)

p3 = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.2,
    mass_overhead = 1.0/m_scale
);

p3.trust_region_factor = 0.05

solve!(p3; 
    fixed_segments = false,
    fixed_rendezvous = false
)


plot_trajectory_paper_2(
    p3;
    # label_text = "Ship 1:\n10 deployments\n10 collections\n781.41 kg returned",
    output_file = "output/plots/scp_example_individual.png"
)


temp = [64328.00,64848.95,64952.82,65137.31,69109.33,69256.70,69413.11,69791.29]


plot_discretization_comparison(p1, p2, p3; output_file = "output/plots/discretization_comparison.png")








id_subset = [3241, 15184, 19702, 46418, 53592]


mip_problem = MixedIntegerProblem(id_subset, [3], [3]; time_parameter_days = 175)


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


mip_problem = MixedIntegerProblem(id_subset, [3, 2], [2, 3]; time_parameter_days = 175)
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






id_subset = [3241, 15184, 19702, 46418, 53592]


mip_problem1 = MixedIntegerProblem(id_subset, [3], [3]; time_parameter_days = 175)


mip_problem1.cost_limit = 6/v_scale


solve!(mip_problem1;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 1000.0,
    solutions_count_maximum = 10000
)


mip_problem2 = MixedIntegerProblem(id_subset, [3], [3]; time_parameter_days = 175)


mip_problem2.cost_limit = 60000/v_scale


solve!(mip_problem2;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 1000.0,
    solutions_count_maximum = 10000
)



plot_bip_solution_values(
    mip_problem1, 
    mip_problem2; 
    output_file = "output/plots/bip_solution_values.png"
)






# id_subset = sort([15184, 3241, 2032, 53592, 46418, 19702, 23056, 46751, 32088, 23987])
id_subset = sort([2032, 3241, 15184, 19702, 23056, 23987, 32088, 46418, 46751, 53592, 3896, 37818, 15083, 5707, 19434, 981, 48748, 40804, 23483, 47817, 2174, 28289, 43836, 39557, 9260, 17983, 13655, 22108, 3302, 57913])
# id_subset = sort([2032, 3241, 15184, 19702, 23056, 23987, 32088, 46418, 46751, 53592, 3896, 37818, 15083, 5707, 19434, 981, 48748, 40804, 23483, 47817])


mip_problem = MixedIntegerProblem(id_subset, [10], [10]; time_parameter_days = 145)
mip_problem.cost_limit = 100/v_scale
# mip_problem.cost_limit = 10/v_scale
# mip_problem.cost_limit = 8/v_scale



# join([@sprintf("%5s ", val) for val in id_subset])

# join([@sprintf("%5i ", val) for val in convert_time_to_mjd.(mip_problem.times_journey[1])])

# join([@sprintf("%5s ", val) for val in mip_problem.id_journey_solutions[5][1]])


# join([@sprintf("%s ", val) for val in round.(convert_time_to_mjd.(scp_problem.times_journey[1]); digits=2)])



mip_problem_objectives = Vector{Float64}[]
scp_problem_objectives = Vector{Float64}[]





solution_number = 1


solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.2,
    solutions_count_maximum = 1,
    # solutions_count_maximum = 3*solution_number,
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

@time solve!(scp_problem)



push!(scp_problem_objectives, [-m_scale*scp_problem.Δm0[n][end] for n in 1:scp_problem.mixing_number])




plot_convergence_comparison(
    [collect(1:min(solution_number, length(mip_problem_objectives[i]))) for i in 1:length(mip_problem_objectives)],
    scp_problem_objectives;
    output_file = "output/plots/nested_loop_progress.png"
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



solution_number = 50


mip_problem = MixedIntegerProblem(id_subset, [10], [10]; time_parameter_days = 145)
mip_problem.cost_limit = 6/v_scale


solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.2,
    solutions_count_maximum = 3*solution_number,
    time_limit_seconds = 300
)


mip_problem.solutions = 50


plot_graph_structure_big(
    mip_problem;
    plot_pruned = true,
    plot_optimal_path = true,
    output_file = "output/plots/bip_big.png",
    figure_size = (800, 400)
)





scp_problem = SequentialConvexProblem(
    [mip_problem.id_journey_solutions[1][1]], 
    [mip_problem.times_journey[1]];
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.01,
    mass_overhead = 0.0/m_scale,
    dynamical_error = 1e-4
);

scp_problem.trust_region_factor = 0.005

solve!(scp_problem)


plot_thrust_information(scp_problem; solution_indices = [1])

plot_trajectory(scp_problem; solution_indices = [1], plot_3d = false, rotating = true)






times_5day = collect(0:5*day_scale:maximum_time)
asteroids_cartesian_5day = ephemeris_cartesian_at(asteroids_classical, times_5day)



check_times_deploy = collect(convert_mjd_to_5day_time(64900):20:convert_mjd_to_5day_time(66500))
check_times_collect = collect(convert_mjd_to_5day_time(68200):20:convert_mjd_to_5day_time(69000))


times_join_check = collect(convert_mjd_to_time(65000):100*day_scale:convert_mjd_to_time(69000))
asteroids_join_check = ephemeris_cartesian_at(asteroids_classical, times_join_check)

id_subset = sort([2032, 3241, 15184, 19702, 23056, 23987, 32088, 46418, 46751, 53592, 3896, 37818, 15083, 5707, 19434, 981, 48748, 40804, 23483, 47817, 2174, 28289, 43836, 39557, 9260, 17983, 13655, 22108, 3302, 57913])


while length(id_subset) < 100
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











scp_iterations = 80

node_time_spacing = 20.0*day_scale


id_journey, times_journey = load_result_files(
    ["data/PART_1.txt"]
)

p = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.02,
    mass_overhead = 0.1/m_scale
);

solve!(p)


plot_thrust_profile_paper(p; output_file = "output/plots/individual_profile.png")



# plot_trajectory_paper(p; rotating = false, plot_3d = false)

plot_trajectory_paper(p; rotating = true, plot_3d = false, output_file = "output/plots/individual_trajectory.png")

# plot_trajectory_paper(p; rotating = true, plot_3d = true)

plot_trajectory_paper_2(
    p;
    output_file = "output/plots/individual_trajectory.png"
)

# plot_trajectory_and_thrust_profile_paper(
#     p;
#     label_text = "Ship 1:\n10 deployments\n10 collections\n781.41 kg returned",
#     output_file = "output/plots/individual_trajectory.png"
# )








id_subset = sort([2032, 3241, 15184, 19702, 23056, 23987, 32088, 46418, 46751, 53592, 3896, 37818, 15083, 5707, 19434, 981, 48748, 40804, 23483, 47817, 2174, 28289, 43836, 39557, 9260, 17983, 13655, 22108, 3302, 57913])

mip_problem = MixedIntegerProblem(
    id_subset, 
    [9, 8], 
    [8, 9],
    time_parameter_days = 145
)

mip_problem.cost_limit = 6/v_scale
# mip_problem.cost_limit = 8/v_scale



solution_number = 25

solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.2,
    solutions_count_maximum = 3*solution_number,
    time_limit_seconds = 20*60
)


push!(mip_problem_objectives, mip_problem.objective_solutions.*v_scale)

scp_problem = SequentialConvexProblem(
    # [mip_problem.id_journey_solutions[k][:] for k in 1:min(solution_number, mip_problem.solutions)], 
    # [mip_problem.times_journey[1] for k in 1:min(solution_number, mip_problem.solutions)];
    mip_problem.id_journey_solutions[1][:], 
    mip_problem.times_journey[:];
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.05,
    mass_overhead = 0.0/m_scale,
    dynamical_error = 1e-4
);

scp_problem.trust_region_factor = 0.02

@time solve!(scp_problem)



# plot_trajectory_and_thrust_profile_paper(
#     scp_problem;
#     # label_text = "Ship 1:\n10 deployments\n10 collections\n781.41 kg returned",
#     output_file = "output/plots/mixed_trajectory.png"
# )

plot_trajectory_paper_2(
    scp_problem;
    # label_text = "Ship 1:\n10 deployments\n10 collections\n781.41 kg returned",
    output_file = "output/plots/mixed_trajectory.png"
)






mip_problem = MixedIntegerProblem(id_subset, [9, 8], [8, 9];
    times_journey = scp_problem.times_journey
)

mip_problem.cost_limit = 6/v_scale







id_journey, times_journey, _, _, _, _ = load_result_folders_grouping(["data/bundled/submitted/JPL.txt"])





scp_problem = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.05,
    mass_overhead = 0.0/m_scale,
    dynamical_error = 1e-4
);



@time solve!(scp_problem)


# 791.3 seconds








id_journey1, times_journey1, _, _, _, _ = load_result_folders_grouping(["data/37_mass_optimal_self_cleaning.txt"])

id_journey2, times_journey2, _, _, _, _ = load_result_folders_grouping(["data/39_mass_optimal.txt"])


scp_iterations = 80

node_time_spacing = 20.0*day_scale


scp_problem1 = SequentialConvexProblem(
    id_journey1, 
    times_journey1;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.02,
    mass_overhead = 0.0/m_scale,
    dynamical_error = 1e-4
);

scp_problem2 = SequentialConvexProblem(
    id_journey2, 
    times_journey2;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.02,
    mass_overhead = 0.0/m_scale,
    dynamical_error = 1e-4
);


solve!(scp_problem1; fixed_rendezvous = false, fixed_segments = false)
solve!(scp_problem2; fixed_rendezvous = false, fixed_segments = false)




scp_problem1.times_journey






plot_trajectory_final(scp_problem1, scp_problem2; output_file = "output/plots/combined_both.png")



# plot_trajectory_paper(
#     scp_problem; 
#     rotating = false, 
#     # solution_indices = collect(13:16)
#     # solution_indices = collect(13:20)
# )





plot_trajectory_showcase(scp_problem2)







id_journey_all, times_journey_all, mined_mass_all, penalised_mass_all, groups_all, result_files_all = load_result_folders_grouping(
    [
        "data/bundled/reoptimized/CUSTOM",
        "data/bundled/reoptimized/CUSTOM_MIXES/JPL_MIX_1",
        "data/bundled/reoptimized/CUSTOM_MIXES/JPL_MIX_2",
        "data/bundled/reoptimized/CUSTOM_MIXES/JPL_MIX_1_REF",
        "data/bundled/reoptimized/CUSTOM_MIXES/JPL_MIX_2_REF",
        "data/bundled/reoptimized/CUSTOM_MIXES/JPL_MIX_3_REF",
        "data/bundled/reoptimized/CUSTOM_MIXES/JPL_MIX_4_REF",
        "data/bundled/reoptimized/CUSTOM_MIXES/JACK_MIX_1",
        "data/bundled/reoptimized/CUSTOM_MIXES/JACK_MIX_1_REF",
        "data/bundled/reoptimized/CUSTOM_EXTRA",
        "data/bundled/reoptimized/CUSTOM_TODO",
        "data/bundled/reoptimized/CUSTOM_REMIX",
        "data/bundled/reoptimized/CUSTOM_REF",
        "data/bundled/reoptimized/NUAA", 
        "data/bundled/reoptimized/ADL", 
        "data/bundled/reoptimized/ATQ", 
        "data/bundled/reoptimized/SIGMA", 
        "data/bundled/reoptimized/NUDT",
        "data/bundled/reoptimized/OptimiCS",
        "data/bundled/reoptimized/TheAntipodes",
        "data/bundled/reoptimized/ESA",
        "data/bundled/reoptimized/JPL",
        "data/bundled/reoptimized/JPL_POST",
        "data/bundled/reoptimized/BIT",
    ]
);

penalised_mass_all = deepcopy(mined_mass_all)


temp = is_self_cleaning.(id_journey_all)

id_journey_all = id_journey_all[temp]
penalised_mass_all = penalised_mass_all[temp]
mined_mass_all = mined_mass_all[temp]
groups_all = groups_all[temp]




@time chosen_ships_all, average_mined_mass_all, average_penalised_mass_all, allowed_mined_mass_all = choose_best_ships(id_journey_all, penalised_mass_all, mined_mass_all, groups_all; max_ships = 50)

plot_solution_ships(
    average_mined_mass_all, 
    allowed_mined_mass_all; 
    max_ships = 50, 
    max_mass = 800.0,
    output_file = "output/plots/ship_selection.png"
)












id_subset = sort([2032, 3241, 15184, 19702, 23056, 23987, 32088, 46418, 46751, 53592, 3896, 37818, 15083, 5707, 19434, 981, 48748, 40804, 23483, 47817, 2174, 28289, 43836, 39557, 9260, 17983, 13655, 22108, 3302, 57913])

mip_problem1 = MixedIntegerProblem(id_subset, [10], [10]; time_parameter_days = 145)


mip_problem1.cost_limit = 6/v_scale


solve!(mip_problem1;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 1000.0,
    solutions_count_maximum = 10000,
    time_limit_seconds = 600
)


mip_problem2 = MixedIntegerProblem(id_subset, [10], [10]; time_parameter_days = 145)


mip_problem2.cost_limit = 60000/v_scale

solve!(mip_problem2;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 1000.0,
    solutions_count_maximum = 10000,
    time_limit_seconds = 300
)



plot_bip_solution_values(
    mip_problem1, 
    mip_problem2; 
    output_file = "output/plots/bip_solution_values_2.png"
)


