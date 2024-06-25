include("header.jl")

GLMakie.activate!()
CairoMakie.activate!()



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
    trust_region_factor = 0.05,
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



plot_trajectory(p; plot_3d = true)




write_solution(p, "output/Result.txt")


for i in 1:length(p.u_nodes[1])
    display(maximum(norm.(norm.(eachcol(p.u_nodes[1][i][1:3, :])) .- p.u_nodes[1][i][4, :])))
end














id_subset = [298 3889 4445 7103 10290 10649 10916 12577 14062 14079 15291 15906 16070 16110 18835 18913 19496 20616 20651 22313 24024 25663 27293 27414 28050 30492 31751 32796 34679 36517 36565 36716 38270 38514 41013 41509 41629 43596 46392 46789 46877 47175 47291 47647 49082 49251 50836 50898 53789 54572 55377 56717 57669 57998 59572 59916][:]

id_subset = sort([15184, 3241, 2032, 53592, 46418, 19702, 23056, 46751, 32088, 23987][:])





mip_problem = MixedIntegerProblem(id_subset, [5], [5])
mip_problem.cost_limit = 6/v_scale



solve!(mip_problem;
    # self_cleaning = true,
    include_intermediate_transfer_cost = true,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 10
)

plot_graph(mip_problem)




scp_problem = SequentialConvexProblem(
    mip_problem.id_journey_solutions[1], 
    mip_problem.times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.05,
    mass_overhead = 1.0/m_scale
);


solve!(scp_problem,
    MixedTimeAdaptive(); 
    adaptive_time = true
)

convert_logged_mass_to_mass!(scp_problem)

solve!(scp_problem,
    MixedTimeAdaptive(); 
    adaptive_time = true
)

scp_problem.mass_overhead = 0.1/m_scale

solve!(scp_problem,
    MixedTimeAdaptive(); 
    adaptive_time = false
)

plot_thrust_information(scp_problem)

plot_trajectory(scp_problem)





mip_problem = MixedIntegerProblem(id_subset, [8], [8];
    times_journey = scp_problem.times_journey
)



write_solution(scp_problem, "output/Result.txt")



for i in 1:length(scp_problem.u_nodes[1])
    display(maximum(norm.(norm.(eachcol(scp_problem.u_nodes[1][i][1:3, :])) .- scp_problem.u_nodes[1][i][4, :])))
end






