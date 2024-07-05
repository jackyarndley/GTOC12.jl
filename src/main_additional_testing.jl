include("header.jl")


g0 = 9.80665e-3/(r_scale/t_scale^2) # LU/TU^2
g0_isp = 19.6133*t_scale/r_scale
isp = g0_isp/g0


thrust_scale = 1.0
thrust = thrust_scale*0.5e-3/(m_scale*r_scale/t_scale^2)



custom_cartesian = [
    -140699693 -51614428 980 9.7746 -28.0783 4.3377e-4
    -172682023 176959469 7948912 -16.4274 -14.8605 9.2149e-2
]'

# custom_cartesian = [
#     -3.01166620e+07 -1.48851968e+08  1.86741232e+02 2.87120370e+01 -6.01882549e+00  7.55087690e-06
#     -1.83865882e+06  2.34962319e+08  4.96736048e+06 -23.3085395 1.86654249 0.61187782
# ]'


custom_cartesian[1:3, :] ./= r_scale
custom_cartesian[4:6, :] ./= v_scale

custom_classical = stack(cartesian_to_classical.(
    eachcol(custom_cartesian)
))












scp_iterations = 20

node_time_spacing = 5.0*day_scale

node_time_spacing = 348.79*day_scale/40
# node_time_spacing = 689.5745203452379*day_scale/40
# node_time_spacing = 689.5745203452379*day_scale/8



id_journey = [
    [60001, 60002]
]

times_journey = [
    [0.0, 348.79*day_scale]
    # [0.0, 689.5745203452379*day_scale]
]




p = SequentialConvexProblem(
    id_journey, 
    times_journey;
    objective_config = LoggedMassConfig(),
    trust_region_factor = 0.01
);


plot_trajectory(p)

solve!(p,
    MixedTimeAdaptive(); 
    adaptive_time = false
)




convert_logged_mass_to_mass!(p)

solve!(p,
    MixedTimeAdaptive(); 
    adaptive_time = false
)





plot_thrust_information(p)

plot_trajectory(p; plot_3d = true)




for i in 1:length(p.u_nodes[1])
    display(maximum(norm.(norm.(eachcol(p.u_nodes[1][i][1:3, :])) .- p.u_nodes[1][i][4, :])))
end




