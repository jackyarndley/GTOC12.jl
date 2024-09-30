include("header.jl")

GLMakie.activate!()
CairoMakie.activate!()

using PairPlots
using DataFrames
using CSV

function plot_asteroid_distribution_2d(asteroids_cartesian; title = "")
    f = Figure(size = (800, 800), backgroundcolor = :white, figure_padding = 5)

    ax = Axis(
        f[1, 1]; 
        xlabel = "x [AU]", 
        ylabel = "y [AU]",
        limits = (-4.0, 4.0, -4.0, 4.0),
        title = title
    )

    ν_plot = collect(LinRange(0.0, 2π, 400)) 

    for (i, planet_classical) in enumerate(eachcol(planets_classical))
        temp = repeat(planet_classical, 1, 400)
        temp[6, :] .= ν_plot
        temp = hcat(classical_to_cartesian.(eachcol(temp))...)

        lines!(ax,
            temp[1, :],
            temp[2, :],
            temp[3, :],
            # linestyle = :dashdot,
            color = Makie.wong_colors()[[2, 1, 6][i]],
            # alpha = 0.35,
            # label = ["Venus Orbit", "Earth Orbit", "Mars Orbit"][i]
        )
    end

    scatter!(
        ax,
        asteroids_cartesian[1, :], 
        asteroids_cartesian[2, :],
        markersize = 2,
        color = :black
    )

    display(f)

    return
end

function plot_asteroid_statistics(asteroids_classical; title = "")
    f = Figure(size = (800, 400), backgroundcolor = :white, figure_padding = 5)

    for i in 1:3
        ax = Axis(
            f[i, 1]; 
            title = i == 1 ? title : "",
            xlabel = ["a", "e", "i"][i],
            yticklabelsvisible = false,
            yticksvisible = false
        )
    
        hist!(
            ax,
            asteroids_classical[i, :], 
            bins = 100
        )
    end
    
    display(f)

    return
end



function plot_transfer_time_statistics(
    times_journey; title = ""
)
    times_first_transfer = [i[2] - i[1] for i in times_journey]/day_scale
    times_last_transfer = [i[end] - i[end-1] for i in times_journey]/day_scale

    times_intermediate_transfer = [i[3:end-1] - i[2:end-2] for i in times_journey]/day_scale
    times_intermediate_transfer_all = vcat(times_intermediate_transfer...)

    f = Figure(size = (800, 400), backgroundcolor = :white, figure_padding = 5)

    for (i, times) in zip(1:3, [times_first_transfer, times_last_transfer, times_intermediate_transfer_all])
        ax = Axis(
            f[i, 1]; 
            title = i == 1 ? title : "",
            xlabel = ["initial tof [days]", "final tof [days]", "intermediate tof [days]"][i],
            yticklabelsvisible = false,
            yticksvisible = false
        )
    
        hist!(
            ax,
            times, 
            bins = 100
        )
    end
    
    display(f)

    return
end


function plot_distribution(
    d
)
    f = Figure(size = (800, 400), backgroundcolor = :white, figure_padding = 5)

    ax = Axis(
        f[1, 1]; 
    )

    hist!(
        ax,
        rand(d, 100000), 
        bins = 100
    )

    display(f)

    return
end



function plot_rtn_state_statistics(
    rtn_start_states,
    rtn_end_states
)
    f = Figure(size = (800, 800), backgroundcolor = :white, figure_padding = 5)

    for (i, rtn_states) in enumerate([rtn_start_states, rtn_end_states])
        for j in 1:6
            temp = [k[j] for k in rtn_states]

            if (maximum(temp) - minimum(temp)) < 0.01
                continue
            end

            ax = Axis(
                f[j, i]; 
                yticklabelsvisible = false,
                yticksvisible = false
            )
        
            hist!(
                ax,
                [k[j] for k in rtn_states], 
                bins = 100
            )
        end
    end
    
    display(f)

    return
end

function plot_tof_dv_comparison(
    tof,
    dv
)
    f = Figure(size = (800, 400), backgroundcolor = :white, figure_padding = 5)

    ax = Axis(
        f[1, 1]; 
        xlabel = "TOF [days]",
        ylabel = "Lambert Δv [km/s]"
    )

    scatter!(
        ax,
        tof,
        dv
    )

    display(f)

    return
end

function plot_mass_use_statistics(
    maximum_initial_masses,
    maximum_final_masses
)
    f = Figure(size = (800, 600), backgroundcolor = :white, figure_padding = 5)

    ax1 = Axis(
        f[1, 1]; 
        xlabel = "maximum initial mass [kg]"
    )

    ax2 = Axis(
        f[2, 1]; 
        xlabel = "maximum final mass [kg]"
    )

    ax3 = Axis(
        f[3, 1]; 
        xlabel = "maximum mass use [kg]"
    )

    hist!(
        ax1,
        maximum_initial_masses, 
        bins = 100
    )
    
    hist!(
        ax2,
        maximum_final_masses, 
        bins = 100
    )

    hist!(
        ax3,
        maximum_initial_masses .- maximum_final_masses, 
        bins = 100
    )

    display(f)

    return
end






plot_asteroid_distribution_2d(asteroids_cartesian; title = "asteroid distibution no pruning")

plot_asteroid_statistics(asteroids_classical; title = "asteroid statistics no pruning")


pruned_asteroids = collect(1:60000)[(asteroids_classical[2, :] .<= 0.2) .& (2.6 .<= asteroids_classical[1, :] .<= 3.0) .& (asteroids_classical[3, :] .<= deg2rad(15))]

plot_asteroid_distribution_2d(asteroids_cartesian[:, pruned_asteroids]; title = "asteroid distibution pruning")

plot_asteroid_statistics(asteroids_classical[:, pruned_asteroids]; title = "asteroid statistics pruning")





# id_journey, times_journey, _, _, _, _ = load_result_folders_grouping(["data/39_mass_optimal.txt"])
id_journey, times_journey, _, _, _, _ = load_result_folders_grouping(["data/37_mass_optimal_self_cleaning.txt"])

selected_asteroids = setdiff(unique(vcat(id_journey...)), [0, -3])




plot_asteroid_distribution_2d(asteroids_cartesian[:, selected_asteroids]; title = "asteroid distibution in solution")

plot_asteroid_statistics(asteroids_classical[:, selected_asteroids]; title = "asteroid statistics in solution")

plot_transfer_time_statistics(times_journey; title = "transfer time statistics in solution")




times_intermediate_transfer = [i[3:end-1] - i[2:end-2] for i in times_journey]/day_scale
times_intermediate_transfer_all = vcat(times_intermediate_transfer...)



using Distributions


distribution_intermediate_times = truncated(
    fit(Normal{Float64}, times_intermediate_transfer_all[times_intermediate_transfer_all .<= 500]),
    lower = 50, 
    upper = 300
)

plot_distribution(distribution_intermediate_times)







f = Figure(size = (800, 800), backgroundcolor = :white, figure_padding = 5)

ax = Axis(
    f[1, 1]; 
    limits = (0.5, 3.5, -0.2, 2.8),
    title = "relative positions start vs. final in start RTN frame",
    xlabel = "r [AU]",
    ylabel = "t [AU]"
    # aspect = 1,
)

cartesian_start_states = []
cartesian_end_states = []

rtn_start_states = []
rtn_end_states = []

transfer_times = []
transfer_lambert_costs = []

for (id_journey, times_journey) in zip(id_journey, times_journey)
    for (id_start, id_end, time_start, time_end) in zip(
        id_journey[2:end-2], id_journey[3:end-1],
        times_journey[2:end-2], times_journey[3:end-1],
    )
        if !(50 <= (time_end - time_start)/day_scale <= 300)
            continue
        end

        if id_start ∉ pruned_asteroids || id_end ∉ pruned_asteroids
            continue
        end


        cartesian_start = ephermeris_cartesian_from_id(id_start, time_start)
        cartesian_end = ephermeris_cartesian_from_id(id_end, time_end)

        r_vector = normalize(cartesian_start[1:3])
        n_vector = normalize(cross(r_vector, cartesian_start[4:6]))
        t_vector = normalize(cross(n_vector, r_vector))

        rtn_proj = hcat(r_vector, t_vector, n_vector)'

        rtn_start = rtn_proj * cartesian_start[1:3]
        rtn_end = rtn_proj * cartesian_end[1:3]

        lambert_cost = get_transfer_dv(id_start, time_start, [id_end], time_end - time_start)[1]

        push!(cartesian_start_states, cartesian_start)
        push!(cartesian_end_states, cartesian_end)

        push!(rtn_start_states, vcat(rtn_start, rtn_proj * cartesian_start[4:6]))
        push!(rtn_end_states, vcat(rtn_end, rtn_proj * cartesian_end[4:6]))

        push!(transfer_times, time_end - time_start)
        push!(transfer_lambert_costs, lambert_cost)

        temp_color = (time_end - time_start)/day_scale

        scatter!(
            ax,
            rtn_start[1],
            rtn_start[2],
            rtn_start[3],
            color = temp_color,
            colorrange = (50, 300),
            colormap = :phase,
        )

        scatter!(
            ax,
            rtn_end[1],
            rtn_end[2],
            rtn_end[3],
            color = temp_color,
            colorrange = (50, 300),
            colormap = :phase,
        )
    end
end

Colorbar(f[1, 2], limits = (50, 300), labelsize=16, colormap = :phase, label = "tof [days]")


display(f)




plot_rtn_state_statistics(rtn_start_states, rtn_end_states)




df = DataFrame(
    "start r [AU]" => [k[1] for k in rtn_start_states],
    "start vᵣ [km/s]" => [k[4]*v_scale for k in rtn_start_states],
    "start vₜ [km/s]" => [k[5]*v_scale for k in rtn_start_states],
    "end r [AU]" => [k[1] for k in rtn_end_states],
    "end t [AU]" => [k[2] for k in rtn_end_states],
    "end n [AU]" => [k[3] for k in rtn_end_states],
    "end vᵣ [km/s]" => [k[4]*v_scale for k in rtn_end_states],
    "end vₜ [km/s]" => [k[5]*v_scale for k in rtn_end_states],
    "end vₙ [km/s]" => [k[6]*v_scale for k in rtn_end_states],
)





fig = Figure(size=(1600,1600))


pairplot(
    fig[1, 1],
    df => (
        PairPlots.Scatter(markersize = 5), 
        PairPlots.Contour(),
        PairPlots.MarginHist(),
        PairPlots.MarginConfidenceLimits(),
        PairPlots.MarginDensity(),
        PairPlots.Correlation()
    )
)

display(fig)




scatter(
    transfer_times./day_scale,
    transfer_lambert_costs.*v_scale,
)



plot_tof_dv_comparison(
    transfer_times./day_scale,
    transfer_lambert_costs.*v_scale
)






df = DataFrame(
    "ID" => Int64[],
    "t_n [days]" => Float64[],
    "x_n [AU]" => Float64[],
    "y_n [AU]" => Float64[],
    "z_n [AU]" => Float64[],
    "vx_n [km/s]" => Float64[],
    "vy_n [km/s]" => Float64[],
    "vz_n [km/s]" => Float64[],
    "ux_n [N]" => Float64[],
    "uy_n [N]" => Float64[],
    "uz_n [N]" => Float64[],
    "m_n [kg]" => Float64[],
)

maximum_initial_masses = Float64[]
maximum_final_masses = Float64[]

id = 1

for (id_journey, times_journey) in zip(id_journey, times_journey)
    for (id_start, id_end, time_start, time_end) in zip(
        id_journey[2:end-2], id_journey[3:end-1],
        times_journey[2:end-2], times_journey[3:end-1],
    )
        if !(50 <= (time_end - time_start)/day_scale <= 300)
            continue
        end

        if id_start ∉ pruned_asteroids || id_end ∉ pruned_asteroids
            continue
        end

        node_time_spacing = (time_end - time_start) / 99

        prob = SequentialConvexProblem(
            [[id_start, id_end]],
            [[time_start, time_end]],
        )

        solve!(prob;
            fixed_segments = true,
            fixed_rendezvous = true,
            maximum_mass = true
        )

        convert_logged_mass_to_mass!(prob)

        solve!(prob;
            fixed_segments = true,
            fixed_rendezvous = true,
            maximum_mass = true
        )

        solve!(prob;
            fixed_segments = true,
            fixed_rendezvous = true,
            maximum_mass = true
        )

        for i in 1:length(prob.t_nodes[1][1])
            control = if i == length(prob.t_nodes[1][1])
                [0.0, 0.0, 0.0]
            else
                prob.u_nodes[1][1][1:3, i]
            end

            push!(df, (
                id,
                prob.t_nodes[1][1][i] / day_scale,
                prob.x_nodes[1][1][1, i],
                prob.x_nodes[1][1][2, i],
                prob.x_nodes[1][1][3, i],
                prob.x_nodes[1][1][4, i] * v_scale,
                prob.x_nodes[1][1][5, i] * v_scale,
                prob.x_nodes[1][1][6, i] * v_scale,
                control[1] * thrust * m_scale * a_scale * 1e3,
                control[2] * thrust * m_scale * a_scale * 1e3,
                control[3] * thrust * m_scale * a_scale * 1e3,
                prob.x_nodes[1][1][7, i] * m_scale,
            ))

           
        end

        push!(maximum_initial_masses, prob.x_nodes[1][1][7, 1]*m_scale)
        push!(maximum_final_masses, prob.x_nodes[1][1][7, end]*m_scale)


        id += 1
    end
end


CSV.write("output/transfer_information.csv", df)




plot_mass_use_statistics(
    maximum_initial_masses,
    maximum_final_masses,
)



scatter(
    transfer_times./day_scale,
    maximum_initial_masses
)


scatter(
    transfer_lambert_costs.*v_scale,
    maximum_initial_masses
)



rocket_equation_dv = g0_isp .* log.(maximum_initial_masses ./ maximum_final_masses)



df = DataFrame(
    "tof [days]" => transfer_times./day_scale,
    "lambert Δv [km/s]" => transfer_lambert_costs.*v_scale,
    "low-thrust maximum Δv [km/s]" => rocket_equation_dv.*v_scale,
    "maximum initial mass [kg]" => maximum_initial_masses,
    "maximum final mass [kg]" => maximum_final_masses,
)





fig = Figure(size=(1000,1000))


pairplot(
    fig[1, 1],
    df => (
        PairPlots.Scatter(markersize = 5), 
        PairPlots.Contour(),
        PairPlots.MarginHist(),
        PairPlots.MarginConfidenceLimits(),
        PairPlots.MarginDensity(),
        PairPlots.Correlation()
    )
)

display(fig)











df = DataFrame(
    "x0 [AU]" => [k[1] for k in cartesian_start_states],
    "y0 [AU]" => [k[2] for k in cartesian_start_states],
    "z0 [AU]" => [k[3] for k in cartesian_start_states],
    "vx0 [km/s]" => [k[4] for k in cartesian_start_states],
    "vy0 [km/s]" => [k[5] for k in cartesian_start_states],
    "vz0 [km/s]" => [k[6] for k in cartesian_start_states],
    "x1 [AU]" => [k[1] for k in cartesian_end_states],
    "y1 [AU]" => [k[2] for k in cartesian_end_states],
    "z1 [AU]" => [k[3] for k in cartesian_end_states],
    "vx1 [km/s]" => [k[4] for k in cartesian_end_states],
    "vy1 [km/s]" => [k[5] for k in cartesian_end_states],
    "vz1 [km/s]" => [k[6] for k in cartesian_end_states],
    "r0 [AU]" => [k[1] for k in rtn_start_states],
    "t0 [AU]" => [k[2] for k in rtn_start_states],
    "n0 [AU]" => [k[3] for k in rtn_start_states],
    "vr0 [km/s]" => [k[4] for k in rtn_start_states],
    "vt0 [km/s]" => [k[5] for k in rtn_start_states],
    "vn0 [km/s]" => [k[6] for k in rtn_start_states],
    "r1 [AU]" => [k[1] for k in rtn_end_states],
    "t1 [AU]" => [k[2] for k in rtn_end_states],
    "n1 [AU]" => [k[3] for k in rtn_end_states],
    "vr1 [km/s]" => [k[4] for k in rtn_end_states],
    "vt1 [km/s]" => [k[5] for k in rtn_end_states],
    "vn1 [km/s]" => [k[6] for k in rtn_end_states],
    "tof [days]" => transfer_times./day_scale,
    "Δv_lambert [km/s]" => transfer_lambert_costs.*v_scale,
    "Δv_lowthrust_maximum [km/s]" => rocket_equation_dv.*v_scale,
    "m0_maximum [kg]" => maximum_initial_masses,
    "m1_maximum [kg]" => maximum_final_masses,
)


CSV.write("output/transfer_statistics.csv", df)









df = DataFrame(CSV.File("output/reachability.csv"))

starting_locations = Matrix(df[40000:end, 1:6])'
ending_locations = Matrix(df[40000:end, 7:12])'

starting_locations[4:6, :] ./= v_scale
ending_locations[4:6, :] ./= v_scale

maximum_initial_masses = Float64[]
maximum_final_masses = Float64[]


for (starting_location, ending_location) in zip(eachcol(starting_locations), eachcol(ending_locations))
    custom_cartesian[:, 1] = starting_location
    custom_cartesian[:, 2] = ending_location

    prob = SequentialConvexProblem(
        [[60001, 60002]],
        [[0.0, tof]],
    )

    solve!(prob;
        fixed_segments = true,
        fixed_rendezvous = true,
        maximum_mass = true
    )

    convert_logged_mass_to_mass!(prob)

    solve!(prob;
        fixed_segments = true,
        fixed_rendezvous = true,
        maximum_mass = true
    )

    solve!(prob;
        fixed_segments = true,
        fixed_rendezvous = true,
        maximum_mass = true
    )

    push!(maximum_initial_masses, prob.x_nodes[1][1][7, 1]*m_scale)
    push!(maximum_final_masses, prob.x_nodes[1][1][7, end]*m_scale)

    break
end







































# starting_locations = Matrix(df[:, 1:6])'
# ending_locations = Matrix(df[:, 7:12])'

# starting_locations[4:6, :] ./= v_scale
# ending_locations[4:6, :] ./= v_scale

node_time_spacing = 10.0*day_scale
scp_iterations = 20

# Circular Orbit
starting_location = [
    3.0,
    0.0,
    0.0,
    0.0,
    sqrt(1.0/3.0),
    0.0
]

tof = 500.0*day_scale

ending_location = integrate_arc(starting_location, [0.0, 0.0, 0.0, 0.0], [0.0, tof])



location_delta_x = LinRange(-1.0, 1.0, 50)
location_delta_y = LinRange(-1.0, 1.0, 50)


maximum_initial_transfer_masses = zeros(50, 50)
maximum_final_transfer_masses = zeros(50, 50)

custom_cartesian = zeros(6, 2)


for (i, delta_x) in enumerate(location_delta_x), (j, delta_y) in enumerate(location_delta_y)
    custom_cartesian[:, 1] = starting_location
    custom_cartesian[:, 2] = ending_location

    adjusted_ending_location = copy(ending_location)
    adjusted_ending_location[1] += delta_x
    adjusted_ending_location[2] += delta_y

    # Calculate the new velocity components for a circular orbit
    r = sqrt(adjusted_ending_location[1]^2 + adjusted_ending_location[2]^2)
    v = sqrt(1.0 / r)
    adjusted_ending_location[4] = -v * (adjusted_ending_location[2] / r)
    adjusted_ending_location[5] = v * (adjusted_ending_location[1] / r)
    adjusted_ending_location[6] = 0.0

    custom_cartesian[:, 2] = adjusted_ending_location

    prob = SequentialConvexProblem(
        [[60001, 60002]],
        [[0.0, tof]],
    )

    solve!(prob;
        fixed_segments = true,
        fixed_rendezvous = true,
        maximum_mass = true
    )

    convert_logged_mass_to_mass!(prob)

    solve!(prob;
        fixed_segments = true,
        fixed_rendezvous = true,
        maximum_mass = true
    )

    maximum_initial_transfer_masses[i, j] = prob.x_nodes[1][1][7, 1]*m_scale
    maximum_final_transfer_masses[i, j] = prob.x_nodes[1][1][7, end]*m_scale

    break
end


# Plot the contour of the mass use
f = Figure(size = (800, 800), backgroundcolor = :white, figure_padding = 0)

ax = Axis(
    f[1, 1]; 
    xlabel = "Δx [AU]", 
    ylabel = "Δy [AU]",
    title = "maximum initial mass for 3.0 AU circular orbit to circular orbit transfers, tof = 800 days"
)

cont = contourf!(
    ax,
    location_delta_x,
    location_delta_y,
    maximum_final_transfer_masses,
    levels = 20,
    colormap = :viridis
)

# Add Colorbar
Colorbar(f[1, 2], cont, labelsize=16, label = "maximum final mass [kg]", tickformat = "{:.0f}")

display(f)





node_time_spacing = 10.0*day_scale
scp_iterations = 50



# Plotting provided data
df = DataFrame(CSV.File("output/reachability_01_nonrotated.csv"))

# Convert row 13 to Float64
df[!, 13] = Float64.(df[!, 13])

# Add row for maximum initial and final mass to df
df[!, "maximum_initial_mass"] .= 0.0
df[!, "maximum_final_mass"] .= 0.0

custom_cartesian = zeros(6, 2)


temp = []

for row in eachrow(df)
    custom_cartesian[:, 1] = Vector(row[1:6])
    custom_cartesian[:, 2] = Vector(row[7:12])

    custom_cartesian[4:6, 1] ./= v_scale
    custom_cartesian[4:6, 2] ./= v_scale

    push!(temp, get_transfer_dv(60001, 0.0, [60002], 50.0*day_scale)[1]*v_scale)
end

argmin(temp)



get_transfer_dv(60001, 0.0, [60002], 50.0*day_scale)*v_scale


plot_trajectory(prob)


for row in eachrow(df)
    custom_cartesian[:, 1] = Vector(row[1:6])
    custom_cartesian[:, 2] = Vector(row[7:12])

    custom_cartesian[4:6, 1] ./= v_scale
    custom_cartesian[4:6, 2] ./= v_scale

    prob = SequentialConvexProblem(
        [[60001, 60002]],
        [[0.0, row[13]*day_scale]],
    )

    solve!(prob;
        fixed_segments = true,
        fixed_rendezvous = true,
        maximum_mass = true
    )

    convert_logged_mass_to_mass!(prob)

    solve!(prob;
        fixed_segments = true,
        fixed_rendezvous = true,
        maximum_mass = true
    )

    row[14] = prob.x_nodes[1][1][7, 1]*m_scale
    row[15] = prob.x_nodes[1][1][7, end]*m_scale

    # break
end

CSV.write("output/reachability_01_nonrotated_mass.csv", df)

# Plot the contour of the mass use
f = Figure(size = (800, 800), backgroundcolor = :white, figure_padding = 0)

ax = Axis(
    f[1, 1]; 
    xlabel = "x1 [AU]", 
    ylabel = "y1 [AU]",
    title = "maximum initial mass for 3.0 AU circular orbit to circular orbit transfers, tof = 50 days"
)

cont = tricontourf!(
    ax,
    Vector(df[!, 7]),
    Vector(df[!, 8]),
    Vector(df[!, 14]),
    levels = 20,
    colormap = :viridis
)

# Add Colorbar
Colorbar(f[1, 2], cont, labelsize=16, label = "maximum final mass [kg]", tickformat = "{:.0f}")

display(f)

