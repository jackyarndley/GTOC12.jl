include("header.jl")

GLMakie.activate!()
CairoMakie.activate!()

using PairPlots
using DataFrames

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





maximum_initial_masses = Float64[]
maximum_final_masses = Float64[]


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

        prob = SequentialConvexProblem(
            [[id_start, id_end]],
            [[time_start, time_end]],
        )

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

        push!(maximum_initial_masses, exp(prob.x_nodes[1][1][7, 1])*m_scale)
        push!(maximum_final_masses, exp(prob.x_nodes[1][1][7, end])*m_scale)

    end
end



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















