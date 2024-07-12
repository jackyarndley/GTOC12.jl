


mutable struct MixedIntegerProblem{T <: Real}
    id_journey::Union{Nothing, Vector{Vector{T}}}
    times_journey::Vector{Vector{T}}
    id_journey_solutions::Union{Nothing, Vector{Vector{Vector{Int64}}}}
    id_subset::Vector{Int64}
    deployment_nums::Vector{Int64}
    deployment_arcs::Vector{Int64}
    collection_nums::Vector{Int64}
    collection_arcs::Vector{Int64}
    solutions::Int64
    mixing_number::Int64
    deployment_cost::Vector{Array{T, 3}}
    collection_cost::Vector{Array{T, 3}}
    intermediate_cost::Vector{Array{T, 2}}
    cost_limit::Float64
end


function MixedIntegerProblem(
    id_subset,
    deployment_nums,
    collection_nums;
    times_journey = nothing,
    cost_limit = 6.0/v_scale,
)
    if length(id_subset) > 250
        error("Subset size is too big")
    end

    if length(deployment_nums) != length(collection_nums)
        error("Size of deployments should be same as collections")
    end

    mixing_number = length(deployment_nums)

    if isnothing(times_journey)
        times_journey = generate_default_times(
            11, 
            2, 
            id_subset;
            algorithm = 3,
            # time_parameter_days = 135
            time_parameter_days = 145
            # time_parameter_days = 175
        )

        times_journey = [
            vcat(
                times_journey[n][1],
                times_journey[n][2:(deployment_nums[n] + 1)],
                times_journey[n][(end - collection_nums[n]):(end-1)],
                times_journey[n][end],
            ) for n in 1:mixing_number
        ]
    end


    times_deployment = [times_journey[n][2:(deployment_nums[n] + 1)] for n in 1:mixing_number]
    times_collection = [times_journey[n][end-collection_nums[n]:end-1] for n in 1:mixing_number]

    deployment_arcs = deployment_nums .- 1
    collection_arcs = collection_nums .- 1

    deployment_cost = [
        stack(stack(get_transfer_dv(id, times_deployment[n][i], id_subset, times_deployment[n][i + 1] - times_deployment[n][i]) for id in id_subset) for i in 1:deployment_arcs[n])
        for n in 1:mixing_number
    ]

    collection_cost = [
        stack(stack(get_transfer_dv(id, times_collection[n][i], id_subset, times_collection[n][i + 1] - times_collection[n][i]) for id in id_subset) for i in 1:collection_arcs[n])
        for n in 1:mixing_number
    ]

    intermediate_cost = 1.0.*[
        stack(stack(get_transfer_dv(id, times_deployment[n][end] + max(0.0, times_collection[n][1] - times_deployment[n][end] - 500*day_scale), id_subset, times_collection[n][1] - times_deployment[n][end] - max(0.0, times_collection[n][1] - times_deployment[n][end] - 500*day_scale)) for id in id_subset))
        for n in 1:mixing_number
    ]

    return MixedIntegerProblem(
        nothing,
        times_journey,
        nothing,
        id_subset,
        deployment_nums,
        deployment_arcs,
        collection_nums,
        collection_arcs,
        0,
        mixing_number,
        deployment_cost,
        collection_cost,
        intermediate_cost,
        cost_limit
    )
end


function solve!(
    p::MixedIntegerProblem;
    time_limit_seconds = 60.0,
    fix_start_asteroids = nothing,
    self_cleaning = false,
    permit_intermediate_transfer = true,
    include_intermediate_transfer_cost = false,
    solutions_relative_allowance = 0.1,
    solutions_count_maximum = 100,
)
    subset_size = length(p.id_subset)

    model = Model(Gurobi.Optimizer)
    set_attribute(model, "TimeLimit", time_limit_seconds)
    set_attribute(model, "MIPFocus", 1)
    set_attribute(model, "Heuristics", 0.1)
    set_attribute(model, "ConcurrentMIP", 3)
    set_attribute(model, "Method", 1)
    set_attribute(model, "ImproveStartTime", 0.8*time_limit_seconds)

    # For larger subset size
    if subset_size > 60
        set_attribute(model, "NoRelHeurTime", 0.1*time_limit_seconds)
    end

    deployment = []
    collection = []
    intermediate = []

    for n in 1:p.mixing_number
        if isnothing(p.id_journey)
            # Create variables for each permitted deployment arc
            push!(deployment, @variable(model, [i=1:subset_size, j=1:subset_size, k=1:p.deployment_arcs[n]; (i != j) && ((p.deployment_cost[n][i, j, k] <= p.cost_limit))], Bin))

            # Create variables for each permitted collection arc
            push!(collection, @variable(model, [i=1:subset_size, j=1:subset_size, k=1:p.collection_arcs[n]; (i != j) && ((p.collection_cost[n][i, j, k] <= p.cost_limit))], Bin))

            if include_intermediate_transfer_cost
                push!(intermediate, @variable(model, [i=1:subset_size, j=1:subset_size; (p.intermediate_cost[n][i, j] <= p.cost_limit)], Bin))
            end
        else
            # Get the indices of each of the deployment asteroid ID in id_subset
            deployment_ids = [findfirst(==(p.id_journey[n][2:end-1][k]), id_subset) for k in 1:p.deployment_nums[n]]
        
            # Get the indices of each of the collection asteroid ID in id_subset
            collection_ids = [findfirst(==(p.id_journey[n][end-p.collection_nums[n]:end-1][k]), id_subset) for k in 1:p.collection_nums[n]]

            # Create variables for each permitted deployment arc, also allowing those in the initial solution
            push!(deployment, @variable(model, [i=1:subset_size, j=1:subset_size, k=1:p.deployment_arcs[n]; (i != j) && ((p.deployment_cost[n][i, j, k] <= p.cost_limit) || ((deployment_ids[k]) == i && (deployment_ids[k+1] == j)))], Bin))

            # Create variables for each permitted collection arc, also allowing those in the initial solution
            push!(collection, @variable(model, [i=1:subset_size, j=1:subset_size, k=1:p.collection_arcs[n]; (i != j) && ((p.collection_cost[n][i, j, k] <= p.cost_limit) || ((collection_ids[k]) == i && (collection_ids[k+1] == j)))], Bin))

            for k in 1:p.deployment_arcs[n]
                set_start_value(deployment[n][deployment_ids[k], deployment_ids[k+1], k], 1)
            end

            for k in 1:p.collection_arcs[n]
                set_start_value(collection[n][collection_ids[k], collection_ids[k+1], k], 1)
            end

            if include_intermediate_transfer_cost
                push!(intermediate, @variable(model, [i=1:subset_size, j=1:subset_size; (p.intermediate_cost[n][i, j] <= p.cost_limit) || ((deployment_ids[n][end] == i) && (collection_ids[n][1] == j))], Bin))

                set_start_value(intermediate[n][deployment_ids[end], collection_ids[1]], 1)
            end
        end

        # Only one arc selected per time segment
        @constraint(model, [k=1:p.deployment_arcs[n]], sum(deployment[n][:, :, k]) == 1)
        @constraint(model, [k=1:p.collection_arcs[n]], sum(collection[n][:, :, k]) == 1)

        # If an arc goes into a node it must come out
        @constraint(model, [i=1:subset_size, k=1:(p.deployment_arcs[n]-1)], sum(deployment[n][:, i, k]) == sum(deployment[n][i, :, k + 1]))
        @constraint(model, [i=1:subset_size, k=1:(p.collection_arcs[n]-1)], sum(collection[n][:, i, k]) == sum(collection[n][i, :, k + 1]))
    
        if !isnothing(fix_start_asteroids)
            @constraint(model, sum(deployment[findfirst(==(fix_start_asteroids[n]), subset_ids)][i, :, 1]) == 1)
        end

        if include_intermediate_transfer_cost
            @constraint(model, sum(intermediate[n]) == 1)

            [@constraint(model, 2*intermediate[n][i, j] .<= sum(deployment[n][:, i, p.deployment_arcs[n]]) + sum(collection[n][j, :, 1])) for (i, j) in eachindex(intermediate[n])]
        end

        # Only permit low cost first transfers
        # @constraint(model, [i=1:subset_size, j=1:subset_size; (i != j) && (p.deployment_cost[n][i, j, 1] <= p.cost_limit) && p.deployment_cost[n][i, j, 1] >= 3.8/v_scale], deployment[n][i, j, 1] == 0)

        # Force staying at the asteroid between the deploy and collect
        if !permit_intermediate_transfer
            @constraint(model, [i=1:subset_size], sum(deployment[n][:, i, p.deployment_arcs[n]]) == sum(collection[n][i, :, 1]))
        end
    end

    # Prevent more than 1 deploy at an asteroid
    @constraint(model, [i=1:subset_size], 
        sum([sum(deployment[n][i, :, :]) for n in 1:p.mixing_number]) + sum([sum(deployment[n][:, i, p.deployment_arcs[n]]) for n in 1:p.mixing_number]) <= 1
    )

    # Prevent more than 1 collect at an asteroid
    @constraint(model, [i=1:subset_size], 
        sum([sum(collection[n][i, :, :]) for n in 1:p.mixing_number]) + sum([sum(collection[n][:, i, p.collection_arcs[n]]) for n in 1:p.mixing_number]) <= 1
    )

    # Must collect if deployed at
    if self_cleaning
        [@constraint(model, [i=1:subset_size], 
            sum(deployment[n][i, :, :]) + sum(deployment[n][:, i, p.deployment_arcs[n]]) == sum(collection[n][i, :, :]) + sum(collection[n][:, i, p.collection_arcs[n]])
        ) for n in 1:p.mixing_number]
    else
        @constraint(model, [i=1:subset_size], 
            sum([sum(deployment[n][i, :, :]) for n in 1:p.mixing_number]) + sum([sum(deployment[n][:, i, p.deployment_arcs[n]]) for n in 1:p.mixing_number]) == sum([sum(collection[n][i, :, :]) for n in 1:p.mixing_number]) + sum([sum(collection[n][:, i, p.collection_arcs[n]]) for n in 1:p.mixing_number])
        )
    end

    if include_intermediate_transfer_cost
        @objective(model, Min, 
            sum([p.deployment_cost[n][i, j, k] * deployment[n][i, j, k] for n in 1:p.mixing_number for (i, j, k) in eachindex(deployment[n])]) + sum([p.collection_cost[n][i, j, k] * collection[n][i, j, k] for n in 1:p.mixing_number for (i, j, k) in eachindex(collection[n])]) + sum([p.intermediate_cost[n][i, j] * intermediate[n][i, j] for n in 1:p.mixing_number for (i, j) in eachindex(intermediate[n])])
        )
    else
        @objective(model, Min, 
            sum([p.deployment_cost[n][i, j, k] * deployment[n][i, j, k] for n in 1:p.mixing_number for (i, j, k) in eachindex(deployment[n])]) + sum([p.collection_cost[n][i, j, k] * collection[n][i, j, k] for n in 1:p.mixing_number for (i, j, k) in eachindex(collection[n])])
        )
    end

    if solutions_count_maximum > 1
        set_attribute(model, "PoolSearchMode", 2)
        set_attribute(model, "PoolGap", solutions_relative_allowance)
        set_attribute(model, "PoolSolutions", solutions_count_maximum)
    end

    optimize!(model)

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        println("\nNo solution found.")
        return
    end

    p.id_journey_solutions = []

    # Add the default solution if it is provided
    if !isnothing(p.id_journey)
        push!(p.id_journey_solutions, p.id_journey)
        p.solutions += 1
    end

    for i in 1:min(result_count(model), solutions_count_maximum)
        deployment_value = [value.(x; result = i) for x in deployment]
        collection_value = [value.(x; result = i) for x in collection]

        deployment_ids = [sort(filter(vals -> x[vals[1], vals[2], vals[3]] >= 0.5, collect(eachindex(x))); by=x->x[3]) for x in deployment_value]
        collection_ids = [sort(filter(vals -> x[vals[1], vals[2], vals[3]] >= 0.5, collect(eachindex(x))); by=x->x[3]) for x in collection_value]

        id_journey_solution = [
            vcat(
                [0], 
                p.id_subset[vcat([x[1] for x in deployment_ids], [deployment_ids[end][2]])],
                p.id_subset[vcat([x[1] for x in collection_ids], [collection_ids[end][2]])],
                [-3]
            ) for (deployment_ids, collection_ids) in zip(deployment_ids, collection_ids)
        ]

        if id_journey_solution ∉ p.id_journey_solutions
            push!(p.id_journey_solutions, id_journey_solution)

            p.solutions += 1
        end
    end
end


function add_lowest_cost_intermediate_transfer(
    times_journey_mixed,
    id_journey_mixed,
    ids_to_choose,
    deploy_add_index = 1,
    collect_add_index = 1;
    transfer_Δt = 10.0
)
    deploy_lengths, collect_lengths = get_deploy_collect_lengths(
        id_journey_mixed, 
        times_journey_mixed
    )

    dv_deploy = get_transfer_dv(
        id_journey_mixed[deploy_add_index][deploy_lengths[deploy_add_index] + 1], 
        times_journey_mixed[deploy_add_index][deploy_lengths[deploy_add_index] + 1],
        ids_to_choose,
        transfer_Δt
    )
    
    dv_collect = [get_transfer_dv(
        id, 
        times_journey_mixed[collect_add_index][end - collect_lengths[collect_add_index]] - transfer_Δt,
        id_journey_mixed[collect_add_index][end - collect_lengths[collect_add_index]],
        transfer_Δt
    )[1] for id in ids_to_choose]
    
    println("Mininum $(v_scale*minimum(dv_deploy + dv_collect))")

    new_id = ids_to_choose[argmin(dv_deploy + dv_collect)]


    times_journey_mixed[deploy_add_index] = vcat(
        times_journey_mixed[deploy_add_index][1:deploy_lengths[deploy_add_index] + 1], 
        times_journey_mixed[deploy_add_index][deploy_lengths[deploy_add_index] + 1] + transfer_Δt, 
        times_journey_mixed[deploy_add_index][deploy_lengths[deploy_add_index] + 2:end]
    )

    times_journey_mixed[collect_add_index] = vcat(
        times_journey_mixed[collect_add_index][1:end-collect_lengths[collect_add_index] - 1], 
        times_journey_mixed[collect_add_index][end-collect_lengths[collect_add_index]] - transfer_Δt, 
        times_journey_mixed[collect_add_index][end-collect_lengths[collect_add_index]:end]
    )

    id_journey_mixed[deploy_add_index] = vcat(
        id_journey_mixed[deploy_add_index][1:deploy_lengths[deploy_add_index] + 1], 
        new_id, 
        id_journey_mixed[deploy_add_index][deploy_lengths[deploy_add_index] + 2:end]
    )

    id_journey_mixed[collect_add_index] = vcat(
        id_journey_mixed[collect_add_index][1:end-collect_lengths[collect_add_index] - 1], 
        new_id, 
        id_journey_mixed[collect_add_index][end-collect_lengths[collect_add_index]:end]
    )


    return times_journey_mixed, id_journey_mixed
end