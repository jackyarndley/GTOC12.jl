function mip_ordering_mixing(
    journey_length, 
    subset_ids, 
    times_journey_mixed;
    time_limit_seconds = 180,
    transfer_dv_limit = 3.5,
    allow_middle_transfer = true,
    include_middle_transfer_dv = true,
    self_cleaning = false,
)
    # Subtract 1 because there are n-1 transfers for n asteroid tours
    journey_length -= 1

    times_deploy = [times_journey[2:(journey_length + 2)] for times_journey in times_journey_mixed]
    times_collect = [times_journey[(journey_length + 3):(end-1)] for times_journey in times_journey_mixed]

    # model = Model(HiGHS.Optimizer)
    model = Model(Gurobi.Optimizer)
    set_attribute(model, "TimeLimit", time_limit_seconds)
    set_attribute(model, "MIPFocus", 1)

    if length(subset_ids) >= 200
        error("Subset size is quite big!")
    end

    # Binary variables for arcs
    deployer = [@variable(model, [1:journey_length, 1:length(subset_ids), 1:length(subset_ids)], Bin) for _ in times_journey_mixed]
    collector = [@variable(model, [1:journey_length, 1:length(subset_ids), 1:length(subset_ids)], Bin) for _ in times_journey_mixed]

    # Only one arc selected per time
    [@constraint(model, [i=1:journey_length], sum(deployer[i, :, :]) == 1) for deployer in deployer]
    [@constraint(model, [i=1:journey_length], sum(collector[i, :, :]) == 1) for collector in collector]

    # Cannot have both to and from for start and end arcs
    @constraint(model, [i=1:length(subset_ids)], sum([sum(deployer[end, :, i]) + sum(deployer[1, i, :]) for deployer in deployer]) <= 1)
    @constraint(model, [i=1:length(subset_ids)], sum([sum(collector[end, :, i]) + sum(collector[1, i, :]) for collector in collector]) <= 1)

    # Do not allow self arcs
    @constraint(model, [i=1:length(subset_ids)], sum([sum(deployer[:, i, i]) for deployer in deployer]) == 0)
    @constraint(model, [i=1:length(subset_ids)], sum([sum(collector[:, i, i]) for collector in collector]) == 0)

    # Max one from/to arc selected per id
    @constraint(model, [i=1:length(subset_ids)], sum([sum(deployer[:, :, i]) for deployer in deployer]) <= 1)
    @constraint(model, [i=1:length(subset_ids)], sum([sum(deployer[:, i, :]) for deployer in deployer]) <= 1)
    @constraint(model, [i=1:length(subset_ids)], sum([sum(collector[:, :, i]) for collector in collector]) <= 1)
    @constraint(model, [i=1:length(subset_ids)], sum([sum(collector[:, i, :]) for collector in collector]) <= 1)

    # Continuity constraint between arcs
    [@constraint(model, [i=1:(journey_length-1), j=1:length(subset_ids)], sum(deployer[i, :, j]) == sum(deployer[i + 1, j, :])) for deployer in deployer]
    [@constraint(model, [i=1:(journey_length-1), j=1:length(subset_ids)], sum(collector[i, :, j]) == sum(collector[i + 1, j, :])) for collector in collector]

    # Force staying at the asteroid between the deploy and collect
    if !allow_middle_transfer
        [@constraint(model, [i=1:length(subset_ids)], sum(deployer[n][end, :, i]) == sum(collector[n][1, i, :])) for n in 1:length(deployer)]
    end

    # Force arcs to match between deploy and collect
    if self_cleaning
        [@constraint(model, [i=1:length(subset_ids)], 0.1* (sum(deployer[:, i, :]) + sum(deployer[:, :, i])) <= sum(collector[:, :, i]) + sum(collector[:, i, :])) for (deployer, collector) in zip(deployer, collector)]
    else
        @constraint(model, [i=1:length(subset_ids)], 0.1*sum([sum(deployer[:, i, :]) + sum(deployer[:, :, i]) for deployer in deployer]) <= sum([sum(collector[:, i, :]) + sum(collector[:, :, i]) for collector in collector]))
    end

    # Set the correct Δv
    Δv_deploy = [@expression(model, [i=1:journey_length], sum([sum(deployer[n][i, j, :] .* get_transfer_dv(subset_ids[j], times_deploy[n][i], subset_ids, times_deploy[n][i + 1] - times_deploy[n][i])) for j in 1:length(subset_ids)])) for n in 1:length(deployer)]

    Δv_collect = [@expression(model, [i=1:journey_length], sum([sum(collector[n][i, j, :] .* get_transfer_dv(subset_ids[j], times_collect[n][i], subset_ids, times_collect[n][i + 1] - times_collect[n][i])) for j in 1:length(subset_ids)])) for n in 1:length(collector)]

    cutoff = transfer_dv_limit/v_scale

    [@constraint(model, [i=1:journey_length, j=1:length(subset_ids), k=1:length(subset_ids); get_transfer_dv(subset_ids[j], times_deploy[n][i], subset_ids[k], times_deploy[n][i + 1] - times_deploy[n][i])[1] >= cutoff], deployer[n][i, j, k] == 0) for n in 1:length(deployer)]

    [@constraint(model, [i=1:journey_length, j=1:length(subset_ids), k=1:length(subset_ids); get_transfer_dv(subset_ids[j], times_collect[n][i], subset_ids[k], times_collect[n][i + 1] - times_collect[n][i])[1] >= cutoff], collector[n][i, j, k] == 0) for n in 1:length(collector)]

    if include_middle_transfer_dv
        midpoint = [@variable(model, [i=1:length(subset_ids), j=1:length(subset_ids)], Bin) for _ in times_journey_mixed]

        [@constraint(model, sum(midpoint) == 1) for midpoint in midpoint]
        [@constraint(model, [i=1:length(subset_ids), j=1:length(subset_ids)], midpoint[n][i, j] .>= -1.0 + sum(deployer[n][end, :, i]) + sum(collector[n][1, j, :])) for n in 1:length(midpoint)]

        Δv_midpoint = [@expression(model, 0.5*sum(sum([midpoint[n][i, :] .* get_transfer_dv(subset_ids[i], times_deploy[n][end] + max(0.0, times_collect[n][1] - times_deploy[n][end] - 500*day_scale), subset_ids, times_collect[n][1] - times_deploy[n][end] - max(0.0, times_collect[n][1] - times_deploy[n][end] - 500*day_scale)) for i in 1:length(subset_ids)]))) for n in 1:length(midpoint)]

        @objective(model, Min, sum([sum(Δv_deploy[i]) + sum(Δv_collect[i]) + Δv_midpoint[i] for i in 1:length(Δv_deploy)]))
    else
        @objective(model, Min, sum([sum(Δv_deploy[i]) + sum(Δv_collect[i]) for i in 1:length(Δv_deploy)]))
    end

    # Δv_initial = [@expression(model, sum(sum.(eachrow(deployer[n][1, :, :])) .* get_transfer_dv_low_thrust(0, times_journey_mixed[n][1], subset_ids, times_deploy[n][1] - times_journey_mixed[n][1]))) for n in 1:length(deployer)]
        
    # Δv_final = [@expression(model, sum(sum.(eachcol(collector[n][end, :, :])) .* [get_transfer_dv_low_thrust(subset_id, times_collect[n][end], [-3], times_journey_mixed[n][end] - times_collect[n][end]; start_mass=1400/m_scale)[1] for subset_id in subset_ids]
    # )) for n in 1:length(deployer)]
    
    # @objective(model, Min, sum([sum(Δv_deploy[i]) + sum(Δv_collect[i]) + Δv_initial[i] + Δv_final[i] for i in 1:length(Δv_deploy)]))

    optimize!(model)

    deployer_ids = [sort(findall(>=(0.9999), value.(deployer)); by=x->x[1]) for deployer in deployer]
    collector_ids = [sort(findall(>=(0.9999), value.(collector)); by=x->x[1]) for collector in collector]

    id_journey_mixed = [
        vcat(
            [0], 
            subset_ids[vcat([x[2] for x in deployer_ids], [deployer_ids[end][3]])],
            subset_ids[vcat([x[2] for x in collector_ids], [collector_ids[end][3]])],
            [-3]
        ) for (deployer_ids, collector_ids) in zip(deployer_ids, collector_ids)
    ]

    return id_journey_mixed
end