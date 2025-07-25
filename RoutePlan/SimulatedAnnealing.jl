using Distributions
using Random

function simulated_annealing_timed(init_state, time_limit, prob_function, cost_function, temp_function, neighbors_function)
    current_state = init_state
    current_cost = cost_function(current_state)
    best_state = current_state
    best_cost = current_cost
    visited_list = []

    start = time()
    while time() - start < time_limit
        neighbors = neighbors_function(current_state)
        neighbor = rand(neighbors)
        neighbor_cost = cost_function(neighbor)
        temperature = temp_function(1 - (time() - start) / time_limit)

        if prob_function(current_cost, neighbor_cost, temperature) > rand()
            current_state = neighbor
            current_cost = neighbor_cost
            if current_cost < best_cost
                best_state = current_state
                best_cost = current_cost
            end
            push!(visited_list, [current_state, time() - start, current_cost])
        end
    end
    return best_state, best_cost, visited_list
end

function simulated_annealing_iter(init_state, iter_limit, prob_function, cost_function, temp_function, neighbors_function)
    current_state = init_state
    current_cost = cost_function(current_state)
    best_state = current_state
    best_cost = current_cost
    visited_list = []

    iters = 0
    while iters < iter_limit
        neighbors = neighbors_function(current_state)
        neighbor = rand(neighbors)
        neighbor_cost = cost_function(neighbor)
        temperature = temp_function(1 - (iters / iter_limit))

        if prob_function(current_cost, neighbor_cost, temperature) > rand()
            current_state = neighbor
            current_cost = neighbor_cost
            if current_cost < best_cost
                best_state = current_state
                best_cost = current_cost
            end
            push!(visited_list, [current_state, iters, current_cost])
        end
        iters += 1
    end
    return best_state, best_cost, visited_list
end

function run_simulated_annealing(limit, cost_dist, divisions, neighborhood_size, type="timed")
    function temp_function(r)
        return exp(-r)
    end

    function prob_function(current_cost, neighbor_cost, temperature)
        if neighbor_cost < current_cost
            return 1.0
        else
            return exp((current_cost - neighbor_cost) / temperature)
        end
    end

    function neighbors_function(state)
        return range(Base.max(state - neighborhood_size, 1), Base.min(state + neighborhood_size, divisions))
    end

    function cost_function(state)
        return cost_dist[state]
    end

    function init_state()
        return rand(range(1, divisions))
    end

    init = init_state()
    if type == "timed"
        best_state, best_cost, visited_list = simulated_annealing_timed(init, limit, prob_function, cost_function, temp_function, neighbors_function)
    else
        start = time()
        best_state, best_cost, visited_list = simulated_annealing_iter(init, limit, prob_function, cost_function, temp_function, neighbors_function)
        println("Took $(time() - start) seconds")
    end
    
    return best_state, best_cost, visited_list
end

min_cost = 0.0
max_cost = 1.0
divisions = 100000000
time_limit = 0.5  # seconds
iter_limit = 10000000000  # iterations
neighborhood_size = 100000000000
cost_dist = min_cost .+ rand(divisions) .* (max_cost - min_cost)

start = time()
Random.seed!(1234)  # For reproducibility
println("Running simulated annealing with iteration limit: $iter_limit iterations, divisions: $divisions, neighborhood size: $neighborhood_size")
best_state, best_cost, visited_list = run_simulated_annealing(time_limit, cost_dist, divisions, neighborhood_size, "iter")

println("Best state: $best_state")
println("Best cost: $best_cost")

println("Real best state: $(argmin(cost_dist))")
println("Real best cost: $(minimum(cost_dist))")
println("Number of visited states: $(iter_limit) or $(iter_limit / divisions * 100)%")
