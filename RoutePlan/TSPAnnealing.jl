using Distributions
using Random
using Printf

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

function travelling_salesman(iter_limit, num_cities, distances)
    function temp_function(r)
        sigmoid(z::Real) = 1.0 / (1.0 + exp(-z))
        return 1 - sigmoid(r)
    end

    function prob_function(current_cost, neighbor_cost, temperature)
        if neighbor_cost < current_cost
            return 1.0
        else
            return exp((current_cost - neighbor_cost) / temperature)
        end
    end

    function neighbors_function(state)
        neighbors = []
        for i in 1:(length(state)-1)
            for j in (i+1):length(state)
                new_state = copy(state)
                new_state[i], new_state[j] = new_state[j], new_state[i]
                push!(neighbors, new_state)
            end
        end
        return neighbors
    end

    function cost_function(state)
        return sum(distances[state[i], state[i+1]] for i in 1:(length(state)-1)) + distances[state[end], state[1]]
    end

    function init_state()
        return shuffle(1:num_cities)
    end

    init = init_state()
    start = time()
    best_state, best_cost, _ = simulated_annealing_iter(init, iter_limit, prob_function, cost_function, temp_function, neighbors_function)
    println("Took $(time() - start) seconds")

    return best_state, best_cost
end

iter_limit = 1000000
num_cities = 30
min_distance = 1  # Minimum distance between cities
max_distance = 100  # Maximum distance between cities

Random.seed!(123)  # For reproducibility
distances = min_distance .+ rand(num_cities, num_cities) .* (max_distance - min_distance)
for i in 1:num_cities
    distances[i, i] = 0  # Distance to self is zero
end
for i in 1:num_cities
    for j in (i+1):num_cities
        distances[j, i] = distances[i, j]  # Ensure symmetry
    end
end

best_state, best_cost = travelling_salesman(iter_limit, num_cities, distances)
println("Best state: $best_state")
println("Best cost: $best_cost")

distances_to_string(distances) = join(
    [join([@sprintf("%.2f", distances[i, j]) for j in 1:size(distances, 2)], " ") for i in 1:size(distances, 1)],
    "\n"
)

println("Distances matrix:\n$(distances_to_string(distances))")