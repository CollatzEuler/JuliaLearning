using Distributions
using Random
using Printf
using StatsBase
using Profile
using ProfileView
using Base.Threads

function GA(fitness_function, crossover, mutator, parents_per_child, mutation_rate, population_size, time_limit, num_cities, distances)
    population = [randperm(num_cities) for _ in 1:population_size]  # Initialize population with random permutations
    best_state = population[1]
    best_cost = fitness_function(best_state, distances)

    iters = 0

    time_start = time()

    while time() - time_start < time_limit
        new_population = []
        for i in 1:population_size
            parent_ids = sample(1:population_size, parents_per_child, replace=false)
            parents = population[parent_ids]
            child = crossover(parents, num_cities)
            push!(new_population, child)
        end

        population = vcat(population, new_population)
        for i in 1:length(population)
            if rand() < mutation_rate
                population[i] = mutator(population[i])  # Mutate with probability mutation_rate
            end
        end

        population = partialsort(population, 1:population_size, by = x -> fitness_function(x, distances)) # Sort by fitness and only keep the best individuals

        if fitness_function(population[1], distances) < best_cost
            best_state = population[1]
            best_cost = fitness_function(best_state, distances)
        end
        iters += 1
    end

    return best_state, best_cost, iters
end

function fitness(state, distances)
    cost = 0.0
    # @show state
    for i in 1:(length(state) - 1)
        cost += distances[state[i], state[i + 1]]
    end
    cost += distances[state[end], state[1]]  # Return to the starting city
    return cost
end

# --- Partition-Based Crossover ---
function crossover_partition(parents, num_cities)
    num_parents = length(parents)
    partition = randperm(num_cities - 1)[1:num_parents - 1]
    partition = vcat(1, partition, num_cities)

    child = fill(-1, num_cities)
    assigned = falses(num_cities)

    # Step 1: Fill each segment from a different parent without overwriting
    for i in 1:num_parents
        parent = parents[i]
        seg_start = partition[i]
        seg_end = partition[i + 1]
        for j in seg_start:seg_end
            val = parent[j]
            if !assigned[val] && child[j] == -1  # <- prevent overwriting
                child[j] = val
                assigned[val] = true
            end
        end
    end

    # Step 2: Fill remaining -1s with unused cities
    unused = Int[]
    remaining_slots = Int[]

    for i in 1:num_cities
        if child[i] == -1
            push!(remaining_slots, i)
        end
        if !assigned[i]
            push!(unused, i)
        end
    end

    for (i, idx) in enumerate(remaining_slots)
        child[idx] = unused[i]
    end

    return child
end

function mutate_swap(state)
    mutated = copy(state)  # Create a copy to avoid modifying the original state
    idx1, idx2 = sample(1:length(state), 2, replace=false)  # Select two distinct indices
    mutated[idx1], mutated[idx2] = mutated[idx2], mutated[idx1]  # Swap the values at these indices
    return mutated
end

function main(time_limit, population_size, num_cities, min_distance, max_distance, parents_per_child, mutation_rate)
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

    GA(fitness, crossover_partition, mutate_swap, parents_per_child, mutation_rate, population_size, 0.5, num_cities, distances)
    Profile.clear()

    start = time()
    @profile best_state, best_cost, iters = GA(fitness, crossover_partition, mutate_swap, parents_per_child, mutation_rate, population_size, time_limit, num_cities, distances)
    elapsed = time() - start

    start_idx = findfirst(x -> x == 1, best_state)
    if start_idx === nothing
        error("City 1 not found in the best state")
    elseif start_idx == 1
        best_state = best_state  # Already starts with city 1
    else
        best_state = vcat(best_state[start_idx:end], best_state[1:start_idx-1])  # Rotate to start from city 1
    end

    println("Elapsed time: ", elapsed, " seconds")
    println("Best tour found: ", join(best_state, ", "))
    @printf("Best cost: %.2f\n", best_cost)
    println("Ran for $iters iterations")

    ProfileView.view()
end

time_limit = 30.0  # Time limit in seconds
population_size = 200  # Size of the population
num_cities = 250
min_distance = 1  # Minimum distance between cities
max_distance = 100  # Maximum distance between cities
parents_per_child = 3 # Number of parents to select for crossover
mutation_rate = 0.2  # Probability of mutation

main(time_limit, population_size, num_cities, min_distance, max_distance, parents_per_child, mutation_rate)