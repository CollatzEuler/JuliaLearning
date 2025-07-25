using Distributions
using Random
using Printf
using StatsBase
using Base.Threads

function GA(fitness_function, crossover, mutator, parents_per_child, mutation_rate, population_size, time_limit, num_cities, distances)
    # Initial population with cached fitness
    population = [(randperm(num_cities), 0.0) for _ in 1:population_size]
    for i in 1:population_size
        population[i] = (population[i][1], fitness_function(population[i][1], distances))
    end

    best_state = population[1][1]
    best_cost = population[1][2]
    iters = 0
    time_start = time()

    while time() - time_start < time_limit
        # Preallocate thread-safe buffer
        new_population = Vector{Tuple{Vector{Int}, Float64}}(undef, population_size)

        Threads.@threads for i in 1:population_size
            parent_ids = sample(1:population_size, parents_per_child, replace=false)
            parents = [population[j][1] for j in parent_ids]
            child = crossover(parents, num_cities)
            if rand() < mutation_rate
                child = mutator(child)
            end
            fitness = fitness_function(child, distances)
            new_population[i] = (child, fitness)
        end

        # Combine and keep best individuals
        combined_population = vcat(population, new_population)
        survivors = partialsort(combined_population, 1:population_size, by = x -> x[2])
        population = survivors

        # Track best state
        if survivors[1][2] < best_cost
            best_state = survivors[1][1]
            best_cost = survivors[1][2]
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

function crossover_ox_multiparent(parents::Vector{Vector{Int}}, num_cities::Int)
    num_parents = length(parents)
    parent1 = parents[1]
    others = parents[2:end]

    # Step 1: Random slice from parent1
    start = rand(1:num_cities)
    stop = rand(start:num_cities)

    child = fill(-1, num_cities)
    used = falses(num_cities)

    # Copy slice from parent1
    for i in start:stop
        val = parent1[i]
        child[i] = val
        used[val] = true
    end

    # Step 2: Fill the rest using round-robin from other parents
    p_idx = 1
    pos = 1
    while pos <= num_cities
        if child[pos] != -1
            pos += 1
            continue
        end
        for attempt in 1:(num_parents - 1)
            p = others[p_idx]
            for val in p
                if !used[val]
                    child[pos] = val
                    used[val] = true
                    pos += 1
                    break
                end
            end
            p_idx = mod1(p_idx + 1, num_parents - 1)
            if child[pos - 1] != -1
                break
            end
        end
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

    start = time()
    best_state, best_cost, iters = GA(fitness, crossover_ox_multiparent, mutate_swap, parents_per_child, mutation_rate, population_size, time_limit, num_cities, distances)
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
end

time_limit = 90.0  # Time limit in seconds
population_size = 200  # Size of the population
num_cities = 250
min_distance = 1  # Minimum distance between cities
max_distance = 100  # Maximum distance between cities
parents_per_child = 5 # Number of parents to select for crossover
mutation_rate = 0.1  # Probability of mutation

main(time_limit, population_size, num_cities, min_distance, max_distance, parents_per_child, mutation_rate)