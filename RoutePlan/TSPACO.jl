using Random
using Printf
using Distributions

function get_next_city(pheremones, distances, current_city, visited, α, β)
    num_cities = size(pheremones, 1)
    probabilities = zeros(num_cities)

    for city in 1:num_cities
        if city in visited
            probabilities[city] = 0.0
        else
            pheromone = pheremones[current_city, city]
            heuristic = 1 / distances[current_city, city]
            probabilities[city] = pheromone^α * heuristic^β
        end
    end

    probabilities /= sum(probabilities)  # Normalize probabilities
    return rand(Categorical(probabilities))  # Select next city based on probabilities
end

function construct_tour(pheremones, distances, num_cities, α, β)
    visited = Int[]
    current_city = 1 # Doesn't matter which city to start from
    push!(visited, current_city)

    while length(visited) < num_cities
        next_city = get_next_city(pheremones, distances, current_city, visited, α, β)
        push!(visited, next_city)
        current_city = next_city
    end

    return visited
end

function calculate_cost(state, distances)
    cost = 0.0
    for i in 1:(length(state) - 1)
        cost += distances[state[i], state[i + 1]]
    end
    cost += distances[state[end], state[1]]  # Return to the starting city
    return cost
end

function ACO(α, β, ρ, Q, time_limit, num_ants, num_cities, distances)
    pheremones = ones(num_cities, num_cities)  # Initialize pheromone levels

    start = time()
    iters = 0

    best_state = Int[]
    best_cost = Inf
    while time() - start < time_limit
        ants = []

        for i in 1:num_ants
            state = construct_tour(pheremones, distances, num_cities, α, β) # Construct a tour for each ant
            cost = calculate_cost(state, distances)  # Complete the tour
            push!(ants, (state, cost))
        end

        # Evaporate pheromones
        pheremones *= (1 - ρ)

        # Update pheromone levels based on the ants' tours
        for (state, cost) in ants
            for i in 1:num_cities
                for j in (i+1):num_cities
                    pheremones[state[i], state[j]] += Q / cost
                    pheremones[state[j], state[i]] += Q / cost  # Ensure symmetry
                end
            end
        end

        for (state, cost) in ants
            if iters == 0 || cost < best_cost
                best_state = state
                best_cost = cost
            end
        end

        iters += 1
    end

    return best_state, best_cost, iters
end

time_limit = 30.0  # Time limit in seconds
num_cities = 250
num_ants = 40
α = 1.0  # Pheromone importance
β = 2.0  # Heuristic importance
ρ = 0.1  # Evaporation rate
Q = 100.0  # Pheromone deposit factor
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

start = time()
best_state, best_cost, iters = ACO(α, β, ρ, Q, time_limit, num_ants, num_cities, distances)
elapsed = time() - start

println("Elapsed time: ", elapsed, " seconds")
println("Best tour found: ", join(best_state, ", "))
@printf("Best cost: %.2f\n", best_cost)
println("Ran for $iters iterations")