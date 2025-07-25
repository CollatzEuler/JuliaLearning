using Random

# --- Original Crossover (PMX-like) ---
function crossover_pmx(parent1, parent2)
    size = length(parent1)
    start = rand(1:size)
    stop = rand(start:size)
    child = fill(-1, size)

    child[start:stop] = parent1[start:stop]
    @show "Crossover segment: ", child[start:stop], start, stop

    idx = 1
    for i in 1:size
        if child[i] == -1
            while parent2[idx] in child
                idx += 1
            end
            child[i] = parent2[idx]
        end
    end

    return child
end

# --- Fixed Partition-Based Crossover ---
function crossover_partition(parents, num_cities)
    num_parents = length(parents)
    partition = sort(randperm(num_cities - 1)[1:num_parents - 1])
    partition = vcat(1, partition, num_cities)

    child = fill(-1, num_cities)
    assigned = Set{Int}()
    assigned_positions = Set{Int}()

    # Step 1: Assign segments, skip duplicates
    for i in 1:num_parents
        parent = parents[i]
        seg_start = partition[i]
        seg_end = partition[i + 1]
        for j in seg_start:seg_end
            val = parent[j]
            if !(val in assigned) && !(j in assigned_positions)
                child[j] = val
                push!(assigned, val)
                push!(assigned_positions, j)
            end
        end
    end

    # Step 2: Fill remaining positions
    remaining = collect(setdiff(1:num_cities, assigned))
    missing_indices = findall(==( -1), child)

    # Safe fallback: only assign up to min(len(remaining), len(missing))
    fill_count = min(length(missing_indices), length(remaining))

    for i in 1:fill_count
        child[missing_indices[i]] = remaining[i]
    end

    # Final pass: if any positions are still -1, randomly fill with unused values
    if any(x -> x == -1, child)
        unused = setdiff(1:num_cities, Set(child))
        rest = findall(x -> x == -1, child)
        for (i, idx) in enumerate(rest)
            child[idx] = unused[i]
        end
    end

    return child
end


# --- Test Runner ---
function test_crossovers()
    num_cities = 10
    parent1 = shuffle(1:num_cities)
    parent2 = shuffle(1:num_cities)
    parent3 = shuffle(1:num_cities)

    println("Parent 1: ", parent1)
    println("Parent 2: ", parent2)
    println("Parent 3: ", parent3)

    println("\n--- PMX-Like Crossover ---")
    child_pmx = crossover_pmx(parent1, parent2)
    println("Child:    ", child_pmx)

    println("\n--- Partition-Based Crossover ---")
    child_part = crossover_partition([parent1, parent2, parent3], num_cities)
    println("Child:    ", child_part)
end

# Run the test
test_crossovers()
