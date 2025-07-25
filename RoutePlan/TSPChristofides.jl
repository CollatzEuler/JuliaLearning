function christofides(distances)
    # This function should return a tour that is at most 1.5 times the optimal length
    # Need triangle inequality for Christofides algorithm
    for 1:num_cities
        for j in 1:num_cities
            for k in 1:num_cities
                if distances[i, j] + distances[j, k] < distances[i, k]
                    distances[i, k] = distances[i, j] + distances[j, k]  # Update distance if shorter path found
                    distances[j, k] = 0  # Distance to self is zero
                end
            end
        end
    end
end