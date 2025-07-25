using DifferentialEquations
using Plots

# Constants
const G = 6.67430e-11
const M_sun = 1.989e30
const M_earth = 5.972e24
const M_moon = 7.34767309e22

const AU = 1.496e11           # Earth-Sun distance (m)
const earth_moon_distance = 384400e3 *10  # Earth-Moon distance (m)
const year = 3.154e7          # One year (s)

# Newtonian N-body interaction
function n_body_orbit(du, u, p, t)
    # N-body gravitational interaction
    n = length(u) ÷ 4  # Number of bodies
    for i in 1:n
        x_i = u[4*i - 3]
        y_i = u[4*i - 2]
        vx_i = u[4*i - 1]
        vy_i = u[4*i]

        du[4*i - 3] = vx_i  # dx/dt = vx
        du[4*i - 2] = vy_i  # dy/dt = vy

        du[4*i - 1] = 0.0  # Initialize x acceleration
        du[4*i] = 0.0  # Initialize y acceleration

        for j in 1:n
            if i == j
                continue  # Skip self-interaction
            end
            x_j = u[4*j - 3]
            y_j = u[4*j - 2]
            
            dx = x_j - x_i  # Difference in x positions
            dy = y_j - y_i  # Difference in y positions
            r_ij = sqrt(dx^2 + dy^2)  # Distance between bodies i and j
            
            if r_ij < R_earth + R_moon
                @show r_ij, R_earth + R_moon, t, i, j, x_i, y_i, x_j, y_j, vx_i, vy_i
                error("Bodies are too close to each other!")
            end

            force = G * p[j] / r_ij^3  # Gravitational force magnitude
            du[4*i - 1] += force * dx  # Gravitational force in x-direction
            du[4*i] += force * dy  # Gravitational force in y-direction
        end
    end
end

function setup_body(params; origin_pos=[0.0, 0.0], origin_vel=[0.0, 0.0])
    speed, angle, dist = params  # Unpack [speed, angle, distance]
    
    # Position in global coordinates
    x = origin_pos[1] + dist * cos(angle)
    y = origin_pos[2] + dist * sin(angle)
    
    # Velocity relative to origin (tangential)
    vx = origin_vel[1] + speed * sin(angle)
    vy = origin_vel[2] + speed * cos(angle)
    
    return [x, y, vx, vy]  # Return absolute [x, y, vx, vy]
end


# Setup initial positions and velocities
function sun_earth_moon_orbit(moon_params, earth_params, sun_params)
    sun_init = setup_body(sun_params)
    sun_pos = sun_init[1:2]
    sun_vel = sun_init[3:4]

    earth_init = setup_body(earth_params; origin_pos=sun_pos, origin_vel=sun_vel)
    earth_pos = earth_init[1:2]
    earth_vel = earth_init[3:4]

    moon_init = setup_body(moon_params; origin_pos=earth_pos, origin_vel=earth_vel)

    total_init = vcat(moon_init, earth_init, sun_init)
    tspan = (0.0, year * 3)

    masses = [M_moon, M_earth, M_sun]
    prob = ODEProblem(n_body_orbit, total_init, tspan, masses)
    soln = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)

    locs = Dict(
        "Moon" => moon_init[1:2],
        "Earth" => earth_init[1:2],
        "Sun" => sun_init[1:2]
    )

    return soln, locs
end

function plot_rel_to_sun(soln, locs)
    plt = plot(soln, idxs = (1, 2), color=:gray, xlabel="X Position (m)", ylabel="Y Position (m)", title="Sun-Earth-Moon System (Sun-Centered)", legend=false)

    plot!(soln, idxs = (5, 6), color=:green, markersize=5, label="Earth Orbit", style=:dash)
    plot!(soln, idxs = (9, 10), color=:yellow, markersize=5, label="Sun Orbit")
    scatter!([locs["Sun"][1]], [locs["Sun"][2]], color=:yellow, markersize=10, label="Sun", marker=:circle)
    scatter!([locs["Earth"][1]], [locs["Earth"][2]], color=:blue, markersize=5, label="Initial Earth Position", marker=:circle)
    scatter!([locs["Moon"][1]], [locs["Moon"][2]], color=:red, markersize=5, label="Initial Moon Position", marker=:circle)
    savefig(plt, "3_body_orbit_rel_sun.png")

    return plot
end

function plot_rel_to_earth(soln, locs)
    # Extract full trajectories
    moon_x = [u[1] for u in soln.u]
    moon_y = [u[2] for u in soln.u]
    earth_x = [u[5] for u in soln.u]
    earth_y = [u[6] for u in soln.u]
    sun_x = [u[9] for u in soln.u]
    sun_y = [u[10] for u in soln.u]

    # Shift all positions relative to Earth
    moon_x_rel = moon_x .- earth_x
    moon_y_rel = moon_y .- earth_y
    earth_x_rel = zeros(length(earth_x))  # Earth is at origin
    earth_y_rel = zeros(length(earth_y))
    sun_x_rel = sun_x .- earth_x
    sun_y_rel = sun_y .- earth_y

    plt = plot(xlabel="X (m, relative to Earth)", ylabel="Y (m, relative to Earth)", title="Sun-Earth-Moon System (Earth-Centered)", legend=false)

    plot!(moon_x_rel, moon_y_rel, label="Moon", color=:gray, markersize=0.1)
    plot!(earth_x_rel, earth_y_rel, label="Earth", color=:blue, markersize=0.1)
    plot!(sun_x_rel, sun_y_rel, label="Sun", color=:orange, markersize=3)
    # scatter!([0], [0], markersize=5, color=:blue, label="Earth (centered)")
    savefig(plt, "3_body_orbit_rel_earth.png")

    return plot
end

moon_speed = sqrt(G * M_earth / earth_moon_distance)
earth_speed = sqrt(G * M_sun / AU)
soln, locs = sun_earth_moon_orbit([moon_speed, pi/2, earth_moon_distance], [earth_speed, 0.0, AU], [30000.0, 0.0, 0.0])  # Initial parameters for Moon, Earth, and Sun

plot_rel_to_earth(soln, locs)
plot_rel_to_sun(soln, locs)