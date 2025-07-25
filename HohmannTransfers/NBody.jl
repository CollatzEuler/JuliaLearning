using DifferentialEquations
using Plots

# Constants
const G = 6.67430e-11
const M_sun = 1.989e30
const M_earth = 5.972e24
const M_mars = 6.4171e23
const M_sat = 1e5  # Mass of the satellite (kg)
const R_earth = 6.371e6  # Radius of Earth (m)
const R_sat = 5.0  # Radius of the satellite (m)

const sun_earth_distance = 1.496e11           # Earth-Sun distance (m)
const sun_mars_distance = 227.9e9  # Sun-Mars distance (m)
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
            
            if r_ij < R_earth + R_sat
                @show r_ij, R_earth + R_sat, t, i, j, x_i, y_i, x_j, y_j, vx_i, vy_i
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
function sun_earth_mars_orbit(sat_params, mars_params, earth_params, sun_params)
    sun_init = setup_body(sun_params)
    sun_pos = sun_init[1:2]
    sun_vel = sun_init[3:4]

    earth_init = setup_body(earth_params; origin_pos=sun_pos, origin_vel=sun_vel)
    earth_pos = earth_init[1:2]
    earth_vel = earth_init[3:4]

    mars_init = setup_body(mars_params; origin_pos=earth_pos, origin_vel=earth_vel)
    mars_pos = mars_init[1:2]
    mars_vel = mars_init[3:4]

    sat_init = setup_body(sat_params; origin_pos=earth_pos, origin_vel=earth_vel)

    total_init = vcat(sat_init, mars_init, earth_init, sun_init)
    tspan = (0.0, year * 3)

    masses = [M_sat, M_mars, M_earth, M_sun]
    prob = ODEProblem(n_body_orbit, total_init, tspan, masses)
    soln = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)

    locs = Dict(
        "Satellite" => sat_init[1:2],
        "Mars" => mars_init[1:2],
        "Earth" => earth_init[1:2],
        "Sun" => sun_init[1:2]
    )

    return soln, locs
end

# Plotting
function plot_orbits_relative(soln, locs::Dict{String, Vector{Float64}}, names::Vector{String}, reference_body::String; filename::String="")
    n = length(names)
    steps = length(soln.u)

    # Get reference trajectory
    ref_idx = findfirst(isequal(reference_body), names)
    ref_x = [soln.u[t][4ref_idx - 3] for t in 1:steps]
    ref_y = [soln.u[t][4ref_idx - 2] for t in 1:steps]

    plt = plot(title="N-body Orbit (Relative to $reference_body)", xlabel="X (m)", ylabel="Y (m)", aspect_ratio=1, legend=false)

    for (i, name) in enumerate(names)
        x = [soln.u[t][4i - 3] for t in 1:steps]
        y = [soln.u[t][4i - 2] for t in 1:steps]

        x_rel = x .- ref_x
        y_rel = y .- ref_y

        plot!(x_rel, y_rel, label=name)
        scatter!([x_rel[1]], [y_rel[1]], label="Init $name", markersize=4)
    end

    scatter!([0], [0], label="$reference_body (centered)", color=:black, markersize=5, marker=:cross)

    if filename != ""
        savefig(plt, filename)
    end

    return plt
end

sat_speed = 7.8e3  # Satellite speed (m/s)
mars_speed = sqrt(G * M_sun / sun_mars_distance)
earth_speed = sqrt(G * M_sun / sun_earth_distance)
soln, locs = sun_earth_mars_orbit([sat_speed, -pi/2 , 1e9], [mars_speed, pi/2, sun_mars_distance], [earth_speed, 0.0, sun_earth_distance], [0.0, 0.0, 0.0])  # Initial parameters for mars, Earth, and Sun
names = ["Satellite", "Mars", "Earth", "Sun"]
plot_orbits_relative(soln, locs, names, "Sun"; filename="n_body_orbit_sun.png")
plot_orbits_relative(soln, locs, names, "Earth"; filename="n_body_orbit_earth.png")
plot_orbits_relative(soln, locs, names, "Mars"; filename="n_body_orbit_mars.png")
plot_orbits_relative(soln, locs, names, "Satellite"; filename="n_body_orbit_satellite.png")
