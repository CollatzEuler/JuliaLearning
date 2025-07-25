using DifferentialEquations
using Plots

const M_sun = 1.989e30  # Mass of the Sun in kg
const M_earth = 5.972e24  # Mass of the Earth in kg
const M_moon = 7.34767309e22  # Mass of the Moon in kg
const R_earth = 6.371e6  # Radius of the Earth in meters
const R_moon = 1.7374e6  # Radius of the Moon in meters
const AU = 1.496e11    # Astronomical unit in meters
const earth_moon_distance = 384400e3  # Average distance from Earth to Moon in meters
const G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
const year = 3.154e7  # One year in seconds

function n_body_orbit(du, u, p, t)
    # N-body gravitational interaction
    n = length(u) ÷ 4  # Number of bodies
    for i in 1:n
        x_i = u[4*i - 3]
        y_i = u[4*i - 2]
        vx_i = u[4*i - 1]
        vy_i = u[4*i]
        r_i = sqrt(x_i^2 + y_i^2)  # Distance from the origin
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

            force = G * p[i] * p[j] / r_ij^3  # Gravitational force magnitude
            du[4*i - 1] += force * dx  # Gravitational force in x-direction
            du[4*i] += force * dy  # Gravitational force in y-direction
        end
    end
end

function sun_earth_moon_orbit(moon_params, earth_params, sun_params) # params = [speed, angle]
    
    sun_init_speed = sun_params[1]  # The Sun is stationary in this model
    sun_init_x = 0.0  # Initial x position of the Sun
    sun_init_y = 0.0  # Initial y position of the Sun
    sun_init_vx = 0.0  # Initial x velocity of the Sun
    sun_init_vy = 0.0  # Initial y velocity of the Sun
    sun_init = [sun_init_x, sun_init_y, sun_init_vx, sun_init_vy]  # Initial state vector of the Sun

    earth_init_speed = earth_params[1]  # Initial speed of the Earth in m/s
    earth_init_angle = earth_params[2]  # Initial angle of the Earth in radians (0 degrees)
    earth_init_x = sun_init_x + AU * cos(earth_init_angle)  # Initial x position of the Earth (1 AU from the Sun)
    earth_init_y = sun_init_y + AU * sin(earth_init_angle)  # Initial y position of the Earth (0 AU)
    earth_init_vx = sun_init_vx + earth_init_speed * sin(earth_init_angle)  # Initial x velocity of the Earth
    earth_init_vy = sun_init_vy + earth_init_speed * cos(earth_init_angle)  # Initial y velocity of the Earth
    earth_init = [earth_init_x, earth_init_y, earth_init_vx, earth_init_vy]  # Initial state vector of the Earth

    moon_init_speed = moon_params[1] # Initial speed of the Moon in m/s
    moon_init_angle = moon_params[2]  # Initial angle of the Moon in radians (90 degrees)
    moon_init_x = earth_init_x + earth_moon_distance * cos(moon_init_angle)  # Initial x position of the Moon
    moon_init_y = earth_init_y + earth_moon_distance * sin(moon_init_angle)  # Initial y position of the Moon
    moon_init_vx = earth_init_vx + moon_init_speed * sin(moon_init_angle)  # Initial x velocity of the Moon
    moon_init_vy = earth_init_vy + moon_init_speed * cos(moon_init_angle)  # Initial y velocity of the Moon
    moon_init = [moon_init_x, moon_init_y, moon_init_vx, moon_init_vy]  # Initial state vector of the Moon

    total_init = vcat(moon_init, earth_init, sun_init)  # Combine Moon, Earth, and Sun initial states

    tspan = (0.0, year* 30/365)  # Time span for one year

    prob = ODEProblem(n_body_orbit, total_init, tspan, [M_moon, M_earth, M_sun])
    soln = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)

    return soln
end

function lagrange_points(mu)
    # Calculate the Lagrange points for a two-body system
    
end

soln = sun_earth_moon_orbit([9000, π/2], [29780, 0π], [0, 0])

plot(soln, idxs = (1, 2), xlabel="X Position (m)", ylabel="Y Position (m)", title="Three Body Orbit", legend=false)
plot!(soln, idxs = (1, 2), color=:blue, markersize=5, label="Moon Orbit")
plot!(soln, idxs = (5, 6), color=:green, markersize=5, label="Earth Orbit")
plot!(soln, idxs = (9, 10), color=:yellow, markersize=5, label="Sun Orbit")
# scatter!([0], [0], color=:yellow, markersize=10, label="Sun", marker=:circle)
# scatter!([AU], [0], color=:blue, markersize=5, label="Initial Earth Position", marker=:circle)
# scatter!([earth_moon_distance * cos(π / 2)], [earth_moon_distance * sin(π / 2)], color=:red, markersize=5, label="Initial Moon Position", marker=:circle)

savefig("3_body_orbit.png")

