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

function is_possible(r1, r2, t)
    # Calculate the semi-major axis of the transfer orbit
    a = (r1 + r2) / 2

    # Calculate the period of the transfer orbit
    T = 2 * π * sqrt(a^3 / (G * (M_sun + M_earth)))

    # Calculate the time to reach the destination
    t_transfer = T / 2

    # Check if the transfer time is within the specified time limit
    if t_transfer > t
        return "Transfer not possible within the specified time."
    else
        return "Transfer possible in $(t_transfer / year) years."
    end
end

function circular_orbit_velocity(radius)
    return sqrt(G * M_earth / radius)
end

function elliptical_orbit_velocity(radius, semi_major_axis)
    return sqrt(G * M_earth * (2 / radius - 1 / semi_major_axis))
end

function plot_orbit(r1, r2)
    theta = LinRange(0, 2π, 100)
    x1 = r1 * cos.(theta)
    y1 = r1 * sin.(theta)
    x2 = r2 * cos.(theta)
    y2 = r2 * sin.(theta)

    plot(x1, y1, label="Initial Orbit", xlabel="x (m)", ylabel="y (m)", aspect_ratio=:equal, color=:lightgreen)
    plot!(x2, y2, label="Destination Orbit", color=:lightgreen)
    scatter!([0], [0], label="Earth", color=:blue, markersize=10)
    scatter!([r1], [0], label="Start (Moon)", color=:gray, markersize=5)
    scatter!([r2], [0], label="End (Destination)", color=:red, markersize=5)

    title!("Hohmann Transfer Orbit")
    savefig("hohmann_ideal_orbits.png")
end

function orbit(du, u, p, t)
    Mass = p[1]  # Mass of the central body (e.g., Sun)
    r = sqrt(u[1]^2 + u[2]^2)  # Distance from the central body
    if r < R_earth + R_moon
        @show r, t, u[1], u[2], u[3], u[4]
        error("The Moon is too close to the Earth!")
    end
    du[1] = u[3]  # dx/dt = vx
    du[2] = u[4]  # dy/dt = vy
    # Gravitational acceleration
    du[3] = -G * Mass * u[1] / r^3  # Gravitational force in x-direction
    du[4] = -G * Mass * u[2] / r^3  # Gravitational force in y-direction
end


function first_boost_velocity(r1, r2, v_curr)
    # Calculate the velocity needed for the first boost
    a = (r1 + r2) / 2  # Semi-major axis of the elliptical orbit
    v_peri = elliptical_orbit_velocity(r1, a)  # Velocity at periapsis
    Δv1 = v_peri - v_curr  # Required boost velocity
    @show v_peri, v_curr, Δv1
    return Δv1
end

function second_boost_velocity(r1, r2, v_curr)
    # Calculate the velocity needed for the second boost
    v_circ2 = circular_orbit_velocity(r2)  # Velocity at apoapsis
    Δv2 = v_circ2 - v_curr  # Required boost velocity
    @show v_circ2, v_curr, Δv2
    return Δv2
end

function hohmann_transfer_time(r1, r2, Mass)
    μ = G * Mass
    a = (r1 + r2) / 2
    return π * sqrt(a^3 / μ)
end

function hohmann_transfer(r1, r2, time_limit)
    # Initial position and velocity vectors
    ϕ = π / 2  # Initial angle of the satellite in radians (90 degrees)
    ψ = π / 2  # Initial angle of the destination in radians (0 degrees)

    vel1 = circular_orbit_velocity(r1)  # Velocity of the satellite in circular orbit
    vel2 = circular_orbit_velocity(r2)  # Velocity of the destination in circular orbit

    r0 = [r1 * cos(ϕ), r1 * sin(ϕ)]  # Initial position in polar coordinates]
    r_dest = [r2 * cos(ψ), r2 * sin(ψ)]  # Destination position in polar coordinates
    v0 = [vel1 * sin(ϕ), vel1 * cos(ϕ)]  # Initial velocity in circular orbit
    v_dest = [vel2 * sin(ψ), vel2 * cos(ψ)]  # Destination velocity in circular orbit
    init_state_start = [r0[1], r0[2], v0[1], v0[2]]  # Initial state vector
    init_state_dest = [r_dest[1], r_dest[2], v_dest[1], v_dest[2]]  # Destination state vector
    orbit_params = [M_earth]  # Parameters for the orbit function (mass of Earth)
    hohmann_time = hohmann_transfer_time(r1, r2, M_earth)  # Time for the Hohmann transfer

    # plot_orbit(r1, r2)

    one_period_small = 2π * r1 / vel1  # Period of the small orbit
    one_period_large = 2π * r2 / vel2  # Period of the large orbit

    @show is_possible(r1, r2, one_period_small)
    @show init_state_start, init_state_dest
    @show one_period_small, one_period_large

    plot_radius = 1.1 * r2  # 10% padding
    plt = scatter([0], [0], title="Moon Orbit around Earth", xlabel="X Position (m)", ylabel="Y Position (m)", 
        xlims=(-plot_radius, plot_radius), ylims=(-plot_radius, plot_radius), color=:blue, markersize=10, legend=false)

    # One period of the small orbit
    prob_small = ODEProblem(orbit, init_state_start, (0.0, one_period_small * 3/4), orbit_params)
    soln_small = solve(prob_small, Tsit5(), abstol=1e-8, reltol=1e-8)
    plot!(soln_small, idxs = (1, 2), color=:gray)

    # First boost
    first_boost_init_state = soln_small.u[end]  # Last state of the small Orbit
    first_satellite_velocity = (first_boost_init_state[3]^2 + first_boost_init_state[4]^2)^(1/2)  # Velocity of the satellite
    first_satellite_direction = atan(first_boost_init_state[4], first_boost_init_state[3])  # Direction of the satellite
    first_boosted_velocity = first_boost_velocity(r1, r2, first_satellite_velocity)  # Velocity needed for the elliptical orbit
    first_hohmann_boost = [ 
        0.0, 0.0,
        first_boosted_velocity * cos(first_satellite_direction),
        first_boosted_velocity * sin(first_satellite_direction),
    ]  # Boost velocity to reach the destination
    @show first_boost_init_state, first_hohmann_boost

    first_state_boost = first_boost_init_state .+ first_hohmann_boost  # Initial state after the boost
    first_prob_boost = ODEProblem(orbit, first_state_boost, (0.0, hohmann_time), orbit_params)
    first_soln_boost = solve(first_prob_boost, Tsit5(), abstol=1e-8, reltol=1e-8)
    plot!(first_soln_boost, idxs = (1, 2), color=:orange)

    # Second boost
    second_boost_init_state = first_soln_boost.u[end]  # Last state of the small Orbit
    second_satellite_velocity = (second_boost_init_state[3]^2 + second_boost_init_state[4]^2)^(1/2)  # Velocity of the satellite
    second_satellite_direction = atan(second_boost_init_state[4], second_boost_init_state[3])  # Direction of the satellite
    second_boosted_velocity = second_boost_velocity(r1, r2, second_satellite_velocity)  # Velocity needed for the elliptical orbit
    second_hohmann_boost = [ 
        0.0, 0.0,
        second_boosted_velocity * cos(second_satellite_direction),
        second_boosted_velocity * sin(second_satellite_direction),
    ]  # Boost velocity to reach the destination
    @show second_boost_init_state, second_hohmann_boost

    second_state_boost = second_boost_init_state .+ second_hohmann_boost  # Initial state after the boost
    second_prob_boost = ODEProblem(orbit, second_state_boost, (0.0, one_period_large), orbit_params)
    second_soln_boost = solve(second_prob_boost, Tsit5(), abstol=1e-8, reltol=1e-8)
    plot!(second_soln_boost, idxs = (1, 2), color=:green, alpha=0.7)

    # One period of the ideal large orbit
    prob_ideal = ODEProblem(orbit, init_state_dest, (0.0, one_period_large), orbit_params)
    soln_ideal = solve(prob_ideal, Tsit5(), abstol=1e-8, reltol=1e-8)
    plot!(soln_ideal, idxs = (1, 2), color=:red, style=:dash)

    savefig(plt, "hohmann_transfer_orbit.png")
end

hohmann_transfer(earth_moon_distance, 2 * earth_moon_distance, year)
