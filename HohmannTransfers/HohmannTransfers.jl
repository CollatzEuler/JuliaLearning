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

function hohmann_velocities(r1, r2, μ)
    a = (r1 + r2) / 2
    v_circ1 = sqrt(μ / r1)
    v_circ2 = sqrt(μ / r2)
    v_peri = sqrt(μ * (2/r1 - 1/a))
    v_apo = sqrt(μ * (2/r2 - 1/a))

    Δv1 = v_peri - v_circ1  # Burn at r1
    Δv2 = v_circ2 - v_apo   # Burn at r2

    @show v_circ1, v_circ2, v_peri, v_apo, Δv1, Δv2
    return (v_circ1, v_circ2, v_peri, v_apo, Δv1, Δv2)
end

function hohmann_transfer(r1, r2, time_limit)
    tspan = (0.0, time_limit)

    # Initial position and velocity vectors
    ϕ = π / 2  # Initial angle of the satellite in radians (90 degrees)
    ψ = π / 2  # Initial angle of the destination in radians (0 degrees)

    vel1 = 1000#circular_orbit_velocity(r1)  # Velocity of the satellite in circular orbit
    vel2 = 500#circular_orbit_velocity(r2)  # Velocity of the destination in circular orbit

    r0 = [r1 * cos(ϕ), r1 * sin(ϕ)]  # Initial position in polar coordinates]
    r_dest = [r2 * cos(ψ), r2 * sin(ψ)]  # Destination position in polar coordinates
    v0 = [vel1 * sin(ϕ), vel1 * cos(ϕ)]  # Initial velocity in circular orbit
    v_dest = [vel2 * sin(ψ), vel2 * cos(ψ)]  # Destination velocity in circular orbit
    init_state_start = [r0[1], r0[2], v0[1], v0[2]]  # Initial state vector
    init_state_dest = [r_dest[1], r_dest[2], v_dest[1], v_dest[2]]  # Destination state vector
    orbit_params = [M_earth]  # Parameters for the orbit function (mass of Earth)

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

    # Calculate the time to reach the destination
    hohmann_init_state = soln_small.u[end]  # Last state of the small Orbit
    satellite_velocity = (hohmann_init_state[3]^2 + hohmann_init_state[4]^2)^(1/2)  # Velocity of the satellite
    satellite_direction = atan(hohmann_init_state[4], hohmann_init_state[3])  # Direction of the satellite
    boosted_velocity = elliptical_orbit_velocity(r2, (r1 + r2) / 2) - satellite_velocity  # Velocity needed for the elliptical orbit
    hohmann_boost = [ 
        0.0, 0.0,
        boosted_velocity * cos(satellite_direction),
        boosted_velocity * sin(satellite_direction),
    ]  # Boost velocity to reach the destination
    @show soln_small.u[end], hohmann_boost

    init_state_boost = hohmann_init_state .+ hohmann_boost  # Initial state after the boost
    prob_boost = ODEProblem(orbit, init_state_boost, (0.0, one_period_large / 4), orbit_params)
    soln_boost = solve(prob_boost, Tsit5(), abstol=1e-8, reltol=1e-8)
    plot!(soln_boost, idxs = (1, 2), color=:orange)

    # Period of the large orbit after the boost
    prob_large = ODEProblem(orbit, soln_boost.u[end], (0.0, one_period_large * 3/4), orbit_params)
    soln_large = solve(prob_large, Tsit5(), abstol=1e-8, reltol=1e-8)
    plot!(soln_large, idxs = (1, 2), color=:green, style=:dash)

    # One period of the ideal large orbit
    prob_ideal = ODEProblem(orbit, init_state_dest, (0.0, one_period_large), orbit_params)
    soln_ideal = solve(prob_ideal, Tsit5(), abstol=1e-8, reltol=1e-8)
    plot!(soln_ideal, idxs = (1, 2), color=:red)

    savefig(plt, "hohmann_transfer_orbit.png")
end

hohmann_transfer(earth_moon_distance, 2 * earth_moon_distance, year)
