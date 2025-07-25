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

function circular_orbit_velocity(radius)
    return sqrt(G * M_earth / radius)
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

function hohmann_transfer(r1, r2, time_limit)

    # Initial position and velocity vectors
    ϕ = π / 2 # Initial angle of the satellite in radians (90 degrees)
    ψ = π / 2 # Initial angle of the destination in radians (0 degrees)

    vel1 = circular_orbit_velocity(r1)  # Velocity of the satellite in circular orbit
    vel2 = 2000#circular_orbit_velocity(r2)  # Velocity of the destination in circular orbit

    r0 = [r1 * cos(ϕ), r1 * sin(ϕ)]  # Initial position in polar coordinates]
    r_dest = [r2 * cos(ψ), r2 * sin(ψ)]  # Destination position in polar coordinates
    v0 = [vel1 * sin(ϕ), vel1 * cos(ϕ)]  # Initial velocity in circular orbit
    v_dest = [vel2 * sin(ψ), vel2 * cos(ψ)]  # Destination velocity in circular orbit
    init_state_start = [r0[1], r0[2], v0[1], v0[2]]  # Initial state vector
    init_state_dest = [r_dest[1], r_dest[2], v_dest[1], v_dest[2]]  # Destination state vector
    orbit_params = [M_earth]  # Parameters for the orbit function (mass of Earth)

    one_period_small = 2π * r1 / vel1  # Period of the small orbit
    one_period_large = 2π * r2 / vel2  # Period of the large orbit

    @show is_possible(r1, r2, one_period_small)
    @show init_state_start, init_state_dest
    @show one_period_small, one_period_large

    plt = scatter([0], [0], title="Moon Orbit around Earth", xlabel="X Position (m)", ylabel="Y Position (m)", color=:blue, markersize=10, legend=false)

    # One period of the small orbit
    prob_small = ODEProblem(orbit, init_state_start, (0.0, one_period_small * 3/4), orbit_params)
    soln_small = solve(prob_small, Tsit5(), abstol=1e-8, reltol=1e-8)
    plot!(soln_small, idxs = (1, 2), color=:gray)

    savefig(plt, "hohmann_transfer_orbit.png")
end

hohmann_transfer(earth_moon_distance, 2 * earth_moon_distance, year)
