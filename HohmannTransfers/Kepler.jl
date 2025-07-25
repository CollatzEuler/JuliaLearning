using DifferentialEquations
using Plots
# using Unitful #NOTE THIS DOES NOT WORK DUE TO A BUG IN Unitful WITH DifferentialEquations
# using CairoMakie


const M_sun = 1.989e30  # Mass of the Sun in kg
const M_earth = 5.972e24  # Mass of the Earth in kg
const M_moon = 7.34767309e22  # Mass of the Moon in kg
const R_earth = 6.371e6  # Radius of the Earth in meters
const R_moon = 1.7374e6  # Radius of the Moon in meters
const AU = 1.496e11    # Astronomical unit in meters
const earth_moon_distance = 384400e3  # Average distance from Earth to Moon in meters
const G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
const year = 3.154e7  # One year in seconds

function orbit(du, u, p, t)
    Mass = p[1]  # Mass of the central body (e.g., Sun)
    r = sqrt(u[1]^2 + u[2]^2)  # Distance from the central body
    if r < R_earth + R_moon
        error("The Moon is too close to the Earth!")
    end
    du[1] = u[3]  # dx/dt = vx
    du[2] = u[4]  # dy/dt = vy
    # Gravitational acceleration
    du[3] = -G * Mass * u[1] / r^3  # Gravitational force in x-direction
    du[4] = -G * Mass * u[2] / r^3  # Gravitational force in y-direction
end

moon_init_speed = 900 # Initial speed of the Moon in m/s
moon_init_angle = π / 2  # Initial angle of the Moon in radians (90 degrees)
moon_init_x = earth_moon_distance * cos(moon_init_angle)  # Initial x position of the Moon
moon_init_y = earth_moon_distance * sin(moon_init_angle)  # Initial y position of the Moon
moon_init_vx = -moon_init_speed * sin(moon_init_angle)  # Initial x velocity of the Moon
moon_init_vy = moon_init_speed * cos(moon_init_angle)  # Initial y velocity of the Moon
moon_init = [moon_init_x, moon_init_y, moon_init_vx, moon_init_vy]  # Initial state vector of the Moon

tspan = (0.0, year)  # Time span for one year

prob = ODEProblem(orbit, moon_init, tspan, [M_earth])
soln = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)

# Using Plots for a static plot
# plot(soln, idxs = (1, 2), xlabel="X Position (m)", ylabel="Y Position (m)", title="Moon Orbit around Earth", legend=:topright)
# savefig("moon_orbit.png")

# Using CairoMakie for a more interactive plot
# fig = Figure(size = (800, 600))
# ax = Axis(fig[1, 1], xlabel = "X Position (m)", ylabel = "Y Position (m)", title = "Moon Orbit around Earth")
# lines!(ax, soln[1, :], soln[2, :], color = :gray, label = "Moon Orbit")
# scatter!(ax, [0], [0], color = :blue, markersize = 10, label = "Earth", marker = :circle)
# scatter!(ax, [moon_init_x], [moon_init_y], color = :red, markersize = 5, label = "Initial Moon Position", marker = :circle)
# # axislegend(ax, position = :topright)
# save("moon_orbit_cairomakie.png", fig)

# Animated plot
nframes = round(Int64, length(soln.t)/10) 
framerate = 30  # Frames per second
loc_iterator = range(0, stop = nframes - 1, length = nframes)
@show Int(nframes)
ENV["GKSwstype"] = "100"

function time_change(time)
    # Calculate the current position of the Moon
    x = soln.u[Int(time + 1)][1]
    y = soln.u[Int(time + 1)][2]
    return (x, y)
end


anim = @animate for i in loc_iterator
    if i % 100 == 0
        println("Frame: ", i)
    end
    x, y = time_change(i)
    plt = plot(soln, idxs = (1, 2), color = :green, xlabel = "X Position (m)", ylabel = "Y Position (m)", title = "Moon Orbit around Earth", legend=false)
    scatter!(plt, [0], [0], color = :blue, markersize = 10, legend=false, marker = :circle)
    scatter!(plt, [x], [y], color = :gray, markersize = 3, legend=false, marker = :circle)
end
gif(anim, "moon_orbit_animation.gif", fps = framerate)
