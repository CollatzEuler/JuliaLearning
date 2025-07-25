using StaticArrays
using Plots
using BenchmarkTools
using Printf

function ThreeDimLorenz!(u, σ, ρ, β, dt)
    dx = σ * (u[2] - u[1])
    dy = u[1] * (ρ - u[3]) - u[2]
    dz = u[1] * u[2] - β * u[3]
    u[1] += dx * dt
    u[2] += dy * dt
    u[3] += dz * dt
end

function lorenz_plot(u0, σ, ρ, β, dt, steps, frame_step)
    u = copy(u0)
    sol = zeros(SVector{3, Float64}, steps)
    sol[1] = SVector(u...)
    for i in 1:(steps - 1)
        ThreeDimLorenz!(u, σ, ρ, β, dt)
        sol[i + 1] = SVector(u...)
    end

    xs = getindex.(sol, 1)
    ys = getindex.(sol, 2)
    zs = getindex.(sol, 3)

    # Compute bounds with 10% padding
    xlims = extrema(xs)
    ylims = extrema(ys)
    zlims = extrema(zs)

    pad = (a, b) -> (a - 0.1 * abs(b - a), b + 0.1 * abs(b - a))
    xlim = pad(xlims...)
    ylim = pad(ylims...)
    zlim = pad(zlims...)

    anim = @animate for i in 1:frame_step:steps
        plot(
            xs[1:i], ys[1:i], zs[1:i],
            color = :blue,
            linewidth = 0.5,
            xlim = xlim,
            ylim = ylim,
            zlim = zlim,
            legend = false,
            xlabel = "X", ylabel = "Y", zlabel = "Z",
            title = "σ=$(σ), ρ=$(ρ), β=$(round(β, digits=2))"
        )
    end

    filename = @sprintf "ParamPlots/lorenz_σ%.1f_ρ%.1f_β%.2f.gif" σ ρ β
    gif(anim, filename, fps = 30)
end

function sweep_parameters()
    u0 = [2.0, 1.0, 1.0]
    dt = 0.01
    steps = 5000
    frame_step = 50

    σ_vals = [10.0, 20.0, 30.0, 40.0]
    ρ_vals = [14.0, 28.0, 42.0, 56.0]
    β_vals = 2.4:0.2:3.2

    for σ in σ_vals, ρ in ρ_vals, β in β_vals
        println("Running σ=$σ, ρ=$ρ, β=$β...")
        t = @elapsed lorenz_plot(u0, σ, ρ, β, dt, steps, frame_step)
        @printf("Completed in %.2f seconds\n\n", t)
    end
end

sweep_parameters()
