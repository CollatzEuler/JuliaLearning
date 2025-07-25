"""
using Plots
x = range(0, 2π, length=100)
y = sin.(x)
y2 = tanh.(x)
p = plot(x, [y, y2], label=["sin(x)" "tanh(x)"], title="Sine and Tanh Functions", xlabel="x", ylabel="y", legend=:topright, grid=true)
plot!(p, x, cos.(x), label="cos(x)", color=:red) # Add cosine function to the plot
plot!(p, x, log.(x .+ 1), label="log(x+1)", color=:blue) # Add logarithmic function to the plot
savefig(p, "sine_wave.png") # Save the plot as a PNG file

p2 = plot(x, y, y2)
savefig(p2, "3D.png") # Save the 3D plot as a PNG file

x = range(1, 8π+1, length=100)
p3 = plot(x, log.(x .+ 1).^2, xlims=[1, 1e+1], ylims=[1, 1e+7], xscale=:ln, yscale=:ln, label="log(x+1)^2", title="Logarithmic Function Squared")
savefig(p3, "logarithmic_function_squared.png") # Save the logarithmic function squared plot as a PNG file

k = 3
z = range(0, 4pi, length=200)
y = @. sin(z*k) + im* cos(z*k) # Create a complex function
p4 = plot(real(y), imag(y), z, title="Complex Function", xlabel="Re(y)", ylabel="Im(y)", zlabel="x", grid=true)
y2 = @. -77sin(z*k) + im* cos(z*k) # Create a complex function
plot!(p4, real(y2), imag(y2), z, label="Conjugate Complex", color=:green) # Add negative complex part to the plot
savefig(p4, "complex_function.png") # Save the complex function plot as a PNG file

x = range(0, 20, length=100)
y = @. 2x^2 + 3x + 1 # Create a polynomial function
y_noisy = @. y + 20 * randn() # Add some noise
p5 = plot(x, y, label="Polynomial", title="Polynomial Function with Noise", xlabel="x", ylabel="y", grid=true)
scatter!(p5, x, y_noisy, label="Noisy Data", mc=:red, ms=2, ma=0.5, lw=4) # Add noisy data points to the plot
savefig(p5, "polynomial_function_with_noise.png") # Save the polynomial function plot as a PNG file

x = range(0, 50, length=200)
y1 = @. exp(-0.1x) * sin(2π*x / 10) # Create a damped sine wave
y2 = @. exp(-0.1x) * cos(2π*x / 10) # Create a damped cosine wave
p6 = plot(x, [y1, y2], label="Damped Sine Wave", title="Damped Sine and Cosine Waves", xlabel="x", ylabel="y", grid=true, layout=(2, 1))
savefig(p6, "damped_sine_cosine_waves.png") # Save the damped sine and cosine waves plot as a PNG file
"""

using StatsPlots
using DataFrames
df = DataFrame(x=1:10, y=rand(10), z=rand(10))
s = @df df scatter(:x, :y, label="Random Y", title="Scatter Plot of Random Data", xlabel="X", ylabel="Y")
savefig(s, "scatter_plot_random_data.png") # Save the scatter plot of random data as a PNG file

using Distributions
plot(Normal(3, 2), label="Normal Distribution", title="Normal Distribution Plot", xlabel="x", ylabel="Density", grid=true)
plot!(Normal(3, 4), label="Normal Distribution (σ=4)", color=:red) # Add another normal distribution with different standard deviation
plot!(Normal(3, 6), label="Normal Distribution (σ=6)", color=:green) # Add another normal distribution with different standard deviation
plot!(Normal(3, 8), label="Normal Distribution (σ=8)", color=:blue) # Add another normal distribution with different standard deviation
savefig("normal_distribution_plot.png") # Save the normal distribution plot as a PNG file
