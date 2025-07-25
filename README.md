# Julia Learning Projects

This repository is a curated collection of Julia projects I'm building to explore computational science, algorithms, and visualization. These range from physics simulations to machine learning experiments — each designed to deepen my understanding of Julia’s strengths in performance and expressiveness.

---

## Completed Projects

### Lorenz Attractor Explorer
Visualizes the Lorenz system in 3D using Plots.jl, with a parameter sweep over σ, ρ, and β. Includes timing and animations for chaotic trajectories.

### Hohmann Transfer + Orbital Mechanics
Simulates Earth–Mars transfers using realistic ODEs and visualizes orbits in heliocentric and planet-relative frames. Models velocity corrections and orbital insertion.

### Route Plan: Genetic Algorithm for TSP
Parallel implementation of a genetic algorithm to solve the Traveling Salesman Problem for up to 250 cities. Uses multiparent crossover, mutation, and OpenMP threading.

---

## In Progress

### Fractals in Julia
A visual and interactive exploration of fractals.
- [x] Lorenz Attractor (parameter space exploration)
- [ ] Mandelbrot set (core plotting structure)
- [ ] Interactive controls (pan, zoom)
- [ ] Controlled zoom animations
- [ ] Julia sets and beyond

### Titanic Data Analysis
End-to-end analysis using the Titanic dataset with MLJ.jl. Includes feature engineering, model training, and exploratory data visualization.

---

## Planned Projects

### Spiking Neural Networks
Implement biologically inspired spiking neuron models using `snnTorch.jl` or a custom framework. Apply to audio classification using MFCC features.

### Quantum Mechanics Simulations
Solve the Schrödinger equation in 1D for systems like infinite square well and harmonic oscillator. Visualize energy eigenstates and probability densities.

### Simulated Annealing
Implement simulated annealing on benchmark optimization problems. Explore continuous and discrete variants, with tunable cooling schedules.

### COVID-19 Data Modeling
Analyze and forecast COVID-19 case trends using time-series models and differential equations. Focus on parameter inference and policy effects.

### Crime Pattern Analyzer
Use public crime datasets for geospatial and temporal analysis. Includes clustering, seasonal trends, and heatmaps.

### Portfolio Optimization
Use Convex.jl and historical financial data to build efficient frontier visualizations and risk-return balanced portfolios.

### Kaggle Challenge Replications
Solve popular Kaggle problems using Julia:
- Image classification
- Tabular classification
- NLP sentiment analysis

---
