function plot_fidelities_fourlevel(μ::Vector{Float64}, σ::Vector{Float64})
    f = Figure(size=(600, 400), fonts=(; regular="Times New Roman"))
    
    ax = Axis(
        f[1, 1],
        limits=((0, 100), (0.4, 1)),
        xticks=0:25:101,
        yticks=0.4:0.1:1,
        xminorticks=IntervalsBetween(5),
        xminorticksvisible=true,
        yminorticks=IntervalsBetween(2),
        yminorticksvisible=true,
        xgridvisible=false,
        ygridvisible=false,
        xlabel="Iterations",
        ylabel="Fidelity",
        xlabelsize=18,
        ylabelsize=18,
        xticklabelsize=16,
        yticklabelsize=16
    )

    upper, lower = μ + σ, μ - σ

    band!(ax, 1:length(μ), upper, lower, color=(:blue, 0.2), label="Standard Deviation")
    lines!(ax, 1:length(μ), μ, color=:blue, label="Average Fidelity")
    axislegend(ax, position=:rb, labelsize=18)
    
    f
end

function plot_fidelities_twolevel(μ::Vector{Float64}, σ::Vector{Float64})
    f = Figure(size=(600, 400), fonts=(; regular="Times New Roman"))
    
    ax = Axis(
        f[1, 1],
        limits=((1, 10), (0.4, 1)),
        xticks=1:1:10,
        yticks=0.4:0.1:1,
        xminorticks=IntervalsBetween(5),
        xminorticksvisible=true,
        yminorticks=IntervalsBetween(2),
        yminorticksvisible=true,
        xgridvisible=false,
        ygridvisible=false,
        xlabel="Iterations",
        ylabel="Fidelity",
        xlabelsize=18,
        ylabelsize=18,
        xticklabelsize=16,
        yticklabelsize=16
    )

    upper, lower = μ + σ, μ - σ

    band!(ax, 1:length(μ), upper, lower, color=(:blue, 0.2), label="Standard Deviation")
    lines!(ax, 1:length(μ), μ, color=:blue, label="Average Fidelity")
    axislegend(ax, position=:rb, labelsize=18)
    
    f
end


# define sphere and splines

function sphere(radius::T, center::NTuple{3, T}; n::Int=100) where {T <: Real}
    u = range(-π, π; length=n)
    v = range(0, π; length=n)
    x = center[1] .+ radius*cos.(u) * sin.(v)'
    y = center[2] .+ radius*sin.(u) * sin.(v)'
    z = center[3] .+ radius*ones(n) * cos.(v)'
    return x, y, z
end


function splines(radius::T, center::NTuple{3, T}; n::Int=100) where {T <: Real}
    u = range(-π, π; length=n)
    y, x = radius*sin.(u), radius*cos.(u)

    offseted(x::Float64, y::Float64, z::Float64) = Point3f(x + center[1], y + center[2], z + center[3])

    [offseted.(0.0, x, y)..., Point3f(NaN), offseted.(x,0.0,y)...,  Point3f(NaN), offseted.(x,y,0.0)...]
end


# define plotting function for sphere and splines

function plot_sphere(radius::T, center::NTuple{3, T}; n::Int=100) where {T <: Real}
    x, y, z = sphere(radius, center; n=n)
    f = Figure()
    ax = Axis3(f[1, 1])
    surface!(ax, sphere(1.0, (0.0, 0.0, 0.0))..., transparency=true, color=fill((:gray, 0.1), 100, 100), shading=NoShading)
    f, ax
end

function plot_sphere!(ax::Axis3, radius::T, center::NTuple{3, T}; n::Int=100) where {T <: Real}
    x, y, z = sphere(radius, center; n=n)
    surface!(ax, sphere(1.0, (0.0, 0.0, 0.0))..., transparency=true, color=fill((:gray, 0.1), 100, 100), shading=NoShading)
    ax
end

function plot_splines(radius::T=1.0, center::NTuple{3, T}=(0.0,0.0,0.0); n::Int=100) where {T <: Real}
    f = Figure()
    ax = Axis3(f[1, 1])
    lines!(ax, splines(radius, center), color=(:gray, 0.5), linewidth=1.0, overdraw=true)
    f, ax
end

function plot_splines!(ax::Axis3, radius::T, center::NTuple{3, T}; n::Int=100) where {T <: Real}
    lines!(ax, splines(radius, center), color=(:gray, 0.5), linewidth=1.0, overdraw=true)
end


# define bloch_sphere

function blochsphere(radius::T=1.0, center::NTuple{3, T}=(0.0,0.0,0.0); n::Int=100) where {T <: Real}
    f = Figure()

    ax = Axis3(
        f[1, 1],
        aspect=(1.0, 1.0, 1.0),
        xgridvisible=false,
        ygridvisible=false,
        zgridvisible=false,
        xspinesvisible=false,
        yspinesvisible=false,
        zspinesvisible=false,
        xticksvisible=false,
        yticksvisible=false,
        zticksvisible=false,
        xlabelvisible=false,
        ylabelvisible=false,
        zlabelvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        zticklabelsvisible=false,
    )

    surface!(ax, sphere(radius, center, n=n)..., transparency=true, color=fill((:gray, 0.1), 100, 100), shading=NoShading)
    lines!(ax, splines(radius, center, n=n), color=(:gray, 0.5), linewidth=1.0, overdraw=true)

    scatter!(ax, Point3f(0), color=:black, markersize=12)
    scatterlines!(ax, [Point3f(0, 0, 1), Point3f(0, 0, -1)], color=:red, markersize=12)
    scatterlines!(ax, [Point3f(0, 1, 0), Point3f(0, -1, 0)], color=:blue, markersize=12)
    scatterlines!(ax, [Point3f(1, 0, 0), Point3f(-1, 0, 0)], color=:green, markersize=12)
    text!(ax, [
        Point3f(0,0,1.3), Point3f(0,0,-1.3), 
        Point3f(0,1.3,0),Point3f(0,-1.3,0), 
        Point3f(1.3,0,0), Point3f(-1.3,0,0)],
    text = [
        L"|0 \rangle", L"|1 \rangle", 
        L"|+ \rangle", L"|- \rangle", 
        L"|R \rangle", L"|L \rangle"],
    fontsize=18,
    align=[
        (:center, :center), (:center, :center), 
        (:center, :center), (:center, :center), 
        (:center, :center), (:center, :center)]
    )

    f, ax
end

function blochvector(state::QuantumState{2})
    ρ̂ = density(state)
    return Point3f(2*real(ρ̂[1,2]), 2*imag(ρ̂[1,2]), real(ρ̂[1,1] - ρ̂[2,2]))
end

state!(ax::Axis3, state::QuantumState{2}; color::Symbol=:black) = scatterlines!(ax, [Point3f(0), blochvector(state)], color=color)
states!(ax::Axis3, states::Vector{QState}; color::Symbol=:gray) where {QState <: QuantumState{2}} = scatterlines!(ax, blochvector.(states), color=color)