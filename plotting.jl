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
    f = Figure()
    ax = Axis3(f[1, 1])
    surface!(ax, sphere(radius, center; n=n)..., transparency=true, color=fill((:gray, 0.1), 100, 100), shading=NoShading)
    f, ax
end

function plot_sphere!(ax::Axis3, radius::T, center::NTuple{3, T}; n::Int=100) where {T <: Real}
    surface!(ax, sphere(radius, center; n=n)..., transparency=true, color=fill((:gray, 0.1), 100, 100), shading=NoShading)
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

function plot_bases()
    f = Figure()
    ax = Axis3(f[1,1])
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

function plot_bases!(ax::Axis3, radius::T, center::NTuple{3, T}; text_spacing::Float64=0.3) where {T <: Real}
    
    text_radius = radius + text_spacing
    scatter!(ax, Point3f(center...), color=:black, markersize=12)
    
    Offset(x::Real, y::Real, z::Real) = Point3f(x + center[1], y + center[2], z + center[3])

    scatterlines!(ax, [Offset(0, 0, radius), Offset(0, 0, -radius)], color=:red, markersize=12)
    scatterlines!(ax, [Offset(0, radius, 0), Offset(0, -radius, 0)], color=:blue, markersize=12)
    scatterlines!(ax, [Offset(1, 0, 0), Offset(-1, 0, 0)], color=:green, markersize=12)
    text!(ax, [
        Offset(0,0,text_radius), Offset(0,0,-text_radius), 
        Offset(0,text_radius,0),Offset(0,-text_radius,0), 
        Offset(text_radius,0,0), Offset(-text_radius,0,0)],
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
    
    plot_sphere!(ax, radius, center, n=n)
    plot_splines!(ax, radius, center, n=n)
    plot_bases!(ax, radius, center)
    

    f, ax
end

# compute bloch vector

blochvector(ρ̂::Operator{2}) = Point3f(2*real(ρ̂[1,2]), 2*imag(ρ̂[1,2]), real(ρ̂[1,1] - ρ̂[2,2]))
blochvector(state::T) where {T <: QuantumState{2}} = blochvector(density(state))

# add state to bloch sphere

state!(ax::Axis3, state::T; color::Symbol=:black) where {T <: QuantumState{2}} = scatterlines!(ax, [Point3f(0), blochvector(state)], color=color)
states!(ax::Axis3, states::Vector{T}; color::Symbol=:gray) where {T <: QuantumState{2}} = isempty(states) || scatterlines!(ax, blochvector.(states), color=color)


struct BlochSphere{T <: QuantumState{2}}
    blochsphere::Function
    states::Vector{T}

    function BlochSphere(blochsphere::Function, states::Vector{T}) where {T <: QuantumState{2}}
        new{T}(blochsphere, states)
    end
end 

BlochSphere(;radius::Float64=1.0, center::NTuple{3, Float64}=(0.0,0.0,0.0), states::Vector{T}=Ket{2}[], n::Int=100) where {T <: QuantumState{2}} = BlochSphere(() -> blochsphere(radius, center; n=n), states)

blochsphere(bs::BlochSphere) = bs.blochsphere
states(bs::BlochSphere) = bs.states
push!(bs::BlochSphere, state::T) where {T <: QuantumState{2}} = (push!(bs.states, state); bs)
push!(bs::BlochSphere, states::T...) where {T <: QuantumState{2}} = (push!(bs.states, states...); bs)
append!(bs::BlochSphere, states::Vector{T}) where {T <: QuantumState{2}} = (append!(bs.states, states); bs)

function plot(bs::BlochSphere; plot_states::Bool=true)
    f, ax = blochsphere(bs)()
    plot_states && states!(ax, states(bs))
    f, ax
end


# plotting density matrix as three-dimensional bar chart
    
function matrixplot(ρ̂::Operator) 
    ρ̂_data = vec(data(ρ̂))
    ρ̂_real = real(ρ̂_data)
    ρ_imag = imag(ρ̂_data)
    
    f = Figure(size=(600, 900))

    xs = map(x -> x-0.5, 1:dims(ρ̂))
    ys = map(y -> y-0.5, 1:dims(ρ̂))
    zs = fill(0, dims(ρ̂), dims(ρ̂))

    kwargs = (;
        limits=((0, dims(ρ̂)), (0, dims(ρ̂)), (-1, 1)),
        xticks=xs,
        yticks=ys,
        zticks=-1:0.5:1,
        zgridvisible=false,
        xgridvisible=false,
        ygridvisible=false,
        aspect=(1.0, 1.0, 0.5),
        xlabel="",
        ylabel="",
        zlabel="",
        yreversed=true,
        xtickformat= _ -> [L"$ | %$(i-1) \rangle $" for i in 1:dims(ρ̂)],
        ytickformat= _ -> [L"$ | %$(i-1) \rangle $" for i in 1:dims(ρ̂)],
        ztickformat= zs -> [L"%$z" for z in zs],
        titlesize=28,
        xticklabelsize=24,
        yticklabelsize=24,
        zticklabelsize=24
    )

    barkwargs(data::Vector{T}) where {T <: Real} = (;
        marker=Rect(Vec3(-0.5, -0.5, 0.0), Vec3(1)),
        markersize=Vec3.(0.7, 0.7, data),
        color=data,
        colormap=:viridis,
        shading=FastShading,
    )

    axis_real = Axis3(f[1, 1]; title=L"Re(\hat{\rho})", kwargs...)
    axis_imag = Axis3(f[2, 1]; title=L"Im(\hat{\rho})", kwargs...)


    meshscatter!(axis_real, xs, ys, zs; barkwargs(ρ̂_real)...)
    meshscatter!(axis_imag, xs, ys, zs; barkwargs(ρ_imag)...)

    f, ax
end

matrixplot(state::T) where {T <: QuantumState} = matrixplot(density(state))
