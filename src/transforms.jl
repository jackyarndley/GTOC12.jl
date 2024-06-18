using LinearAlgebra

function classical_to_cartesian(classical; μ = 1.0)
    a = classical[1, :]
    e = classical[2, :]
    i = classical[3, :]
    Ω = classical[4, :]
    ω = classical[5, :]
    M = classical[6, :]

    v = copy(M)

    for _ in 1:10
        v .-= (v .- e.*sin.(v) .- M)./(1 .- e.*cos.(v))
    end

    # Compute true anomaly v
    v = 2*atan.(sqrt.((1.0 .+ e)./(1.0 .- e)).*tan.(v./2.0))

    i[i .< 0] .+= π

    p = a.*(1 .- e.^2)

    cosω = cos.(ω)
    cosΩ = cos.(Ω)
    cosi = cos.(i)
    cosv = cos.(v)
    sinω = sin.(ω)
    sinΩ = sin.(Ω)
    sini = sin.(i)
    sinv = sin.(v)

    temp1 = [cosω.*cosΩ - sinω.*sinΩ.*cosi cosω.*sinΩ + sinω.*cosΩ.*cosi sinω.*sini]
    temp2 = [-sinω.*cosΩ - cosω.*sinΩ.*cosi -sinω.*sinΩ + cosω.*cosΩ.*cosi cosω.*sini]

    r = p./(1 .+ e.*cos.(v)) .* (temp1.*cosv .+ temp2.*sinv)
    v = sqrt.(μ./p) .* (-temp1.*sinv .+ temp2.*(e .+ cosv))

    return hcat(r, v)'
end

function cartesian_to_classical(cartesian; μ = 1.0)
    r = cartesian[1:3]
    v = cartesian[4:6]

    h = cross(r, v)
    n = cross([0, 0, 1], h)

    r_dot_v = dot(r, v)

    r_norm = norm(r)
    v_norm = norm(v)
    n_norm = norm(n)
    h_norm = norm(h)

    e_vec = ((v_norm^2 - μ/r_norm)*r - r_dot_v*v) / μ
    e = norm(e_vec)

    a = -μ/(v_norm^2 - 2*μ/r_norm)

    i = acos(h[3]/h_norm)
    Ω = acos(n[1]/n_norm)
    ω = acos(dot(n, e_vec)/(e*n_norm))
    ν = acos(clamp(dot(e_vec, r)/(e*r_norm), -1.0, 1.0))

    if n[2] < 0
        Ω = 2π - Ω
    end

    if e_vec[3] < 0
        ω = 2π - ω
    end

    if r_dot_v < 0
        ν = 2π - ν
    end

    M = atan(-sqrt(1-e^2)*sin(ν), -e-cos(ν)) + π - e*(sqrt(1-e^2)*sin(ν))/(1+e*cos(ν))
    
    return [a, e, i, Ω, ω, M]
end

function cartesian_to_spherical(cartesian)
    x, y, z, ẋ, ẏ, ż = cartesian

    r = sqrt(x^2 + y^2 + z^2)
    θ = atan(y, x)
    ϕ = asin(z / r)

    r_dot = (x*ẋ + y*ẏ + z*ż) / r
    θ_dot = (x*ẏ - ẋ*y)/(x^2 + y^2)
    ϕ_dot = (ż*x^2 - ẋ*z*x + ż*y^2 - ẏ*z*y) / (r^2 * sqrt(x^2 + y^2))

    return [r, θ, ϕ, r_dot, θ_dot, ϕ_dot]
end

function spherical_to_cartesian(spherical)
    r, θ, ϕ, r_dot, θ_dot, ϕ_dot = spherical

    sin_θ, cos_θ = sincos(θ)
    sin_ϕ, cos_ϕ = sincos(ϕ)
    
    x = r*cos_ϕ*cos_θ
    y = r*cos_ϕ*sin_θ
    z = r*sin_ϕ

    ẋ = r_dot*cos_ϕ*cos_θ - r*θ_dot*cos_ϕ*sin_θ - ϕ_dot*r*sin_ϕ*cos_θ
    ẏ = r_dot*cos_ϕ*sin_θ - ϕ_dot*r*sin_ϕ*sin_θ + r*θ_dot*cos_ϕ*cos_θ
    ż = r_dot*sin_ϕ + ϕ_dot*r*cos_ϕ

    return [x, y, z, ẋ, ẏ, ż]
end