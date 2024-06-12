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

function spherical_to_orbital_koopman(spherical; length_unit=1.0, μ=1.0, extra=false)
    r, θ, ϕ, r_dot, θ_dot, ϕ_dot = spherical

    p_r = r_dot
    p_θ = r^2 * θ_dot * cos(ϕ)^2
    p_ϕ = r^2 * ϕ_dot

    p_a = sqrt(p_ϕ^2 + p_θ^2 / cos(ϕ)^2)

    Λ = (p_a / r - (μ / p_a)) * sqrt(length_unit / μ)
    η = p_r * sqrt(length_unit / μ)
    s = sin(ϕ)
    γ = (p_ϕ / p_a) * cos(ϕ)
    κ = (1 / p_a) * sqrt(μ * length_unit)

    if p_ϕ >= 0.0
        β = θ - asin(tan(ϕ) * sqrt(p_θ^2 / (p_a^2 - p_θ^2)))
    else
        β = θ + asin(tan(ϕ) * sqrt(p_θ^2 / (p_a^2 - p_θ^2))) + π
    end

    if extra
        χ = p_θ*sqrt(μ*length_unit)^3 / (p_a^4*(s^2 + γ^2))
        ρ = p_θ/p_a

        return [Λ, η, s, γ, κ, β, χ, ρ]
    else
        return [Λ, η, s, γ, κ, β]
    end
end

function orbital_koopman_to_spherical(orbital_koopman; length_unit=1.0, μ=1.0, prograde=true)
    Λ, η, s, γ, κ, β = orbital_koopman

    p_a = sqrt(μ * length_unit) / κ

    r = p_a/(Λ*sqrt(μ/length_unit) + μ/p_a)

    p_r = sqrt(μ / length_unit) * η

    ϕ = asin(s)
    p_ϕ = γ * p_a / cos(ϕ)

    si = sqrt(s^2 + γ^2)

    if p_ϕ >= 0
        u = asin(s / si)
    else
        u = π - asin(s / si)
    end

    # if (s >= 0 && prograde) || (s < 0 && !prograde)
    if s >= 0
        θ = β + acos(cos(u) / cos(ϕ))
    else
        θ = β - acos(cos(u) / cos(ϕ))
    end

    if prograde
        p_θ = sqrt((p_a^2 - p_ϕ^2) * cos(ϕ)^2) 
    else
        p_θ = -sqrt((p_a^2 - p_ϕ^2) * cos(ϕ)^2)
    end
    
    r_dot = p_r
    ϕ_dot = p_ϕ / r^2
    θ_dot = p_θ / (r^2 * cos(ϕ)^2)

    return [r, θ, ϕ, r_dot, θ_dot, ϕ_dot]
end

function cartesian_to_orbital_koopman(cartesian; length_unit = 1.0, μ = 1.0, extra=false)
    spherical = cartesian_to_spherical(cartesian)
    return spherical_to_orbital_koopman(spherical; length_unit, μ, extra)
end

function orbital_koopman_to_cartesian(orbital_koopman; length_unit = 1.0, μ = 1.0, prograde = true)
    spherical = orbital_koopman_to_spherical(orbital_koopman; length_unit, μ, prograde)
    return spherical_to_cartesian(spherical)
end

function cartesian_to_modified_equinoctial(cartesian; μ = 1.0)
    r_dot_v = dot(cartesian[1:3], cartesian[4:6])
    
    r_hat = normalize(cartesian[1:3])
    r_mag = norm(cartesian[1:3])
    
    h_vec = cross(cartesian[1:3], cartesian[4:6])
    h_mag = norm(h_vec)
    h_hat = normalize(h_vec)
    
    v_hat = (r_mag*cartesian[4:6] - r_dot_v*r_hat)/h_mag
    
    p = h_mag^2/μ
    k = h_hat[1]/(1.0 + h_hat[3])
    h = -h_hat[2]/(1.0 + h_hat[3])
    kk = k^2
    hh = h^2
    s2 = 1.0 + hh + kk
    tkh = 2.0 * k * h
    
    ecc = cross(cartesian[4:6], h_vec)/μ - r_hat
    
    f_hat = [1.0 - kk + hh, tkh, -2.0 * k]
    g_hat = [tkh, 1.0 + kk - hh, 2.0 * h]
    f_hat = f_hat / s2
    g_hat = g_hat / s2
    
    f = dot(ecc, f_hat)
    g = dot(ecc, g_hat)
    
    L = atan(r_hat[2] - v_hat[1], r_hat[1] + v_hat[2])
    
    return [p, f, g, h, k, L]
end

function modified_equinoctial_to_cartesian(modified_equinoctial; μ = 1.0)
    p, f, g, h, k, L = modified_equinoctial

    kk = k^2
    hh = h^2
    tkh = 2.0 * k * h
    s2 = 1.0 + hh + kk
    
    cL = cos(L)
    sL = sin(L)
    
    w = 1.0 + f * cL + g * sL
    r = p / w
    
    smp = sqrt(μ / p)
    
    f_hat = [1.0 - kk + hh, tkh, -2.0 * k] / s2
    g_hat = [tkh, 1.0 + kk - hh, 2.0 * h] / s2
    
    x = r * cL
    y = r * sL
    
    x_dot = -smp * (g + sL)
    y_dot = smp * (f + cL)

    return vcat(x*f_hat .+ y*g_hat, x_dot*f_hat .+ y_dot*g_hat)
end



