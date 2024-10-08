function zero_revolution_lambert(r1, r2, tof)
    return multi_revolution_lambert(r1, r2, tof, 0)
end

function multi_revolution_lambert(r1, r2, tof, n_revolutions; is_retrograde=false)
    v1 = copy(r1[1:3, :])
    v2 = copy(r2[1:3, :])

    for i in 1:size(r1, 2)
        tof_i = if length(tof) > 1
            tof[i]
        else
            tof
        end

        n_rev_i = if length(n_revolutions) > 1
            n_revolutions[i]
        else
            n_revolutions
        end

        temp = lambert_problem(r1[1:3, i], r2[1:3, i], tof_i[1], μ, n_rev_i[1]; is_retrograde)

        v1[:, i] = temp[1][:]
        v2[:, i] = temp[2][:]
    end

    return v1, v2
end

# FROM https://github.com/naoyaozaki/LambertProblem.jl/blob/main/src/LambertProblem.jl

function lambert_problem(r1, r2, tof, μ, multi_revs; is_retrograde=false)
    # 0. Sanity Check
    if tof <= 0
        @error "ERROR: Time of flight must be positive!"
    end
    if μ <= 0
        @error "ERROR: Gravity parameter must be positive!"
    end


    # 1. Preparation
    # Variables used in Lambert's Problem
    r1_norm = norm(r1)
    r2_norm = norm(r2)
    c = norm(r2 - r1)
    s = 0.5 * (r1_norm + r2_norm + c)
    λ2 = 1.0 - c / s
    tof_nomdim = sqrt(2.0 * μ / (s^3)) * tof

    # Basis vectors
    ivec_r1 = r1 / r1_norm
    ivec_r2 = r2 / r2_norm
    ivec_h = cross(ivec_r1, ivec_r2)
    ivec_h = ivec_h / norm(ivec_h) # Angular momentum vector
    # NOTE: ivec_h cannot be defined if r1 // r2

    if ivec_h[3] == 0.0
        @error "ERROR: The angular momentum vector has no z component, impossible to define automatically clock or counterclockwise!"
    elseif ivec_h[3] < 0.0
        λ = -sqrt(λ2)
        ivec_t1 = -cross(ivec_h, ivec_r1)
        ivec_t2 = -cross(ivec_h, ivec_r2)
    else
        λ = sqrt(λ2)
        ivec_t1 = cross(ivec_h, ivec_r1)
        ivec_t2 = cross(ivec_h, ivec_r2)
    end

    if is_retrograde # Retrograde Motion
        λ = -λ
        ivec_t1 = -ivec_t1
        ivec_t2 = -ivec_t2
    end


    # 2. Calculate x
    iter_all, x_all = find_xy(λ, tof_nomdim, multi_revs)


    # 3. Calculate Terminal Velocities
    vvec_1_all = zeros(length(x_all), 3)
    vvec_2_all = zeros(length(x_all), 3)

    γ = sqrt(0.5 * μ * s)
    ρ = (r1_norm - r2_norm) / c
    σ = sqrt(1.0 - ρ^2)

    for (i, x) in enumerate(x_all)
        # Calculate Velocity Norm
        y = sqrt(1.0 - λ^2 * (1.0 - x^2))
        v_r1 = γ * ((λ * y - x) - ρ * (λ * y + x)) / r1_norm
        v_r2 = -γ * ((λ * y - x) + ρ * (λ * y + x)) / r2_norm
        v_t1 = γ * σ * (y + λ * x) / r1_norm
        v_t2 = γ * σ * (y + λ * x) / r2_norm

        # Calculate Velocity Vectors
        vvec_1_all[i, :] = (v_r1 * ivec_r1 + v_t1 * ivec_t1)
        vvec_2_all[i, :] = (v_r2 * ivec_r2 + v_t2 * ivec_t2)
    end

    # Output
    return vvec_1_all, vvec_2_all, length(x_all)

end


""" A Method that Calculates tof from x
"""
function x2tof(x, λ, m_max, Δx_battin=0.01, Δx_lagrange=0.2)
    Δx = abs(x - 1.0)

    if (Δx > Δx_battin) && (Δx < Δx_lagrange)
        # Use Lagrange TOF Expression
        a = 1.0 / (1.0 - x^2) # Non dimensional semi-major axis
        if a > 0 # Ellipse
            α = 2.0 * acos(x)
            β = 2.0 * asin(λ / sqrt(a))
            tof = 0.5 * a^(3 / 2) * ((α - sin(α)) - (β - sin(β)) + 2.0 * m_max * π) # Non dimensinal TOF

        else # Hyperbola
            α = 2.0 * acosh(x)
            β = 2.0 * asinh(λ / sqrt(-a))
            tof = -0.5 * (-a)^(3 / 2) * ((α - sinh(α)) - (β - sinh(β))) # Non dimensinal TOF

        end

    elseif (Δx < Δx_battin)
        # Use Battin Series TOF Expression
        ρ = abs(x^2 - 1.0)
        z = sqrt(1.0 + λ^2 * (x^2 - 1.0))
        η = z - λ * x
        s1 = 0.5 * (1.0 - λ - x * η)
        q = 4.0 / 3.0 * hypergeometric_₂F₁(s1)
        tof = (η^3 * q + 4.0 * λ * η) / 2.0 + m_max * π / (ρ^1.5)

    else
        # Use Lancaster TOF Expression
        e = x^2 - 1.0
        z = sqrt(1.0 + λ^2 * e)
        y = sqrt(abs(e))
        if (e < 0.0)
            d = m_max * π + acos(x * z - λ * e)
        else
            d_temp = y * (z - λ * x) + (x * z - λ * e)
            if d_temp > 0.0
                d = log(d_temp)
            else
                @info "Fail to Calculate TOF using Lancaster TOF Expression."
                return NaN
            end
        end
        tof = (x - λ * z - d / y) / e

    end

    return tof
end

function hypergeometric_₂F₁(z, a=3.0, b=1.0, c=2.5, tol=1.0e-11)
    # Initilization
    sj, cj, sj1, cj1, err = 1.0, 1.0, 0.0, 0.0, 1.0

    # Iteration
    j = 0
    while (err > tol)
        cj1 = cj * (a + j) * (b + j) / (c + j) * z / (j + 1)
        sj1 = sj + cj1
        err = abs(cj1)
        sj, cj = sj1, cj1
        j += 1
        if j > 1000
            @error "ERROR: Hypergeometric Function Reaches Maximum Iteration."
            break
        end
    end

    return sj
end

function find_xy(λ, tof, m_multi_revs)
    # Requirements
    if abs(λ) >= 1
        @error "ERROR: Lambda must be more than 1."
    end
    if tof < 0
        @error "ERROR: Non dimensional tof must be a positive number."
    end

    # ----------------
    # 1. Detect m_max
    m_max = floor(tof / π)
    t_00 = acos(λ) + λ * sqrt(1.0 - λ^2) # Eq.(19) in Ref[1]
    t_0m = t_00 + m_max * π # Minimum Energy Transfer Time: Eq.(19) in Ref[1]
    t_1 = 2.0 / 3.0 * (1.0 - λ^3)

    if (m_multi_revs > 0) && (m_max > 0) && (tof < t_0m)
        x_tmin, t_min = find_tof_min_by_halley_method(0.0, t_0m, λ, m_max)

        if (tof < t_min) || isnan(t_min)
            m_max -= 1
        end
    end

    # Crop m_max to m_multi_revs
    m_max = Int8(min(m_multi_revs, m_max))

    # ----------------
    # 2. Calculate All Solutions in x,y
    x_all = zeros(2 * m_max + 1)
    iter_all = zeros(Int16, 2 * m_max + 1)

    # 2.1. Single Revolution Solution
    # Initial guess
    if tof >= t_00
        x_all[1] = -(tof - t_00) / (tof - t_00 + 4.0) #(t_00/tof)^(2/3) - 1.0
    elseif tof <= t_1
        x_all[1] = (5.0 * t_1 * (t_1 - tof)) / (2.0 * tof * (1 - λ^5)) + 1.0
    else
        x_all[1] = (tof / t_00)^(log(2.0) / log(t_1 / t_00)) - 1.0  #(t_00/tof)^(log2(t_1/t_00)) - 1.0
    end

    # Householder iterations
    iter_all[1], x_all[1] = find_x_by_householder(tof, x_all[1], λ, 0, 1.0e-5)

    # 2.2. Multi Revolution Solution
    for i = 1:m_max
        # Left Householder iterations
        temp = ((i * π + π) / (8.0 * tof))^(2.0 / 3.0)
        x_all[2*i] = (temp - 1.0) / (temp + 1.0)
        iter_all[2*i], x_all[2*i] = find_x_by_householder(tof, x_all[2*i], λ, i)

        # Right Householder iterations
        temp = ((8.0 * tof) / (i * π))^(2.0 / 3.0)
        x_all[2*i+1] = (temp - 1.0) / (temp + 1.0)
        iter_all[2*i+1], x_all[2*i+1] = find_x_by_householder(tof, x_all[2*i+1], λ, i)
    end

    return iter_all, x_all
end


"""
Find x that satisfies f(x)=tn(x)-tof=0 by Householder's method
"""
function find_x_by_householder(tof, xn, λ, m, tol_Δx=1.0e-8, max_iter=15)
    iter = 0
    while true
        tn = x2tof(xn, λ, m)

        # Cannot be calculated
        if isnan(tn)
            return iter, NaN
        end

        # Eqs.(22) in Ref[1]
        f(x, t) = t - tof
        df_dx(x, t) = (3.0 * t * x - 2.0 + 2.0 * λ^3 * x / sqrt(1.0 - λ^2 * (1.0 - x^2))) / (1.0 - x^2)
        d2f_dx2(x, t) = (3.0 * t + 5.0 * x * df_dx(x, t) + 2.0 * (1.0 - λ^2) * (λ^3) / (sqrt(1.0 - λ^2 * (1.0 - x^2))^3)) / (1.0 - x^2)
        d3f_dx3(x, t) = (7.0 * x * d2f_dx2(x, t) + 8.0 * df_dx(x, t) - 6.0 * (1.0 - λ^2) * (λ^5) * x / (sqrt(1.0 - λ^2 * (1.0 - x^2))^5)) / (1.0 - x^2)

        # Householder's Method
        xn_new = xn - f(xn, tn) * (
            (df_dx(xn, tn)^2 - 0.5 * f(xn, tn) * d2f_dx2(xn, tn))
            /
            (df_dx(xn, tn) * ((df_dx(xn, tn)^2) - f(xn, tn) * d2f_dx2(xn, tn)) + d3f_dx3(xn, tn) * (f(xn, tn)^2) / 6.0)
        )

        # Break condition
        if abs(xn_new - xn) < tol_Δx
            tn = x2tof(xn_new, λ, m)
            return iter, xn_new
        elseif iter > max_iter
            @info "Householder iteration reaches maximum iteration!"
            tn = x2tof(xn_new, λ, m)
            return iter, xn_new
        end

        # Update the value
        xn = xn_new
        iter += 1

    end

end


"""
Find minimum value of transfer time by Halley's method
"""
function find_tof_min_by_halley_method(xn, tn, λ, m_max, tol_Δx=1.0e-13, max_iter=12)
    iter = 0
    while true
        # Eqs.(22) in Ref[1]
        dt_dx(x, t) = (3.0 * t * x - 2.0 + 2.0 * λ^3 * x / sqrt(1.0 - λ^2 * (1.0 - x^2))) / (1.0 - x^2)
        d2t_dx2(x, t) = (3.0 * t + 5.0 * x * dt_dx(x, t) + 2.0 * (1.0 - λ^2) * (λ^3) / (sqrt(1.0 - λ^2 * (1.0 - x^2))^3)) / (1.0 - x^2)
        d3t_dx3(x, t) = (7.0 * x * d2t_dx2(x, t) + 8.0 * dt_dx(x, t) - 6.0 * (1.0 - λ^2) * (λ^5) * x / (sqrt(1.0 - λ^2 * (1.0 - x^2))^5)) / (1.0 - x^2)

        # Halley's Method
        xn_new = xn - (2.0 * dt_dx(xn, tn) * d2t_dx2(xn, tn)) / (2.0 * (d2t_dx2(xn, tn)^2) - dt_dx(xn, tn) * d3t_dx3(xn, tn))

        # Break condition
        if abs(xn_new - xn) < tol_Δx
            tn = x2tof(xn_new, λ, m_max)
            return xn_new, tn
        elseif iter > max_iter
            @info "Halley iteration reaches maximum iteration!"
            tn = x2tof(xn_new, λ, m_max)
            return xn_new, tn
        end

        # Update the value
        tn = x2tof(xn_new, λ, m_max)
        xn = xn_new
        iter += 1

        # Cannot be calculated
        if isnan(tn)
            return xn_new, tn
        end
    end
end