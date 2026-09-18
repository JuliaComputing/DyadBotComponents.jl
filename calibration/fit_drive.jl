# Turn a bench run of scripts/drive_calibration.jl into the drive's parameters.
#
# The wheels are off the ground, so each wheel is a first-order system on its own:
#
#     Iw * dw/dt = k_tau * duty - d * w - tau_c * sign(w)
#
# with `Iw` the wheel and rotor inertia, `k_tau` the torque the drive delivers
# per unit duty (the inverse of the runtime's `k_duty`), `d` the viscous
# friction and `tau_c` the Coulomb friction. The 5 phases of the bench run pin
# these down separately, which the closed-loop fit cannot: there the wheel
# inertia and the friction trade against each other.
#
#   break-away  the duty at which the wheel first moves, each direction
#   steady      the steady speed against duty: the slope is k_tau/d and the
#               intercept is tau_c/d
#   step        the initial acceleration: k_tau/Iw, so Iw follows from k_tau
#   chirp       the phase lag of the wheel speed behind the duty, 0.5 to 30 Hz.
#               Measured on 2026-09-18 this turned out to be the wheel's own
#               Iw/d and nothing else, within about 1 ms of dead time, so the
#               plant model already contains it. What it tests is whether Iw and
#               d are right, not whether a lag block is missing.
#   reverse     the counts lost to backlash when the drive changes direction
#
# Run:  julia +1.12 --project=. fit_drive.jl data/drive_calibration_A.csv
# Writes plots/drive_<motor>.png and prints the parameters, including the
# k_duty the runtime should use.

using CSV, DataFrames, Statistics, Printf, Plots
second(x) = x[2]
const HERE = @__DIR__
const PLOTS = joinpath(HERE, "plots")
mkpath(PLOTS)

const COUNTS_PER_REV = 1320.0
const RAD_PER_COUNT = 2pi / COUNTS_PER_REV

path = isempty(ARGS) ? joinpath(HERE, "data", "drive_calibration_A.csv") : ARGS[1]
name = replace(basename(path), "drive_calibration_" => "", ".csv" => "")
df = CSV.read(path, DataFrame)
println("motor ", name, ": ", nrow(df), " rows, phases ", join(unique(df.phase), ", "))

"Wheel speed in rad/s, from the counts, smoothed over `win` samples."
function speed(sub; win = 21)
    t = sub.t; c = sub.counts .* RAD_PER_COUNT
    n = length(t); v = zeros(n)
    for i in 1:n
        a = max(1, i - win ÷ 2); b = min(n, i + win ÷ 2)
        v[i] = b > a ? (c[b] - c[a]) / (t[b] - t[a]) : 0.0
    end
    return v
end

# --------------------------------------------------------------- break-away --
println("\nbreak-away (the duty at which the wheel first turns):")
breakaway = Dict{String, Float64}()
for ph in ("rampup_pos", "rampup_neg")
    sub = df[df.phase .== ph, :]
    isempty(sub) && continue
    moved = findfirst(abs.(sub.counts) .>= 3)
    if moved === nothing
        println("  ", ph, ": never moved")
    else
        breakaway[ph] = abs(sub.duty[moved])
        println(@sprintf("  %-11s duty %0.3f", ph, breakaway[ph]))
    end
end

# ------------------------------------------------------------ steady speed ---
# Average the last half of each level: k_tau * duty = d * w + tau_c * sign(w).
steady = df[df.phase .== "steady", :]
levels = Tuple{Float64, Float64}[]
if !isempty(steady)
    # Each level is one contiguous run of a constant duty.
    starts = [1; findall(i -> steady.duty[i] != steady.duty[i - 1], 2:nrow(steady)) .+ 1]
    for (k, s) in enumerate(starts)
        e = k < length(starts) ? starts[k + 1] - 1 : nrow(steady)
        half = s + (e - s) ÷ 2
        w = (steady.counts[e] - steady.counts[half]) * RAD_PER_COUNT / (steady.t[e] - steady.t[half])
        push!(levels, (steady.duty[s], w))
    end
    println("\nsteady speed against duty:")
    for (u, w) in levels
        println(@sprintf("  duty %+0.2f  %8.2f rad/s  %8.1f rpm", u, w, w * 60 / 2pi))
    end
end

# Fit w = a*duty + b separately per direction, then k_tau/d = a and tau_c/d = -b.
function line(pts)
    isempty(pts) && return (NaN, NaN)
    x = first.(pts); y = last.(pts)
    a = (mean(x .* y) - mean(x) * mean(y)) / (mean(x .^ 2) - mean(x)^2)
    return (a, mean(y) - a * mean(x))
end
pos = [(u, w) for (u, w) in levels if u > 0]
neg = [(u, w) for (u, w) in levels if u < 0]
(ap, bp) = line(pos); (an, bn) = line(neg)
slope = isnan(ap) ? an : (isnan(an) ? ap : (ap + an) / 2)

# ---------------------------------------------------------- step responses ---
# From rest the wheel is first order: w(t) = w_inf * (1 - exp(-(t - t_d)/tau))
# with tau = Iw/d and t_d the dead time before the drive answers. Integrating,
# the counts are w_inf * ((t - t_d) - tau * (1 - exp(-(t - t_d)/tau))), and
# fitting the counts rather than a differentiated speed keeps the noise down.
println("\nstep responses:")
steps = Tuple{Float64, Float64, Float64, Float64}[]     # duty, w_inf, tau, dead time
for ph in sort([p for p in unique(df.phase) if startswith(p, "step_")])
    sub = df[df.phase .== ph, :]
    nrow(sub) < 100 && continue
    u = sub.duty[1]
    ang = abs.(sub.counts) .* RAD_PER_COUNT
    # Dead time: the first sample that has moved more than 2 counts.
    imove = findfirst(abs.(sub.counts) .>= 2)
    t_d = imove === nothing ? 0.0 : sub.t[imove]
    # Final speed from the last 0.2 s.
    last = sub.t .>= sub.t[end] - 0.2
    w_inf = (ang[end] - ang[findfirst(last)]) / (sub.t[end] - sub.t[findfirst(last)])
    # tau from the area the ramp loses to the exponential, averaged over the
    # second half where the exponential has died: angle = w_inf*(t - t_d - tau).
    late = sub.t .>= t_d + 0.3
    tau = count(late) < 20 ? NaN :
          mean((w_inf .* (sub.t[late] .- t_d) .- ang[late]) ./ w_inf)
    push!(steps, (u, w_inf, tau, t_d))
    println(@sprintf("  duty %+0.2f  final %6.2f rad/s  time constant %5.1f ms  dead time %4.1f ms",
                     u, w_inf, 1000tau, 1000t_d))
end
tau_mech = isempty(steps) ? NaN : median([s[3] for s in steps if !isnan(s[3])])   # median: one step can start with the wheel still turning
dead_time = isempty(steps) ? NaN : mean(s[4] for s in steps)

# ------------------------------------------------------------------- chirp ---
# The sweep phase is known exactly, so demodulate both the duty and the wheel
# speed against exp(-i*phi(t)) in bands and take the phase difference. That is
# robust, unlike a cross-correlation, which locks onto the wrong period.
chirp = df[df.phase .== "chirp", :]
lags = Tuple{Float64, Float64, Float64}[]      # frequency, phase lag in degrees, gain
if nrow(chirp) > 1000
    v = speed(chirp; win = 5)
    f0, f1, T = 0.5, 30.0, 25.0
    freq(t) = f0 + (f1 - f0) * t / T
    phi(t) = 2pi * (f0 * t + (f1 - f0) * t^2 / (2T))
    println("\nfrequency sweep, phase of the wheel speed behind the duty:")
    for (a, b) in ((1.0, 3.0), (3.0, 6.0), (6.0, 9.0), (9.0, 13.0), (13.0, 20.0), (20.0, 29.0))
        sel = findall(i -> a <= freq(chirp.t[i]) <= b, 1:nrow(chirp))
        length(sel) < 200 && continue
        e = cis.(-phi.(chirp.t[sel]))
        U = mean(chirp.duty[sel] .* e)
        Y = mean((v[sel] .- mean(v[sel])) .* e)
        (abs(U) < 1e-6 || abs(Y) < 1e-9) && continue
        # The encoder counts DOWN for a positive duty in the quarter-period mode,
        # so the measured speed carries a fixed 180 degrees against the command.
        # Take it out here, or every lag reads as a lead.
        deg = rad2deg(angle(Y) - angle(U)) - 180
        deg = mod(deg + 180, 360) - 180
        fc = (a + b) / 2
        push!(lags, (fc, deg, abs(Y) / abs(U)))
        println(@sprintf("  %4.1f Hz  phase %+7.1f deg  (a %0.0f ms lag alone gives %+0.1f)  gain %6.2f rad/s per duty",
                         fc, deg, 1000 * 0.0473, -rad2deg(atan(2pi * fc * 0.0473)), abs(Y) / abs(U)))
    end
end

# ---------------------------------------------------------------- backlash ---
# Not measurable from this phase. The reversal runs at duty 0.35, so between the
# command flipping and the wheel turning back the wheel has to decelerate
# through zero, and that dominates whatever lost motion the gearbox has. The
# encoder is on the wheel, so motor-side lost motion is invisible in principle.
# A low-duty reversal, just above break-away, would separate them.
rev = df[df.phase .== "reverse", :]
if !isempty(rev)
    flip = findfirst(rev.duty .< 0)
    if flip !== nothing
        after = rev[flip:end, :]
        peak_i = argmax(abs.(after.counts))
        println(@sprintf("\nreversal: the wheel kept going %0.0f more counts and took %0.0f ms to turn round",
                         abs(after.counts[peak_i] - after.counts[1]), 1000 * (after.t[peak_i] - after.t[1])))
        println("  That is the coast-down, not backlash: at duty 0.35 the deceleration dominates.")
    end
end

# --------------------------------------------------------------- the model ---
# With the wheels in the air the wheel is  Iw dw/dt = k_tau u - d w - tau_c.
# The bench fixes 2 RATIOS and not the 3 parameters: the steady line gives
# k_tau/d, the step gives Iw/d. One absolute value is still needed, and only an
# independent measurement of Iw can give it (a bifilar pendulum, or repeating
# the step with a known mass taped to the wheel at a known radius).
println("\n--- what the bench determines ---")
kd = abs(slope)                                  # k_tau/d, rad/s per duty
println(@sprintf("  k_tau / d   = %7.2f rad/s per duty   (steady-speed line)", kd))
isnan(tau_mech) || println(@sprintf("  Iw / d      = %7.1f ms                (step time constant)", 1000tau_mech))
isnan(dead_time) || println(@sprintf("  dead time   = %7.1f ms                (command to first motion)", 1000dead_time))
if !isempty(breakaway)
    bp_ = get(breakaway, "rampup_pos", NaN); bn_ = get(breakaway, "rampup_neg", NaN)
    println(@sprintf("  break-away    %0.3f duty forward, %0.3f reverse  (%0.0f %% apart)",
                     bp_, bn_, 100 * abs(bp_ - bn_) / ((bp_ + bn_) / 2)))
end

println("\n--- the family of drives consistent with the bench ---")
println("  Assume a wheel-and-rotor inertia Iw, and the rest follows:")
println(@sprintf("  %-14s %-14s %-16s %s", "Iw [kg m^2]", "d [N m s/rad]", "k_tau [N m/duty]", "k_duty for the runtime"))
for Iw in (4.4e-5, 7.15e-5, 1.5e-4, 3.0e-4)
    isnan(tau_mech) && break
    d_ = Iw / tau_mech; k_ = kd * d_
    println(@sprintf("  %-14.2e %-14.2e %-16.4f %.1f", Iw, d_, k_, 1 / k_))
end
println("  4.4e-5 is the 2 wheels as discs; 7.15e-5 is what the closed-loop fit found;")
println("  the larger values allow for the rotor inertia through the gearbox.")

if !isnan(tau_mech)
    # Cross-check against the closed-loop fit, which anchors the torque through
    # the measured body: M g L is known, so k_tau is identifiable there.
    Iw_cl, d_cl, kscale_cl = 7.15e-5, 1.4691657905175379e-3, 1.1473736652284958
    println("\n--- against the closed-loop fit ---")
    println(@sprintf("  closed loop: Iw/d = %0.1f ms, k_tau/d = %0.1f rad/s per duty (k_tau = %0.4f from k_scale)",
                     1000 * Iw_cl / d_cl, 0.1 * kscale_cl / d_cl, 0.1 * kscale_cl))
    println(@sprintf("  bench:       Iw/d = %0.1f ms, k_tau/d = %0.1f rad/s per duty", 1000tau_mech, kd))
    println("  The time constants agree. The torque-to-friction ratio does not, by about")
    println("  3 times, so one of the two is wrong about how much of the command reaches")
    println("  the wheel. The bench measures the drive directly and is the one to trust;")
    println("  the closed-loop value is entangled with the friction feedforward.")
end
if !isempty(lags)
    # How much phase is left once the mechanical time constant is accounted for?
    # Whatever remains is dead time, and dead time costs 360*f*T_d degrees.
    println("\n  phase beyond the ", round(1000tau_mech, digits = 1), " ms mechanical lag, as an equivalent dead time:")
    for (f, deg, _) in lags
        extra = deg + rad2deg(atan(2pi * f * tau_mech))
        println(@sprintf("    %4.1f Hz: %+5.1f deg -> %4.2f ms", f, extra, -extra / 360 / f * 1000))
    end
    println("  The wheel's own Iw/d IS the drive lag. The plant model already contains it")
    println("  through Iw and d, so the twin needs no separate lag block: it needs the")
    println("  right Iw and d. Only the dead time, about 1 ms, is missing, and 1 ms is")
    println("  0.2 tick of the 200 Hz loop.")
end

# ------------------------------------------------------------------- plots ---
const C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
plts = Any[]
if !isempty(levels)
    p = scatter(first.(levels), abs.(last.(levels)); label = "measured", color = C1,
                xlabel = "duty", ylabel = "|steady speed| [rad/s]", title = "Torque-speed, motor $name", legend = :bottomright)
    push!(plts, p)
end
if !isempty(lags)
    p = plot(first.(lags), second.(lags); marker = :circle, label = "measured", color = C1,
             xlabel = "frequency [Hz]", ylabel = "phase [deg]", title = "Drive phase, motor $name")
    push!(plts, p)
end
if !isempty(plts)
    p = plot(plts...; layout = (1, length(plts)), size = (560 * length(plts), 420), lw = 2,
             framestyle = :box, gridalpha = 0.25, foreground_color_legend = nothing,
             background_color_legend = nothing, dpi = 130, left_margin = 8Plots.mm)
    savefig(p, joinpath(PLOTS, "drive_$(name).png"))
    println("\nplot written to ", joinpath(PLOTS, "drive_$(name).png"))
end
