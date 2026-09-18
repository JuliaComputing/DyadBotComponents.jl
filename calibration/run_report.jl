# Compare 2 robot runs with the same measurements, so a control change shows up
# as a number and not as an impression.
#
# Run:  julia +1.12 --project=. run_report.jl data/run_clean.arrow data/run_permotor.arrow
#
# It reports, for each run:
#
#   step response   overshoot and settled error, SEPARATED BY DIRECTION, because
#                   the drives break away at different duties forward and in
#                   reverse and a one-sided overshoot is the signature
#   yaw             the wheel count difference, which the yaw loop holds at 0.
#                   A smaller spread means the drives themselves are matched.
#   9 Hz bursts     the fraction of time the tilt rate rings between 7 and 12 Hz.
#                   The 2026-09-18 limit cycle sat there.
#   timing          late ticks and the compute time.

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "robot_log_to_releases.jl"))
using Arrow, DataFrames, Statistics, Printf, DSP, Plots

const TS = CFG["robot_log"]["Ts"]

"Steps of the reference, as (time, from, to) with the settled levels."
function steps(t, ref)
    out = NTuple{3, Float64}[]
    # The reference is a filtered square, so it moves slowly and has no jump to
    # find. Take the level changes from the sign of the plateau instead. `to` is
    # then the settled level, +-amplitude, not the reference at the transition
    # sample, which the 2 first-order filters hold near zero.
    lvl = sign.(round.(ref ./ max(maximum(abs, ref), 1e-9), digits = 1))
    for k in 2:length(lvl)
        if lvl[k] != lvl[k - 1] && lvl[k] != 0
            push!(out, (t[k], ref[k - 1], lvl[k] * maximum(abs, ref)))
        end
    end
    return out
end

"""
    envelope(x, lo, hi, fs)

The amplitude of `x` between `lo` and `hi` Hz, as a slow envelope. A band-pass
followed by the absolute value and a 0.2 s smoother.
"""
function envelope(x, lo, hi, fs)
    bp = digitalfilter(Bandpass(lo / (fs / 2), hi / (fs / 2)), Butterworth(4))
    y = abs.(filtfilt(bp, x .- mean(x)))
    n = round(Int, 0.2 * fs)
    return [mean(@view y[max(1, k - n):k]) for k in eachindex(y)]
end

function report(path)
    raw = read_log(path)
    conv = convert_log(raw)
    t = conv.timestamp
    x = conv[!, "plant.x"]
    ref = conv.pos_ref
    println("\n=== ", basename(path), " ===")
    @printf("  %.1f s, %d ticks\n", t[end], length(t))

    # --- timing ---
    dt = diff(Float64.(raw.now)) ./ 1e9
    late = count(>(1.5TS), dt)
    @printf("  late ticks (>%.0f ms): %d   worst gap %.0f ms   compute mean %.2f ms, max %.2f ms\n",
            1500TS, late, 1000maximum(dt), mean(raw.compute_ns) / 1e6, maximum(raw.compute_ns) / 1e6)

    # --- steps, by direction ---
    # Overshoot is the peak PAST THE LOCAL SETTLED VALUE, not past the
    # reference. The robot holds a slowly drifting offset from the reference, so
    # a peak measured against the reference reports that drift and not the step
    # response. Skip the first 30 s: that is the release and the first settling.
    st = steps(t, ref)
    over = Dict("forward" => Float64[], "reverse" => Float64[])
    travel = Dict("forward" => Float64[], "reverse" => Float64[])
    drift = Float64[]
    for (ts, from, to) in st
        (ts < 30.0 || ts + 5.0 > t[end]) && continue
        dir = to > from ? "forward" : "reverse"
        w = (t .>= ts) .& (t .<= ts + 4.0)
        s = (t .>= ts + 3.5) .& (t .<= ts + 5.0)
        (count(w) < 100 || count(s) < 50) && continue
        x0 = mean(x[(t .>= ts - 1.0) .& (t .< ts)])
        xf = mean(x[s])
        pk = to > from ? maximum(x[w]) : minimum(x[w])
        push!(over[dir], abs(pk - xf))
        push!(travel[dir], abs(xf - x0))
        push!(drift, xf - to)   # `to` is the settled level from steps(), not ref at the transition
    end
    for dir in ("forward", "reverse")
        isempty(over[dir]) && continue
        @printf("  %-7s steps %2d   overshoot %4.2f +- %4.2f cm   travel %4.1f cm   (%3.0f %% of travel)\n",
                dir, length(over[dir]), 100mean(over[dir]), 100std(over[dir]), 100mean(travel[dir]),
                100mean(over[dir]) / 100mean(travel[dir]) * 100)
    end
    if !isempty(over["forward"]) && !isempty(over["reverse"])
        @printf("  overshoot asymmetry (forward - reverse): %+.2f cm\n",
                100 * (mean(over["forward"]) - mean(over["reverse"])))
    end
    isempty(drift) || @printf("  settled offset from the reference: mean %+.1f cm, spread %.1f cm, worst %.1f cm\n",
                              100mean(drift), 100std(drift), 100maximum(abs, drift))

    # --- settled tracking ---
    err = x .- ref
    hold = trues(length(t))
    for (ts, _, _) in st
        hold[(t .>= ts) .& (t .<= ts + 3.0)] .= false
    end
    @printf("  settled tracking RMS %.1f cm over %.0f s\n", 100sqrt(mean(abs2, err[hold])), sum(hold) * TS)

    # --- yaw ---
    d = Float64.(raw.wheel_a) .- Float64.(raw.wheel_b)
    d .-= d[1]
    @printf("  wheel difference: range %.0f counts (%.0f to %.0f), final %.0f, yaw duty max %.3f\n",
            maximum(d) - minimum(d), minimum(d), maximum(d), d[end], maximum(abs, raw.yaw_duty))

    # --- 9 Hz ring ---
    e = envelope(Float64.(raw.gyro_y), 7.0, 12.0, 1 / TS)
    for thr in (0.05, 0.1)
        @printf("  ring 7-12 Hz above %.2f rad/s: %.1f %% of the run\n", thr, 100count(>(thr), e) / length(e))
    end
    @printf("  ring envelope: median %.3f, p99 %.3f, max %.3f rad/s\n",
            median(e), quantile(e, 0.99), maximum(e))
    return (; t, x, ref, e, d = d, name = basename(path))
end

runs = [report(p) for p in (isempty(ARGS) ? ["data/run_clean.arrow", "data/run_permotor.arrow"] : ARGS)]

C = ["#8a8985", "#eb6834", "#2a78d6"]
p1 = plot(; ylabel = "position [cm]", title = "Run comparison")
plot!(p1, runs[end].t, 100 .* runs[end].ref; label = "reference", color = "#1baf7a", ls = :dash)
p2 = plot(; ylabel = "wheel difference [counts]")
p3 = plot(; ylabel = "7-12 Hz ring [rad/s]", xlabel = "time [s]")
for (i, r) in enumerate(runs)
    plot!(p1, r.t, 100 .* r.x; label = r.name, color = C[i])
    plot!(p2, r.t, r.d; label = r.name, color = C[i])
    plot!(p3, r.t, r.e; label = r.name, color = C[i])
end
p = plot(p1, p2, p3; layout = (3, 1), size = (1200, 850), lw = 1.3, framestyle = :box, gridalpha = 0.25,
         foreground_color_legend = nothing, background_color_legend = nothing, dpi = 130, left_margin = 8Plots.mm)
savefig(p, joinpath(PLOTS, "run_comparison.png"))
println("\nplot written to ", joinpath(PLOTS, "run_comparison.png"))
