# =============================================================================
# Visualization: Parameter Estimation Validation Plots
# =============================================================================

using Plots, Printf

"""
    plot_validation(rn, result, petab_prob, clean_data)

Generate the validation figure: Mass Fraction vs Temperature for Exp 2 (21% O₂)
and Exp 4 (5% O₂), including both noisy synthetic data and the fitted model response.

Stores the figure in results/validation_plot.png
"""
function plot_validation(rn, result, petab_prob, clean_data)
    # Recover optimized parameter values
    recovered = Dict(zip(petab_prob.xnames, result.xmin))
    opt_params = Dict{Symbol, Float64}()
    for (k, v) in recovered
        opt_params[k] = v
    end
    # Include fixed parameter
    opt_params[:nu_char] = 0.2

    # Run simulations with optimized parameters for Exp2 and Exp4
    p = plot(
        xlabel = "Temperature [K]",
        ylabel = "Mass Fraction",
        title  = "TGA Validation: Exp 2 (21% O₂) vs Exp 4 (5% O₂)",
        legend = :bottomleft,
        size   = (800, 500),
        dpi    = 150,
        grid   = true,
        framestyle = :box,
    )

    # --- Exp 2: 21% O₂, β = 5 K/min ---
    exp2_data = clean_data["Exp2"]
    # Display noisy measurements
    scatter!(p, exp2_data.temperatures, exp2_data.mass_noisy,
        label  = "Exp 2 data (21% O₂)",
        color  = :blue,
        alpha  = 0.3,
        markersize = 2,
        markerstrokewidth = 0,
    )
    # Compute fitted response using optimized parameters
    times_opt2, temps_opt2, mass_opt2 = simulate_experiment(
        rn, "Exp2", EXPERIMENTS["Exp2"], opt_params; n_points=300
    )
    plot!(p, temps_opt2, mass_opt2,
        label     = "Exp 2 fit (21% O₂)",
        color     = :blue,
        linewidth = 2.5,
    )

    # --- Exp 4: 5% O₂, β = 5 K/min ---
    exp4_data = clean_data["Exp4"]
    # Display noisy measurements
    scatter!(p, exp4_data.temperatures, exp4_data.mass_noisy,
        label  = "Exp 4 data (5% O₂)",
        color  = :red,
        alpha  = 0.3,
        markersize = 2,
        markerstrokewidth = 0,
    )
    # Compute fitted response using optimized parameters
    times_opt4, temps_opt4, mass_opt4 = simulate_experiment(
        rn, "Exp4", EXPERIMENTS["Exp4"], opt_params; n_points=300
    )
    plot!(p, temps_opt4, mass_opt4,
        label     = "Exp 4 fit (5% O₂)",
        color     = :red,
        linewidth = 2.5,
    )

    # Export figure
    savefig(p, "results/validation_plot.png")
    println("\n  Validation plot saved to: results/validation_plot.png")

    return p
end

"""
    plot_all_experiments(clean_data)

Generate mass-loss curves for all four experiments (clean + noisy) for comparison.
Stores the figure in results/all_experiments.png
"""
function plot_all_experiments(clean_data)
    colors = Dict(
        "Exp1" => :blue,
        "Exp2" => :green,
        "Exp3" => :orange,
        "Exp4" => :red,
    )
    labels = Dict(
        "Exp1" => "Exp 1: 2.5 K/min, 21% O₂",
        "Exp2" => "Exp 2: 5.0 K/min, 21% O₂",
        "Exp3" => "Exp 3: 10.0 K/min, 21% O₂",
        "Exp4" => "Exp 4: 5.0 K/min, 5% O₂",
    )

    p = plot(
        xlabel = "Temperature [K]",
        ylabel = "Mass Fraction",
        title  = "Synthetic TGA Data — All Experiments",
        legend = :bottomleft,
        size   = (800, 500),
        dpi    = 150,
        grid   = true,
        framestyle = :box,
    )

    for exp_id in ["Exp1", "Exp2", "Exp3", "Exp4"]
        data = clean_data[exp_id]
        # Show noisy measurements as scatter points
        scatter!(p, data.temperatures, data.mass_noisy,
            label  = labels[exp_id] * " (noisy)",
            color  = colors[exp_id],
            alpha  = 0.25,
            markersize = 1.5,
            markerstrokewidth = 0,
        )
        # Show clean simulated profile
        plot!(p, data.temperatures, data.mass_clean,
            label     = labels[exp_id] * " (true)",
            color     = colors[exp_id],
            linewidth = 2,
        )
    end

    savefig(p, "results/all_experiments.png")
    println("  All-experiments plot saved to: results/all_experiments.png")

    return p
end