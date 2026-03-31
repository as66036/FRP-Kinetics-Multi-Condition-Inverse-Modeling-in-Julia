# =============================================================================
# Task B: Synthetic Dataset Creation
# =============================================================================
#
# Produce simulated TGA (Thermogravimetric Analysis) data for all four experiments
# using the predefined true parameters. Introduce Gaussian noise to mimic
# experimental measurement uncertainty. The final output is formatted to be
# compatible with PEtab.
# =============================================================================

using OrdinaryDiffEq, DataFrames, Random

"""
    simulate_experiment(rn, exp_id, exp_cond, params; n_points=150)

Run a simulation for a single TGA experiment based on the reaction system.

Arguments:
- `rn`: Fully defined ReactionSystem
- `exp_id`: Identifier for the experiment (e.g., "Exp1")
- `exp_cond`: Named tuple containing (beta, PO2)
- `params`: Dictionary of true parameter values
- `n_points`: Number of sampled time points

Returns: (times, temperatures, mass_values)
"""
function simulate_experiment(rn, exp_id, exp_cond, params; n_points=150)
    # Determine simulation duration from temperature evolution
    # T = T0 + β*t  →  final time = (T_final - T0) / β
    t_final = (T_FINAL - INITIAL_TEMP) / exp_cond.beta

    # Assign initial values to all species
    @unpack M, C, F, Temp, G1, G2 = rn
    u0 = [
        M    => INITIAL_M,
        C    => INITIAL_C,
        F    => INITIAL_F,
        Temp => INITIAL_TEMP,
        G1   => INITIAL_G1,
        G2   => INITIAL_G2,
    ]

    # Assign model parameters (true kinetic values + experimental conditions)
    @unpack log10_A1, E1, n1, nu_char, log10_A2, E2, n2, m2,
            log10_A3, E3, n3, m3, PO2, beta = rn
    pmap = [
        log10_A1 => params[:log10_A1],
        E1       => params[:E1],
        n1       => params[:n1],
        nu_char  => params[:nu_char],
        log10_A2 => params[:log10_A2],
        E2       => params[:E2],
        n2       => params[:n2],
        m2       => params[:m2],
        log10_A3 => params[:log10_A3],
        E3       => params[:E3],
        n3       => params[:n3],
        m3       => params[:m3],
        PO2      => exp_cond.PO2,
        beta     => exp_cond.beta,
    ]

    # Define and solve the ODE system
    oprob = ODEProblem(rn, u0, (0.0, t_final), pmap)
    sol = solve(oprob, Rodas5P(); abstol=1e-10, reltol=1e-8, saveat=t_final/n_points)

    # Retrieve simulation outputs
    times = sol.t
    # Total solid mass = M + C + F (gas products excluded)
    mass_values = sol[M] .+ sol[C] .+ sol[F]
    temperatures = sol[Temp]

    return times, temperatures, mass_values
end

"""
    generate_all_synthetic_data(rn; noise_sigma=0.005, n_points=150, seed=42)

Create noisy synthetic TGA datasets for all experimental conditions.

Arguments:
- `rn`: Fully defined ReactionSystem
- `noise_sigma`: Standard deviation of Gaussian noise (default 0.005 = 0.5%)
- `n_points`: Number of sampled data points per experiment
- `seed`: Seed for reproducibility of random noise

Returns: (measurements_df, clean_data)
- `measurements_df`: DataFrame formatted for PEtab
- `clean_data`: Dictionary storing clean and noisy simulation results
"""
function generate_all_synthetic_data(rn; noise_sigma=0.005, n_points=150, seed=42)
    Random.seed!(seed)

    # Container for storing noise-free results (useful for plotting/validation)
    clean_data = Dict{String, NamedTuple}()

    # Arrays for assembling PEtab measurement dataset
    all_obs_ids    = String[]
    all_sim_ids    = String[]
    all_times      = Float64[]
    all_measurements = Float64[]

    for (exp_id, exp_cond) in sort(collect(EXPERIMENTS))
        println("  Simulating $exp_id: β = $(exp_cond.beta*60) K/min, PO₂ = $(exp_cond.PO2)")

        times, temperatures, mass_clean = simulate_experiment(
            rn, exp_id, exp_cond, TRUE_PARAMS; n_points=n_points
        )

        # Add Gaussian-distributed noise (σ = 0.005)
        noise = noise_sigma .* randn(length(mass_clean))
        mass_noisy = mass_clean .+ noise
        # Ensure values remain within physical limits [0, 1]
        mass_noisy = clamp.(mass_noisy, 0.0, 1.0)

        # Store both clean and noisy results for each experiment
        clean_data[exp_id] = (
            times = times,
            temperatures = temperatures,
            mass_clean = mass_clean,
            mass_noisy = mass_noisy,
        )

        # Populate PEtab measurement arrays
        for i in eachindex(times)
            push!(all_obs_ids, "mass_total")
            push!(all_sim_ids, exp_id)
            push!(all_times, times[i])
            push!(all_measurements, mass_noisy[i])
        end
    end

    # Create PEtab-compatible DataFrame
    measurements_df = DataFrame(
        obs_id        = all_obs_ids,
        simulation_id = all_sim_ids,
        time          = all_times,
        measurement   = all_measurements,
    )

    return measurements_df, clean_data
end