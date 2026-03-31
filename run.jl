# =============================================================================
# CFRP Oxidative Decomposition — Multi-Condition Inverse Modeling
# =============================================================================
# Homework 2: Introduction to Scientific Computing
#
# Overall workflow:
#   Task A — Construct the Catalyst.jl reaction system
#   Task B — Create synthetic TGA datasets (4 experimental cases)
#   Task C — Formulate the PEtab-based inverse problem
#   Task D — Perform parameter estimation, generate comparison, and validate results
#
# How to run:
#   1. Activate environment:   ] activate .
#   2. Install packages:      ] instantiate
#   3. Execute script:        include("run.jl")
# =============================================================================

using Statistics, CSV

# Load all required modules
include("src/model.jl")
include("src/data_generation.jl")
include("src/inverse_problem.jl")
include("src/visualization.jl")

function main()
    println("="^60)
    println("  CFRP Oxidative Decomposition — Inverse Modeling")
    println("="^60)

    # =========================================================================
    # TASK A: Forward Model Construction
    # =========================================================================
    println("\n📌 Task A: Building Catalyst.jl reaction network...")
    rn = build_cfrp_reaction_system()
    println("  ✓ Reaction system built: $(length(reactions(rn))) reactions, " *
            "$(length(species(rn))) species, $(length(parameters(rn))) parameters")

    # =========================================================================
    # TASK B: Synthetic Data Generation
    # =========================================================================
    println("\n📌 Task B: Generating synthetic TGA data...")
    measurements_df, clean_data = generate_all_synthetic_data(rn; noise_sigma=0.005, n_points=150, seed=42)
    println("  ✓ Generated $(nrow(measurements_df)) measurement points across 4 experiments")

    # Export measurement dataset
    CSV.write("results/measurements.csv", measurements_df)
    println("  ✓ Measurements saved to results/measurements.csv")

    # Visualize all experiments
    println("\n  Plotting all experiments...")
    plot_all_experiments(clean_data)

    # =========================================================================
    # TASK C: Inverse Problem Formulation
    # =========================================================================
    println("\n📌 Task C: Setting up PEtab inverse problem...")
    petab_prob = setup_petab_model(rn, measurements_df)
    println("  ✓ PEtab problem created with $(length(petab_prob.xnames)) parameters to estimate")
    println("  Parameters: ", join(string.(petab_prob.xnames), ", "))

    # =========================================================================
    # TASK D: Parameter Estimation and Validation
    # =========================================================================
    println("\n📌 Task D: Running parameter estimation...")
    result = run_parameter_estimation(petab_prob; n_multistarts=10)

    # --- Parameter Comparison ---
    comparison_df = compare_parameters(result, petab_prob)
    print_comparison_table(comparison_df)

    # Save comparison results
    CSV.write("results/parameter_comparison.csv", comparison_df)
    println("  ✓ Comparison table saved to results/parameter_comparison.csv")

    # --- Validation Visualization ---
    println("\n  Generating validation plot...")
    plot_validation(rn, result, petab_prob, clean_data)

    # =========================================================================
    # Discussion: Importance of Experiment 4
    # =========================================================================
    println("\n" * "="^70)
    println("  DISCUSSION")
    println("="^70)
    println("""
  Why Experiment 4 (5% O₂) is essential for identifying m₂ and m₃:

  The oxygen reaction orders m₂ and m₃ appear in the oxidation rate laws
  as r ∝ PO₂^m. If all experiments are conducted at the same oxygen level
  (21%), then PO₂^m becomes a constant multiplier that cannot be separated
  from the pre-exponential factor A. As a result, different combinations of
  A and m can produce identical model behavior, leading to parameter
  identifiability issues.

  Introducing Experiment 4 with a lower oxygen concentration (5%) changes
  this situation. The ratio of reaction rates between experiments at
  different oxygen levels becomes dependent on m:
      r(21%)/r(5%) = (0.21/0.05)^m
  This additional variation allows the model to distinguish between A and m,
  ensuring that both oxygen reaction orders can be uniquely estimated.
    """)

    println("="^60)
    println("  ✓ All tasks completed successfully!")
    println("="^60)

    return rn, measurements_df, clean_data, petab_prob, result, comparison_df
end

# Execute the complete workflow
rn, measurements_df, clean_data, petab_prob, result, comparison_df = main()