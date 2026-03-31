# =============================================================================
# Task A: Forward Simulation — CFRP Thermal Decomposition Reaction Model
# =============================================================================
#
# Reaction network describing CFRP degradation:
#   Reaction 1: Matrix decomposition without oxygen  — M → ν_char·C + (1-ν_char)·G1
#   Reaction 2: Oxidation of char (requires O₂)     — C + O₂ → G2
#   Reaction 3: Oxidation of fiber (requires O₂)    — F + O₂ → G2
#
# Temperature evolves dynamically: dT/dt = β (constant heating rate)
# Reaction rates follow Arrhenius kinetics with power-law dependence on species and oxygen
# =============================================================================

using Catalyst, ModelingToolkit

"""
    build_cfrp_reaction_system()

Constructs the Catalyst.jl ReactionSystem representing CFRP oxidative degradation.

Species included: M (Matrix), C (Char), F (Fiber), Temp (Temperature), G1 (Volatile gas), G2 (Oxidation gas)
Parameters include kinetic constants and experimental conditions.

Returns a fully defined ReactionSystem ready for numerical simulation.
"""
function build_cfrp_reaction_system()
    # Time variable
    t = default_t()

    # Declare chemical species
    @species M(t) C(t) F(t) Temp(t) G1(t) G2(t)

    # Define model parameters
    # Kinetic parameters (log10 form used for numerical robustness)
    @parameters log10_A1 E1 n1 nu_char
    @parameters log10_A2 E2 n2 m2
    @parameters log10_A3 E3 n3 m3
    # External experimental parameters
    @parameters PO2 beta

    # Universal gas constant [J/(mol·K)]
    R_gas = 8.314

    # --- Reaction rate definitions ---
    # max(species, 0) ensures stability during optimization (avoids invalid fractional powers)
    
    # Reaction 1: Matrix thermal decomposition (no oxygen involvement)
    # r1 = A1 * exp(-E1/(R*T)) * M^n1
    r1 = (10.0^log10_A1) * exp(-E1 / (R_gas * Temp)) * max(M, 0.0)^n1

    # Reaction 2: Char oxidation (oxygen-dependent)
    # r2 = A2 * exp(-E2/(R*T)) * C^n2 * PO2^m2
    r2 = (10.0^log10_A2) * exp(-E2 / (R_gas * Temp)) * max(C, 0.0)^n2 * PO2^m2

    # Reaction 3: Fiber oxidation (oxygen-dependent)
    # r3 = A3 * exp(-E3/(R*T)) * F^n3 * PO2^m3
    r3 = (10.0^log10_A3) * exp(-E3 / (R_gas * Temp)) * max(F, 0.0)^n3 * PO2^m3

    # --- Reaction definitions ---
    rxns = [
        # Temperature evolution: linear heating profile (dT/dt = beta)
        Reaction(beta, nothing, [Temp], nothing, [1]),

        # Reaction 1: Matrix → Char + Volatiles
        # only_use_rate=true ensures custom rate expression is used directly
        Reaction(r1, [M], [C, G1], [1], [nu_char, 1 - nu_char]; only_use_rate=true),

        # Reaction 2: Char → Oxidation gas
        Reaction(r2, [C], [G2], [1], [1]; only_use_rate=true),

        # Reaction 3: Fiber → Oxidation gas
        Reaction(r3, [F], [G2], [1], [1]; only_use_rate=true),
    ]

    # Assemble and finalize the reaction system
    @named cfrp = ReactionSystem(rxns, t)
    cfrp = complete(cfrp)

    return cfrp
end

# =============================================================================
# Reference (true) parameter values used for synthetic data generation
# =============================================================================
const TRUE_PARAMS = Dict(
    :log10_A1 => 5.0,       # Pre-exponential factor (log10) for matrix decomposition
    :E1       => 120_000.0, # Activation energy [J/mol]
    :n1       => 1.0,       # Reaction order of matrix
    :nu_char  => 0.2,       # Char yield (fixed parameter)
    :log10_A2 => 7.0,       # Pre-exponential factor (log10) for char oxidation
    :E2       => 160_000.0, # Activation energy [J/mol]
    :n2       => 1.5,       # Reaction order of char
    :m2       => 0.6,       # Oxygen dependency exponent
    :log10_A3 => 9.0,       # Pre-exponential factor (log10) for fiber oxidation
    :E3       => 250_000.0, # Activation energy [J/mol]
    :n3       => 2.0,       # Reaction order of fiber
    :m3       => 0.85,      # Oxygen dependency exponent
)

# =============================================================================
# Initial state values for simulation
# =============================================================================
const INITIAL_M    = 0.3    # Initial fraction of matrix
const INITIAL_F    = 0.7    # Initial fraction of fiber
const INITIAL_C    = 0.0    # No char present initially
const INITIAL_TEMP = 300.0  # Initial temperature [K]
const INITIAL_G1   = 0.0    # No volatile gases at start
const INITIAL_G2   = 0.0    # No oxidation gases at start

# =============================================================================
# Experimental setup (heating rate and oxygen level)
# Note: Heating rate β converted from K/min to K/s for consistency
# =============================================================================
const EXPERIMENTS = Dict(
    "Exp1" => (beta = 2.5 / 60.0,  PO2 = 0.21),  # Slow heating, air condition
    "Exp2" => (beta = 5.0 / 60.0,  PO2 = 0.21),  # Moderate heating, air condition
    "Exp3" => (beta = 10.0 / 60.0, PO2 = 0.21),  # Fast heating, air condition
    "Exp4" => (beta = 5.0 / 60.0,  PO2 = 0.05),  # Moderate heating, low oxygen
)

# Target final temperature for TGA simulation [K]
const T_FINAL = 1200.0