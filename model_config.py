# -*- coding: utf-8 -*-
"""
================================================================================
Lettuce UVA Model Shared Configuration Module
================================================================================
Version: v2.0
Date: 2025-12-05

Description:
---------
This module centralizes all shared configuration settings to avoid
desynchronization caused by repeated definitions across different files.
Includes:
  1. Environmental base settings (plant factory environment)
  2. Treatment group configurations (UVA irradiation time and intensity)
  3. Experimental target data (observed values)
  4. Dataset split (training/validation/test)
  5. ODE solver settings
  6. Simulation settings
  7. Night inhibition parameter prior ranges

Usage:
-----
    from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS

Notes:
---------
  - ODE solver settings (RK45, max_step=60) have been calibrated, do not modify!
  - After modifying any settings, please re-run run_validation.py to verify results

Dependent scripts:
---------
  - simulate_uva_model.py
  - run_validation.py
  - run_scenarios.py
  - optimize_parameters.py
  - sensitivity_analysis.py
================================================================================
"""

# ==============================================================================
# Part 1: Environmental Base Settings (Plant Factory Environment)
# ==============================================================================
ENV_BASE = {
    # --- Light settings ---
    'light_on_hour': 6,       # Light start time (06:00)
    'light_off_hour': 22,     # Light end time (22:00) - 16 hours of light total
    # PPFD 130 μmol/m²/s → PAR power = 130 × 0.219 = 28.5 W/m²
    # Sun model assumes PAR = 0.5 × I, so equivalent shortwave radiation I = 28.5 / 0.5 = 57 W/m²
    'I_day': 57,              # Daytime shortwave radiation [W/m²] (physically calculated value)

    # --- Temperature settings ---
    'T_day': 25,              # Daytime temperature [°C]
    'T_night': 18,            # Nighttime temperature [°C]

    # --- CO2 settings ---
    'CO2_day': 1200,          # Daytime CO2 concentration [ppm] (Sun model calibrated value)
    'CO2_night': 1200,        # Nighttime CO2 concentration [ppm]

    # --- Humidity settings ---
    'RH_day': 0.70,           # Daytime relative humidity [0-1]
    'RH_night': 0.85,         # Nighttime relative humidity [0-1]

    # --- Plant density ---
    'plant_density': 36,      # Plant density [plants/m²]
}


# ==============================================================================
# Part 2: Treatment Group Configurations
# ==============================================================================
# Description:
# - uva_start_day / uva_end_day: Days after sowing
# - uva_hour_on / uva_hour_off: Daily UVA irradiation start/end time
# - Simulation starts from day 14 after sowing (transplant day), lasts 21 days, ends on day 35

TREATMENT_CONFIGS = {
    # --- Control group ---
    'CK': {
        'uva_on': False,
        'description': 'Control group (no UVA)',
    },

    # --- Standard daytime irradiation (optimal treatment) ---
    'L6D6': {
        'uva_on': True,
        'uva_intensity': 11.0,          # UVA intensity [W/m²] (v10.7: actual LED power)
        'uva_start_day': 29,            # Start on day 29 after sowing
        'uva_end_day': 35,              # End on day 35 after sowing (6 days total)
        'uva_hour_on': 10,              # Start at 10:00
        'uva_hour_off': 16,             # End at 16:00 (6 hours daytime)
        'description': 'Low dose daytime (6h/day, 6 days)',
    },

    # --- Nighttime irradiation (circadian rhythm disruption effect) ---
    'L6D6-N': {
        'uva_on': True,
        'uva_intensity': 11.0,          # v10.7: actual LED power
        'uva_start_day': 29,
        'uva_end_day': 35,
        'uva_hour_on': 22,              # Start at 22:00
        'uva_hour_off': 4,              # End at 04:00 (6 hours overnight)
        'description': 'Low dose nighttime (6h/night, 6 days)',
    },

    # --- High dose stress (damage mechanism) ---
    'H12D3': {
        'uva_on': True,
        'uva_intensity': 11.0,          # v10.7: actual LED power
        'uva_start_day': 32,            # Start on day 32 after sowing
        'uva_end_day': 35,              # End on day 35 after sowing (3 days total)
        'uva_hour_on': 6,               # Start at 06:00
        'uva_hour_off': 18,             # End at 18:00 (12 hours)
        'description': 'High dose stress (12h/day, 3 days)',
    },

    # --- Very low dose long-term (adaptation effect) ---
    'VL3D12': {
        'uva_on': True,
        'uva_intensity': 11.0,          # v10.7: actual LED power
        'uva_start_day': 23,            # Start on day 23 after sowing
        'uva_end_day': 35,              # End on day 35 after sowing (12 days total)
        'uva_hour_on': 10,              # Start at 10:00
        'uva_hour_off': 13,             # End at 13:00 (3 hours)
        'description': 'Very low dose long-term (3h/day, 12 days)',
    },

    # --- Low dose long-term (morphological effect dominant) ---
    'L6D12': {
        'uva_on': True,
        'uva_intensity': 11.0,          # v10.7: actual LED power
        'uva_start_day': 23,
        'uva_end_day': 35,
        'uva_hour_on': 10,
        'uva_hour_off': 16,
        'description': 'Low dose long-term (6h/day, 12 days) - morphological effect accumulation',
    },

    # ==========================================================================
    # Validation group (3-day gradient experiment - Day 32-35)
    # ==========================================================================
    # Validation control group (same settings as training CK, but independent observed values)
    'CK_val': {
        'uva_on': False,
        'description': 'Validation control group (no UVA)',
    },

    # Very low daily dose 3 days
    'VL3D3': {
        'uva_on': True,
        'uva_intensity': 11.0,
        'uva_start_day': 32,
        'uva_end_day': 35,
        'uva_hour_on': 10,
        'uva_hour_off': 13,             # 3h/day
        'description': 'Very low daily dose (3h/day, 3 days)',
    },

    # Low daily dose 3 days
    'L6D3': {
        'uva_on': True,
        'uva_intensity': 11.0,
        'uva_start_day': 32,
        'uva_end_day': 35,
        'uva_hour_on': 10,
        'uva_hour_off': 16,             # 6h/day
        'description': 'Low daily dose (6h/day, 3 days)',
    },

    # Medium daily dose 3 days
    'M9D3': {
        'uva_on': True,
        'uva_intensity': 11.0,
        'uva_start_day': 32,
        'uva_end_day': 35,
        'uva_hour_on': 7,
        'uva_hour_off': 16,             # 9h/day
        'description': 'Medium daily dose (9h/day, 3 days)',
    },

    # High daily dose 3 days (validation group - same settings as training H12D3 but independent observed values)
    'H12D3_val': {
        'uva_on': True,
        'uva_intensity': 11.0,
        'uva_start_day': 32,
        'uva_end_day': 35,
        'uva_hour_on': 6,
        'uva_hour_off': 18,             # 12h/day
        'description': 'High daily dose (12h/day, 3 days) - validation group',
    },

    # Very high daily dose 3 days
    'VH15D3': {
        'uva_on': True,
        'uva_intensity': 11.0,
        'uva_start_day': 32,
        'uva_end_day': 35,
        'uva_hour_on': 5,
        'uva_hour_off': 20,             # 15h/day
        'description': 'Very high daily dose (15h/day, 3 days)',
    },

    # --- Sinusoidal gradual irradiation (dynamic mode) ---
    # Irradiation intensity varies sinusoidally with time, from 0 to max then back to 0
    # Period 12 hours (06:00-18:00), total effective dose equivalent to 6 hours continuous irradiation
    'SIN': {
        'uva_on': True,
        'uva_mode': 'sinusoidal',        # Dynamic irradiation mode
        'uva_intensity': 22.0,           # Maximum intensity [W/m²]
        'uva_start_day': 23,             # Start on day 23 after sowing
        'uva_end_day': 35,               # End on day 35 after sowing (12 days total)
        'uva_hour_on': 6,                # Start at 06:00
        'uva_hour_off': 18,              # End at 18:00 (12 hour period)
        'description': 'Sinusoidal gradual irradiation (6h equivalent/day, 12 days)',
    },

    # --- Intermittent irradiation (dynamic mode) ---
    # 30 minutes on/30 minutes off, total 6 hours effective irradiation distributed over 12 hours
    'INT': {
        'uva_on': True,
        'uva_mode': 'intermittent',      # Dynamic irradiation mode
        'uva_intensity': 22.0,           # Intensity during irradiation [W/m²]
        'uva_start_day': 23,             # Start on day 23 after sowing
        'uva_end_day': 35,               # End on day 35 after sowing (12 days total)
        'uva_hour_on': 6,                # Start at 06:00
        'uva_hour_off': 18,              # End at 18:00
        'intermittent_on_min': 30,       # On duration [minutes]
        'intermittent_off_min': 30,      # Off duration [minutes]
        'description': 'Intermittent irradiation (30min on/off, 6h equivalent/day, 12 days)',
    },
}


# ==============================================================================
# Part 3: Experimental Target Data (Observed Values)
# ==============================================================================
# Unit description:
# - FW: Fresh Weight [g/plant]
# - Anth: Anthocyanin concentration [ppm = mg/kg FW]

TARGETS = {
    # Data update (2026-01-09)
    # Fresh weight (g/plant):
    #   Ref=87, 6h/6d(D)=91.4, 6h/6d(N)=80.8, 12h/3d=60.6, 3h/12d=67, 6h/12d=60.4
    # Anthocyanin concentration (ppm = mg/kg FW):
    #   Unit corrected: original mg/100g → now mg/kg (×10)
    'CK': {'FW': 87.0, 'Anth': 433.0},       # Anth STD: 9.6
    'L6D6': {'FW': 91.4, 'Anth': 494.0},     # Anth STD: 6.8
    'L6D6-N': {'FW': 80.8, 'Anth': 493.0},   # Anth STD: 7.4
    'H12D3': {'FW': 60.6, 'Anth': 651.0},    # Anth STD: 14.6
    'VL3D12': {'FW': 67.0, 'Anth': 482.0},   # Anth STD: 2.7
    'L6D12': {'FW': 60.4, 'Anth': 518.0},    # Anth STD: 3.4
    # --- Validation group (3-day gradient experiment) ---
    'CK_val': {'FW': 85.2, 'Anth': 413.0},    # Validation control group
    'VL3D3': {'FW': 89.0, 'Anth': 437.0},     # 3h/day × 3d
    'L6D3': {'FW': 92.2, 'Anth': 468.0},      # 6h/day × 3d
    'M9D3': {'FW': 83.8, 'Anth': 539.0},      # 9h/day × 3d
    'H12D3_val': {'FW': 62.2, 'Anth': 657.0}, # 12h/day × 3d (validation group)
    'VH15D3': {'FW': 51.3, 'Anth': 578.0},    # 15h/day × 3d
}


# ==============================================================================
# Part 4: Dataset Split (Training/Validation/Test)
# ==============================================================================
# Description:
# - Training set: Used for parameter estimation (6 groups)
# - Validation set: Used for 3-day gradient experiment validation (6 groups)
# - Test set: SIN, INT (dynamic modes, not yet implemented)

DATASET_SPLIT = {
    'train': ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12'],
    'validation': ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3'],
    'test': [],
}


# ==============================================================================
# Part 5: ODE Solver Settings (Important: Do not modify!)
# ==============================================================================
# Description:
# - method='RK45' and max_step=60 have been calibrated, modifying will cause inconsistent results
# - LSODA produces unstable results at different step sizes, so RK45 is used

ODE_SETTINGS = {
    'method': 'RK45',         # Solver method (do not change!)
    'max_step': 60,           # Maximum step size [seconds] (do not change!)
}


# ==============================================================================
# Part 6: Simulation Settings
# ==============================================================================
SIMULATION = {
    'days': 21,               # Simulation days (0-21 days after transplant)
    'transplant_offset': 14,  # Transplant offset (transplant on day 14 after sowing)
    'initial_fw_g': 10,       # Initial fresh weight at transplant [g/plant]
}


# ==============================================================================
# Part 7: Night Inhibition Parameter Prior Ranges (Based on Literature)
# ==============================================================================
# References:
# - Bennie et al. (2016) J. Ecology DOI: 10.1111/1365-2745.12551
# - Deng et al. (2025) Biology DOI: 10.3390/biology14050571
# - Harmer (2009) Annu. Rev. Plant Biol.
# - Covington et al. (2008) Genome Biol.

NIGHT_INHIBITION_SCENARIOS = {
    'Low': {
        'night_uva_base_inhibition': 0.08,
        'circadian_inhibition_decay': 3.2e-5,  # Half-life approximately 6 hours
        'description': 'Low inhibition scenario (fast recovery)',
    },
    'Medium': {
        'night_uva_base_inhibition': 0.12,
        'circadian_inhibition_decay': 1.6e-5,  # Half-life approximately 12 hours
        'description': 'Medium inhibition scenario (literature median)',
    },
    'High': {
        'night_uva_base_inhibition': 0.18,
        'circadian_inhibition_decay': 8.0e-6,  # Half-life approximately 24 hours
        'description': 'High inhibition scenario (slow recovery)',
    },
}


# ==============================================================================
# Part 8: Parameter Prior Ranges (For Sensitivity Analysis and Optimization)
# ==============================================================================
PARAM_PRIORS = {
    'night_uva_base_inhibition': {
        'range': (0.08, 0.18),
        'default': 0.12,
        'unit': '[-]',
        'source': 'Bennie et al. 2016; Deng et al. 2025',
    },
    'circadian_inhibition_decay': {
        'range': (8.0e-6, 3.2e-5),
        'default': 1.6e-5,
        'unit': '[1/s]',
        'source': 'Covington et al. 2008; Harmer 2009',
    },
    'adaptation_rate': {
        'range': (1e-7, 5e-7),
        'default': 3e-7,
        'unit': '[1/Ws]',
        'source': 'Hideg et al. 2013',
    },
    'max_adaptation': {
        'range': (0.4, 0.8),
        'default': 0.65,
        'unit': '[-]',
        'source': 'Jansen et al. 1998',
    },
}


# ==============================================================================
# Utility Functions
# ==============================================================================

def get_env_for_treatment(treatment: str) -> dict:
    """
    Get the complete environment settings for a specific treatment group.

    Merges settings from ENV_BASE and TREATMENT_CONFIGS.

    Parameters:
    -----------
    treatment : str
        Treatment group code (e.g., 'CK', 'L6D6', etc.)

    Returns:
    --------
    dict
        Complete environment settings dictionary
    """
    env = ENV_BASE.copy()
    if treatment in TREATMENT_CONFIGS:
        env.update(TREATMENT_CONFIGS[treatment])
    return env


def get_target(treatment: str) -> dict:
    """
    Get the target values for a specific treatment group.

    Parameters:
    -----------
    treatment : str
        Treatment group code

    Returns:
    --------
    dict
        Target values dictionary {'FW': float, 'Anth': float}
    """
    return TARGETS.get(treatment, {'FW': 0, 'Anth': 0})


def get_treatments_by_dataset(dataset: str) -> list:
    """
    Get the list of treatment groups for a specific dataset.

    Parameters:
    -----------
    dataset : str
        Dataset name ('train', 'validation', 'test')

    Returns:
    --------
    list
        List of treatment group codes
    """
    return DATASET_SPLIT.get(dataset, [])


def print_config_summary():
    """
    Print configuration summary.

    Used for quickly reviewing current model settings.
    """
    print("=" * 70)
    print("Lettuce UVA Model Configuration Summary")
    print("=" * 70)

    print("\nEnvironment Settings:")
    print(f"  Light: {ENV_BASE['light_on_hour']}:00-{ENV_BASE['light_off_hour']}:00, "
          f"{ENV_BASE['I_day']} W/m²")
    print(f"  Temperature: Daytime {ENV_BASE['T_day']}°C / Nighttime {ENV_BASE['T_night']}°C")
    print(f"  CO2: Daytime {ENV_BASE['CO2_day']} / Nighttime {ENV_BASE['CO2_night']} ppm")
    print(f"  Plant density: {ENV_BASE['plant_density']} plants/m²")

    print("\nTreatment Groups:")
    for name, config in TREATMENT_CONFIGS.items():
        target = TARGETS.get(name, {})
        print(f"  {name}: {config.get('description', '')}")
        if target:
            print(f"       Target: FW={target['FW']}g, Anth={target['Anth']}ppm")

    print("\nDataset Split:")
    print(f"  Training set: {DATASET_SPLIT['train']}")
    print(f"  Validation set: {DATASET_SPLIT['validation']}")
    print(f"  Test set: {DATASET_SPLIT['test']} (not yet implemented)")

    print("\nODE Settings:")
    print(f"  Method: {ODE_SETTINGS['method']}")
    print(f"  Max step: {ODE_SETTINGS['max_step']} seconds")

    print("\nSimulation Settings:")
    print(f"  Simulation days: {SIMULATION['days']} days")
    print(f"  Initial fresh weight: {SIMULATION['initial_fw_g']} g/plant")

    print("=" * 70)


if __name__ == '__main__':
    print_config_summary()
