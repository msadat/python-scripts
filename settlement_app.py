"""
Settlement Prediction Tool - V1
Author: Dr. Sadat / ChatGPT

Purpose:
Estimate one-dimensional consolidation settlement for normally consolidated
or overconsolidated clay using classical geotechnical settlement equations.

Core assumptions in V1:
- 1D primary consolidation settlement
- Clay layer represented by average properties
- Stress increase is applied as an average increase across the compressible layer
- Units are user-controlled but must be consistent

Recommended US customary workflow:
- Thickness H: ft
- Stresses: psf
- Unit weight: pcf
- Settlement output: ft and inches
"""

from dataclasses import dataclass
from math import log10
from typing import Optional, List, Dict

#import matplotlib.pyplot as plt


@dataclass
class ClayLayer:
    """Represents one compressible clay layer."""
    name: str
    thickness_ft: float
    initial_void_ratio: float
    initial_effective_stress_psf: float
    stress_increase_psf: float
    compression_index_cc: Optional[float] = None
    recompression_index_cr: Optional[float] = None
    preconsolidation_stress_psf: Optional[float] = None

    def validate(self) -> None:
        if self.thickness_ft <= 0:
            raise ValueError("Layer thickness must be greater than zero.")
        if self.initial_void_ratio < 0:
            raise ValueError("Initial void ratio cannot be negative.")
        if self.initial_effective_stress_psf <= 0:
            raise ValueError("Initial effective stress must be greater than zero.")
        if self.stress_increase_psf < 0:
            raise ValueError("Stress increase cannot be negative.")
        if self.compression_index_cc is None and self.recompression_index_cr is None:
            raise ValueError("At least Cc or Cr must be provided.")


def settlement_normally_consolidated(layer: ClayLayer) -> float:
    """
    Settlement for normally consolidated clay:

        S = H * Cc / (1 + e0) * log10((sigma0' + Delta sigma) / sigma0')

    Returns settlement in ft.
    """
    layer.validate()

    if layer.compression_index_cc is None:
        raise ValueError("Compression index Cc is required for normally consolidated settlement.")

    sigma0 = layer.initial_effective_stress_psf
    delta_sigma = layer.stress_increase_psf
    final_stress = sigma0 + delta_sigma

    settlement_ft = (
        layer.thickness_ft
        * layer.compression_index_cc
        / (1 + layer.initial_void_ratio)
        * log10(final_stress / sigma0)
    )

    return settlement_ft


def settlement_overconsolidated(layer: ClayLayer) -> float:
    """
    Settlement for overconsolidated clay.

    Cases:
    1. Final stress <= preconsolidation stress:
       S = H * Cr / (1 + e0) * log10(sigma_f' / sigma0')

    2. Final stress > preconsolidation stress:
       S = H/(1+e0) * [Cr * log10(sigma_p' / sigma0')
                      + Cc * log10(sigma_f' / sigma_p')]

    Returns settlement in ft.
    """
    layer.validate()

    if layer.recompression_index_cr is None:
        raise ValueError("Recompression index Cr is required for overconsolidated settlement.")
    if layer.compression_index_cc is None:
        raise ValueError("Compression index Cc is required when final stress exceeds preconsolidation stress.")
    if layer.preconsolidation_stress_psf is None:
        raise ValueError("Preconsolidation stress is required for overconsolidated settlement.")

    sigma0 = layer.initial_effective_stress_psf
    sigma_p = layer.preconsolidation_stress_psf
    sigma_f = sigma0 + layer.stress_increase_psf
    h_factor = layer.thickness_ft / (1 + layer.initial_void_ratio)

    if sigma_p < sigma0:
        raise ValueError("Preconsolidation stress should be greater than or equal to initial effective stress.")

    if sigma_f <= sigma_p:
        settlement_ft = h_factor * layer.recompression_index_cr * log10(sigma_f / sigma0)
    else:
        settlement_ft = h_factor * (
            layer.recompression_index_cr * log10(sigma_p / sigma0)
            + layer.compression_index_cc * log10(sigma_f / sigma_p)
        )

    return settlement_ft


def settlement_to_inches(settlement_ft: float) -> float:
    """Convert settlement from ft to inches."""
    return settlement_ft * 12


def calculate_layer_settlement(layer: ClayLayer, condition: str = "normally_consolidated") -> Dict[str, float]:
    """
    Calculate settlement for a clay layer.

    condition options:
    - "normally_consolidated"
    - "overconsolidated"
    """
    if condition == "normally_consolidated":
        settlement_ft = settlement_normally_consolidated(layer)
    elif condition == "overconsolidated":
        settlement_ft = settlement_overconsolidated(layer)
    else:
        raise ValueError("condition must be 'normally_consolidated' or 'overconsolidated'.")

    return {
        "settlement_ft": settlement_ft,
        "settlement_in": settlement_to_inches(settlement_ft),
        "final_effective_stress_psf": layer.initial_effective_stress_psf + layer.stress_increase_psf,
    }


def calculate_total_settlement(layers: List[ClayLayer], condition: str = "normally_consolidated") -> Dict[str, float]:
    """Calculate total settlement for multiple clay layers."""
    total_ft = 0.0

    for layer in layers:
        result = calculate_layer_settlement(layer, condition=condition)
        total_ft += result["settlement_ft"]

    return {
        "total_settlement_ft": total_ft,
        "total_settlement_in": settlement_to_inches(total_ft),
    }


def sensitivity_stress_increase(
    base_layer: ClayLayer,
    stress_increases_psf: List[float],
    condition: str = "normally_consolidated",
) -> List[Dict[str, float]]:
    """Run settlement sensitivity for a range of stress increases."""
    results = []

    for delta_sigma in stress_increases_psf:
        test_layer = ClayLayer(
            name=base_layer.name,
            thickness_ft=base_layer.thickness_ft,
            initial_void_ratio=base_layer.initial_void_ratio,
            initial_effective_stress_psf=base_layer.initial_effective_stress_psf,
            stress_increase_psf=delta_sigma,
            compression_index_cc=base_layer.compression_index_cc,
            recompression_index_cr=base_layer.recompression_index_cr,
            preconsolidation_stress_psf=base_layer.preconsolidation_stress_psf,
        )

        result = calculate_layer_settlement(test_layer, condition=condition)
        results.append({
            "stress_increase_psf": delta_sigma,
            "settlement_ft": result["settlement_ft"],
            "settlement_in": result["settlement_in"],
        })

    return results


def plot_stress_sensitivity(results: List[Dict[str, float]]) -> None:
    """Plot stress increase vs settlement."""
    stress = [row["stress_increase_psf"] for row in results]
    settlement = [row["settlement_in"] for row in results]

    plt.figure(figsize=(8, 5))
    plt.plot(stress, settlement, marker="o")
    plt.xlabel("Stress Increase, Δσ (psf)")
    plt.ylabel("Settlement (inches)")
    plt.title("Settlement Sensitivity to Stress Increase")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Example based on a typical normally consolidated clay problem.
    # User should replace these values with project-specific parameters.

    example_layer = ClayLayer(
        name="NC Clay Layer",
        thickness_ft=30.0,
        initial_void_ratio=0.65,
        initial_effective_stress_psf=30.0 * 97.0,  # approximate sigma0' from 30 ft soil column at 97 pcf
        stress_increase_psf=200.0,
        compression_index_cc=0.29,
    )

    result = calculate_layer_settlement(example_layer, condition="normally_consolidated")

    print("Settlement Prediction - Normally Consolidated Clay")
    print(f"Layer: {example_layer.name}")
    print(f"Initial effective stress: {example_layer.initial_effective_stress_psf:,.1f} psf")
    print(f"Final effective stress: {result['final_effective_stress_psf']:,.1f} psf")
    print(f"Settlement: {result['settlement_ft']:.4f} ft")
    print(f"Settlement: {result['settlement_in']:.2f} inches")

    stress_values = [100, 200, 400, 600, 800]
    sensitivity_results = sensitivity_stress_increase(example_layer, stress_values)

    print("\nStress Sensitivity")
    for row in sensitivity_results:
        print(
            f"Δσ = {row['stress_increase_psf']:>5.0f} psf | "
            f"Settlement = {row['settlement_in']:.2f} in"
        )

    #plot_stress_sensitivity(sensitivity_results)
