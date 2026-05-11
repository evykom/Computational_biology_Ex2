"""Run the codon optimization genetic algorithm."""

from __future__ import annotations

import json
from pathlib import Path

from fitness import FitnessEvaluator
from genetic_algorithm import EvolutionStrategy, GeneticAlgorithm, GeneticAlgorithmConfig

GENETIC_CODE_PATH = Path("genetic_code.json")
ECOLI_FREQUENCIES_PATH = Path("ecoli_frequencies.json")
INSULIN_TARGET_PATH = Path("insulin_target.txt")


def run_basic_genetic_algorithm() -> None:
    target_sequence = INSULIN_TARGET_PATH.read_text(encoding="utf-8").strip().upper()
    with GENETIC_CODE_PATH.open(encoding="utf-8") as genetic_code_file:
        genetic_code = json.load(genetic_code_file)
    with ECOLI_FREQUENCIES_PATH.open(encoding="utf-8") as frequencies_file:
        codon_frequencies = {codon: float(score) for codon, score in json.load(frequencies_file).items()}
    fitness_evaluator = FitnessEvaluator(codon_frequencies=codon_frequencies)

    algorithm = GeneticAlgorithm(
        target_sequence=target_sequence,
        genetic_code=genetic_code,
        fitness_evaluator=fitness_evaluator,
        config=GeneticAlgorithmConfig(),
        strategy=EvolutionStrategy.DARWINIAN_BASIC,
    )
    result = algorithm.run()

    best_dna = fitness_evaluator.codons_to_dna(result.best_codons)
    print("Strategy: basic genetic algorithm")
    print(f"Best DNA sequence: {best_dna}")
    print(f"Best fitness score: {result.best_fitness:.2f}")
    print(f"Fitness evaluations: {result.fitness_evaluations}")


if __name__ == "__main__":
    run_basic_genetic_algorithm()
