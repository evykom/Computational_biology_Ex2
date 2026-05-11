"""Run the codon optimization genetic algorithm."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from src.fitness import FitnessEvaluator
from src.genetic_algorithm import (
    EvolutionStrategy,
    GeneticAlgorithm,
    GeneticAlgorithmConfig,
    GeneticAlgorithmResult,
)

GENETIC_CODE_PATH = Path("genetic_code.json")
ECOLI_FREQUENCIES_PATH = Path("ecoli_frequencies.json")
INSULIN_TARGET_PATH = Path("insulin_target.txt")
CSV_OUTPUT_DIR = Path("csv_files")

STRATEGY_CHOICES = {
    "basic": EvolutionStrategy.DARWINIAN_BASIC,
    "baldwin": EvolutionStrategy.BALDWIN,
    "lamarckian": EvolutionStrategy.LAMARCKIAN,
}


def load_problem_data() -> tuple[str, dict[str, list[str]], dict[str, float]]:
    target_sequence = INSULIN_TARGET_PATH.read_text(encoding="utf-8").strip().upper()
    with GENETIC_CODE_PATH.open(encoding="utf-8") as genetic_code_file:
        genetic_code = json.load(genetic_code_file)
    with ECOLI_FREQUENCIES_PATH.open(encoding="utf-8") as frequencies_file:
        codon_frequencies = {
            codon: float(score) for codon, score in json.load(frequencies_file).items()
        }
    return target_sequence, genetic_code, codon_frequencies


def write_generation_csv(csv_path: Path, result: GeneticAlgorithmResult) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            [
                "generation",
                "average_fitness",
                "generation_best_fitness",
                "best_fitness",
            ]
        )
        for record in result.generation_history:
            writer.writerow(
                [
                    record.generation,
                    f"{record.average_fitness:.2f}",
                    f"{record.generation_best_fitness:.2f}",
                    f"{record.best_fitness:.2f}",
                ]
            )


def run_strategy(
    strategy: EvolutionStrategy,
    target_sequence: str,
    genetic_code: dict[str, list[str]],
    codon_frequencies: dict[str, float],
    config: GeneticAlgorithmConfig,
    write_csv: bool,
) -> None:
    fitness_evaluator = FitnessEvaluator(codon_frequencies=codon_frequencies)
    algorithm = GeneticAlgorithm(
        target_sequence=target_sequence,
        genetic_code=genetic_code,
        fitness_evaluator=fitness_evaluator,
        config=config,
        strategy=strategy,
    )
    result = algorithm.run()

    best_dna = fitness_evaluator.codons_to_dna(result.best_codons)
    print(f"Strategy: {strategy.value}")
    print(f"Best DNA sequence: {best_dna}")
    print(f"Best fitness score: {result.best_fitness:.2f}")
    print(f"Fitness evaluations: {result.fitness_evaluations}")
    if write_csv:
        CSV_OUTPUT_DIR.mkdir(exist_ok=True)
        csv_path = CSV_OUTPUT_DIR / f"{strategy.value}_generations.csv"
        write_generation_csv(csv_path, result)
        print(f"Generation results written to: {csv_path}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description="Codon optimization genetic algorithm")
    parser.add_argument(
        "--strategy",
        choices=[*STRATEGY_CHOICES.keys(), "all"],
        default="basic",
        help="Evolution strategy to run",
    )
    parser.add_argument(
        "--local-search-budget",
        type=int,
        default=5,
        help="Number of lifetime local-search fixes for Baldwin and Lamarckian runs",
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Write per-generation results to a CSV file",
    )
    args = parser.parse_args()

    target_sequence, genetic_code, codon_frequencies = load_problem_data()
    config = GeneticAlgorithmConfig(local_search_budget=args.local_search_budget)

    strategies = (
        list(STRATEGY_CHOICES.values())
        if args.strategy == "all"
        else [STRATEGY_CHOICES[args.strategy]]
    )

    for strategy in strategies:
        run_strategy(
            strategy=strategy,
            target_sequence=target_sequence,
            genetic_code=genetic_code,
            codon_frequencies=codon_frequencies,
            config=config,
            write_csv=args.csv,
        )


if __name__ == "__main__":
    main()
