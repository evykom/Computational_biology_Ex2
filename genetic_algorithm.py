"""Basic genetic algorithm for codon optimization."""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from enum import Enum

from fitness import FitnessEvaluator


class EvolutionStrategy(Enum):
    DARWINIAN_BASIC = "darwinian_basic"
    BALDWIN = "baldwin"
    LAMARCKIAN = "lamarckian"


@dataclass
class GeneticAlgorithmConfig:
    population_size: int = 120
    generations: int = 250
    crossover_rate: float = 0.85
    mutation_rate: float = 0.03
    tournament_size: int = 3
    elitism_count: int = 2
    random_seed: int | None = 42


@dataclass
class GeneticAlgorithmResult:
    best_codons: list[str]
    best_fitness: float
    fitness_evaluations: int
    generation_history: list[tuple[float, float]] = field(default_factory=list)


class GeneticAlgorithm:
    def __init__(
        self,
        target_sequence: str,
        genetic_code: dict[str, list[str]],
        fitness_evaluator: FitnessEvaluator,
        config: GeneticAlgorithmConfig | None = None,
        strategy: EvolutionStrategy = EvolutionStrategy.DARWINIAN_BASIC,
        local_search_budget: int = 0,
    ) -> None:
        self.target_sequence = target_sequence
        self.genetic_code = genetic_code
        self.fitness_evaluator = fitness_evaluator
        self.config = config or GeneticAlgorithmConfig()
        self.strategy = strategy
        self.local_search_budget = local_search_budget

        if self.config.random_seed is not None:
            random.seed(self.config.random_seed)

        self._validate_target_sequence()

    def _validate_target_sequence(self) -> None:
        for amino_acid in self.target_sequence:
            if amino_acid not in self.genetic_code:
                raise ValueError(f"No codons available for amino acid '{amino_acid}'")

    def run(self) -> GeneticAlgorithmResult:
        population = [
            self._create_random_individual()
            for _ in range(self.config.population_size)
        ]
        fitnesses = [self._evaluate_individual(individual) for individual in population]

        best_index = max(range(len(population)), key=lambda index: fitnesses[index])
        best_codons = population[best_index][:]
        best_fitness = fitnesses[best_index]
        generation_history: list[tuple[float, float]] = []
        print(f"Generation 0: best score {best_fitness:.2f}")

        for generation in range(1, self.config.generations + 1):
            average_fitness = sum(fitnesses) / len(fitnesses)
            generation_history.append((average_fitness, best_fitness))

            next_population: list[list[str]] = []
            elite_indices = sorted(
                range(len(population)),
                key=lambda index: fitnesses[index],
                reverse=True,
            )[: self.config.elitism_count]
            for elite_index in elite_indices:
                next_population.append(population[elite_index][:])

            while len(next_population) < self.config.population_size:
                parent_a = self._tournament_select(population, fitnesses)
                parent_b = self._tournament_select(population, fitnesses)
                child_a, child_b = self._crossover(parent_a, parent_b)
                child_a = self._mutate(child_a)
                child_b = self._mutate(child_b)
                next_population.append(child_a)
                if len(next_population) < self.config.population_size:
                    next_population.append(child_b)

            population = next_population
            fitnesses = [self._evaluate_individual(individual) for individual in population]

            generation_best_index = max(range(len(population)), key=lambda index: fitnesses[index])
            generation_best_fitness = fitnesses[generation_best_index]
            if generation_best_fitness > best_fitness:
                best_fitness = generation_best_fitness
                best_codons = population[generation_best_index][:]

            print(f"Generation {generation}: best score {best_fitness:.2f}")

        return GeneticAlgorithmResult(
            best_codons=best_codons,
            best_fitness=best_fitness,
            fitness_evaluations=self.fitness_evaluator.evaluation_count,
            generation_history=generation_history,
        )

    def _create_random_individual(self) -> list[str]:
        return [
            random.choice(self.genetic_code[amino_acid])
            for amino_acid in self.target_sequence
        ]

    def _evaluate_individual(self, codons: list[str]) -> float:
        return self.fitness_evaluator.evaluate(codons)

    def _tournament_select(
        self,
        population: list[list[str]],
        fitnesses: list[float],
    ) -> list[str]:
        contenders = random.sample(range(len(population)), self.config.tournament_size)
        winner_index = max(contenders, key=lambda index: fitnesses[index])
        return population[winner_index][:]

    def _crossover(
        self,
        parent_a: list[str],
        parent_b: list[str],
    ) -> tuple[list[str], list[str]]:
        if random.random() >= self.config.crossover_rate:
            return parent_a[:], parent_b[:]

        crossover_point = random.randint(1, len(parent_a) - 1)
        child_a = parent_a[:crossover_point] + parent_b[crossover_point:]
        child_b = parent_b[:crossover_point] + parent_a[crossover_point:]
        return child_a, child_b

    def _mutate(self, codons: list[str]) -> list[str]:
        mutated = codons[:]
        for index, amino_acid in enumerate(self.target_sequence):
            if random.random() < self.config.mutation_rate:
                synonymous_codons = self.genetic_code[amino_acid]
                alternatives = [codon for codon in synonymous_codons if codon != mutated[index]]
                if alternatives:
                    mutated[index] = random.choice(alternatives)
        return mutated
