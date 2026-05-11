"""Fitness evaluation for codon-optimized DNA sequences."""

from __future__ import annotations

from dataclasses import dataclass

GC_MIN = 0.40
GC_MAX = 0.60
GC_PENALTY = 2000.0
REPEAT_PENALTY = 25.0


@dataclass
class FitnessEvaluator:
    """Score a codon sequence using E. coli usage, GC content, and repeat penalties."""

    codon_frequencies: dict[str, float]
    evaluation_count: int = 0

    def evaluate(self, codons: list[str]) -> float:
        self.evaluation_count += 1
        score = sum(self.codon_frequencies[codon] for codon in codons)

        gc_fraction = self._gc_fraction(codons)
        if gc_fraction < GC_MIN or gc_fraction > GC_MAX:
            score -= GC_PENALTY

        for index in range(1, len(codons)):
            if codons[index] == codons[index - 1]:
                score -= REPEAT_PENALTY

        return score

    @staticmethod
    def _gc_fraction(codons: list[str]) -> float:
        sequence = "".join(codons)
        if not sequence:
            return 0.0
        gc_count = sum(nucleotide in {"G", "C"} for nucleotide in sequence)
        return gc_count / len(sequence)

    def codons_to_dna(self, codons: list[str]) -> str:
        return "".join(codons)
