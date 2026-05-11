"""Lifetime local optimization for codon sequences."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class LocalSearchConfig:
    budget: int = 5
    bad_codon_threshold: float = 5.0


class LocalSearch:
    def __init__(
        self,
        target_sequence: str,
        genetic_code: dict[str, list[str]],
        codon_frequencies: dict[str, float],
        config: LocalSearchConfig | None = None,
    ) -> None:
        self.target_sequence = target_sequence
        self.genetic_code = genetic_code
        self.codon_frequencies = codon_frequencies
        self.config = config or LocalSearchConfig()

    def improve(self, codons: list[str]) -> list[str]:
        current = codons[:]
        for _ in range(self.config.budget):
            repeat_index = self._find_repeat_index(current)
            if repeat_index is not None:
                current = self._fix_repeat(current, repeat_index)
                continue

            bad_codon_index = self._find_bad_codon_index(current)
            if bad_codon_index is None:
                break
            current = self._fix_bad_codon(current, bad_codon_index)

        return current

    def _find_repeat_index(self, codons: list[str]) -> int | None:
        for index in range(1, len(codons)):
            if codons[index] == codons[index - 1]:
                return index
        return None

    def _fix_repeat(self, codons: list[str], index: int) -> list[str]:
        amino_acid = self.target_sequence[index]
        neighbor_codon = codons[index - 1]
        alternatives = [
            codon
            for codon in self.genetic_code[amino_acid]
            if codon != neighbor_codon
        ]
        if not alternatives:
            return codons[:]

        best_codon = max(alternatives, key=self.codon_frequencies.__getitem__)
        updated = codons[:]
        updated[index] = best_codon
        return updated

    def _find_bad_codon_index(self, codons: list[str]) -> int | None:
        for index, codon in enumerate(codons):
            if self.codon_frequencies[codon] < self.config.bad_codon_threshold:
                return index
        return None

    def _fix_bad_codon(self, codons: list[str], index: int) -> list[str]:
        amino_acid = self.target_sequence[index]
        best_codon = max(
            self.genetic_code[amino_acid],
            key=self.codon_frequencies.__getitem__,
        )
        updated = codons[:]
        updated[index] = best_codon
        return updated
