import math
import os
import random
from dataclasses import dataclass, field
from typing import List, Tuple

from SB_SupportFunction import read_input
from evaluate_all_sols import evaluate_all_sols_check


VALUE_HEAVY = 40
RANDOM_SEED = 2222


@dataclass
class AHADiscreteParams:
    population_size: int
    max_iterations: int
    migration_period: int | None = None
    tournament_size: int = 3
    guided_intensity: Tuple[float, float] = (0.2, 0.6)
    territorial_intensity: Tuple[float, float] = (0.05, 0.2)


@dataclass
class AHADiscreteState:
    population: List[List[int]] = field(default_factory=list)
    fitness: List[float] = field(default_factory=list)
    evaluations: List[Tuple] = field(default_factory=list)
    visit_table: List[List[int]] = field(default_factory=list)
    best_solution: List[int] = field(default_factory=list)
    best_fitness: float = float("inf")
    best_evaluation: Tuple | None = None
    best_history: List[float] = field(default_factory=list)
    fitness_evaluations: int = 0


def _ensure_directory(path: str) -> None:
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)


class AHADiscrete:
    """Artificial Hummingbird Algorithm tailored for permutation-based problems."""

    def __init__(self, name_path_input: str, params: AHADiscreteParams):
        if params.population_size < 2:
            raise ValueError("population_size must be at least 2")

        self.params = params
        self.name_path_input = name_path_input
        self.df_item_pool, _ = read_input(name_path_input)
        self.num_items = self.df_item_pool.shape[0]
        self.heavy_item_set = set(self.df_item_pool[self.df_item_pool["weight"] >= VALUE_HEAVY].index)
        self.state = AHADiscreteState()
        self.params.migration_period = (
            params.migration_period
            if params.migration_period is not None
            else max(1, 2 * params.population_size)
        )

    # ------------------------------------------------------------------
    # Initialisation and evaluation helpers
    # ------------------------------------------------------------------
    def _random_solution(self) -> List[int]:
        solution = list(range(self.num_items))
        random.shuffle(solution)
        return solution

    def _evaluate(self, solution: List[int]) -> Tuple[float, Tuple]:
        evaluation = evaluate_all_sols_check(solution, self.df_item_pool, self.heavy_item_set, self.name_path_input)
        tardiness = evaluation[1]
        self.state.fitness_evaluations += 1
        return tardiness, evaluation

    def _initialise_population(self) -> None:
        population = []
        fitness = []
        evaluations = []
        visit_table = [[0 for _ in range(self.params.population_size)] for _ in range(self.params.population_size)]

        for _ in range(self.params.population_size):
            candidate = self._random_solution()
            tardiness, evaluation = self._evaluate(candidate)
            population.append(candidate)
            fitness.append(tardiness)
            evaluations.append(evaluation)

            if tardiness < self.state.best_fitness:
                self.state.best_fitness = tardiness
                self.state.best_solution = candidate.copy()
                self.state.best_evaluation = evaluation

        self.state.population = population
        self.state.fitness = fitness
        self.state.evaluations = evaluations
        self.state.visit_table = visit_table
        self.state.best_history = [self.state.best_fitness]

    # ------------------------------------------------------------------
    # Selection helpers
    # ------------------------------------------------------------------
    def _tournament_select(self, exclude_index: int) -> int:
        indices = [idx for idx in range(self.params.population_size) if idx != exclude_index]
        sample_size = min(self.params.tournament_size, len(indices))
        sampled = random.sample(indices, sample_size)
        return min(sampled, key=lambda idx: self.state.fitness[idx])

    # ------------------------------------------------------------------
    # Movement operators
    # ------------------------------------------------------------------
    def _guided_foraging_move(self, current: List[int], target: List[int]) -> List[int]:
        differing_positions = [idx for idx, value in enumerate(current) if value != target[idx]]
        if not differing_positions:
            return current.copy()

        lower, upper = self.params.guided_intensity
        intensity_fraction = random.uniform(lower, upper)
        steps = max(1, int(len(differing_positions) * intensity_fraction))

        neighbour = current.copy()
        random.shuffle(differing_positions)
        for pos in differing_positions[:steps]:
            desired_value = target[pos]
            swap_index = neighbour.index(desired_value)
            neighbour[pos], neighbour[swap_index] = neighbour[swap_index], neighbour[pos]
        return neighbour

    def _territorial_foraging_move(self, current: List[int]) -> List[int]:
        lower, upper = self.params.territorial_intensity
        fraction = random.uniform(lower, upper)
        swaps = max(1, int(self.num_items * fraction))

        neighbour = current.copy()
        for _ in range(swaps):
            i, j = random.sample(range(self.num_items), 2)
            neighbour[i], neighbour[j] = neighbour[j], neighbour[i]
        return neighbour

    # ------------------------------------------------------------------
    # Visit table updates
    # ------------------------------------------------------------------
    def _visit_success(self, hummingbird_index: int) -> None:
        for j in range(self.params.population_size):
            if j != hummingbird_index:
                self.state.visit_table[hummingbird_index][j] += 1

    def _visit_failure(self, hummingbird_index: int) -> None:
        for j in range(self.params.population_size):
            if j != hummingbird_index:
                self.state.visit_table[hummingbird_index][j] += 1

    def _visit_after_migration(self, migrated_index: int) -> None:
        for j in range(self.params.population_size):
            if j == migrated_index:
                continue
            self.state.visit_table[migrated_index][j] = 0
            self.state.visit_table[j][migrated_index] = max(self.state.visit_table[j]) + 1

    # ------------------------------------------------------------------
    # Evaluation and replacement
    # ------------------------------------------------------------------
    def _evaluate_and_replace(self, index: int, candidate: List[int]) -> None:
        tardiness, evaluation = self._evaluate(candidate)
        if tardiness < self.state.fitness[index]:
            self.state.population[index] = candidate
            self.state.fitness[index] = tardiness
            self.state.evaluations[index] = evaluation
            self._visit_success(index)

            if tardiness < self.state.best_fitness:
                self.state.best_fitness = tardiness
                self.state.best_solution = candidate.copy()
                self.state.best_evaluation = evaluation
        else:
            self._visit_failure(index)

    # ------------------------------------------------------------------
    # Main optimisation loop
    # ------------------------------------------------------------------
    def run(self) -> Tuple[List[int], float, Tuple | None]:
        self._initialise_population()

        for iteration in range(self.params.max_iterations):
            # Guided foraging
            for idx in range(self.params.population_size):
                guide_idx = self._tournament_select(idx)
                candidate = self._guided_foraging_move(self.state.population[idx], self.state.population[guide_idx])
                self._evaluate_and_replace(idx, candidate)

            # Territorial foraging
            for idx in range(self.params.population_size):
                candidate = self._territorial_foraging_move(self.state.population[idx])
                self._evaluate_and_replace(idx, candidate)

            # Migration
            if (iteration + 1) % self.params.migration_period == 0:
                worst_idx = max(range(self.params.population_size), key=lambda i: self.state.fitness[i])
                new_solution = self._random_solution()
                tardiness, evaluation = self._evaluate(new_solution)
                self.state.population[worst_idx] = new_solution
                self.state.fitness[worst_idx] = tardiness
                self.state.evaluations[worst_idx] = evaluation
                self._visit_after_migration(worst_idx)

                if tardiness < self.state.best_fitness:
                    self.state.best_fitness = tardiness
                    self.state.best_solution = new_solution.copy()
                    self.state.best_evaluation = evaluation

            self.state.best_history.append(self.state.best_fitness)
            print(
                f"Iteration {iteration + 1}: best tardiness = {self.state.best_fitness:.3f}, "
                f"evaluations = {self.state.fitness_evaluations}"
            )

        return self.state.best_solution, self.state.best_fitness, self.state.best_evaluation

    def plot_best_history(self, output_path: str = "plots/aha_best_history.svg") -> str:
        if not self.state.best_history:
            raise ValueError("No optimisation history recorded; run the solver first.")

        _ensure_directory(output_path)

        width, height = 900, 480
        margin = 60
        values = self.state.best_history
        iterations = list(range(len(values)))

        min_val = min(values)
        max_val = max(values)
        if math.isclose(min_val, max_val):
            max_val = min_val + 1.0

        x_span = max(1, len(iterations) - 1)
        y_span = max_val - min_val

        def scale_x(index: float) -> float:
            return margin + (index / x_span) * (width - 2 * margin)

        def scale_y(value: float) -> float:
            norm = (value - min_val) / y_span
            return height - margin - norm * (height - 2 * margin)

        points = [f"{scale_x(i):.2f},{scale_y(v):.2f}" for i, v in zip(iterations, values)]
        path_commands = [f"M {points[0]}"] + [f"L {point}" for point in points[1:]]

        y_ticks = _generate_ticks(min_val, max_val, 5)
        x_ticks = _generate_ticks(0, len(values) - 1, 5)

        svg_lines = [
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
            "  <style>text { font-family: sans-serif; font-size: 14px; }</style>",
            f"  <rect width='{width}' height='{height}' fill='white' stroke='black' stroke-width='1' />",
            f"  <line x1='{margin}' y1='{height - margin}' x2='{width - margin}' y2='{height - margin}' stroke='black' stroke-width='2' />",
            f"  <line x1='{margin}' y1='{height - margin}' x2='{margin}' y2='{margin}' stroke='black' stroke-width='2' />",
        ]

        # X ticks and labels
        for tick in x_ticks:
            x = scale_x(tick)
            svg_lines.append(
                f"  <line x1='{x:.2f}' y1='{height - margin}' x2='{x:.2f}' y2='{height - margin + 10}' stroke='black' />"
            )
            svg_lines.append(
                f"  <text x='{x:.2f}' y='{height - margin + 30}' text-anchor='middle'>{int(round(tick))}</text>"
            )

        # Y ticks and labels
        for tick in y_ticks:
            y = scale_y(tick)
            svg_lines.append(
                f"  <line x1='{margin - 10}' y1='{y:.2f}' x2='{margin}' y2='{y:.2f}' stroke='black' />"
            )
            svg_lines.append(
                f"  <text x='{margin - 15}' y='{y + 5:.2f}' text-anchor='end'>{tick:.2f}</text>"
            )

        svg_lines.append(
            f"  <path d='{' '.join(path_commands)}' fill='none' stroke='#1976d2' stroke-width='3' stroke-linejoin='round' stroke-linecap='round' />"
        )
        svg_lines.append(
            f"  <text x='{width / 2}' y='{margin / 2}' text-anchor='middle'>Best tardiness per iteration</text>"
        )
        svg_lines.append(
            "</svg>"
        )

        svg_content = "\n".join(svg_lines) + "\n"
        with open(output_path, "w", encoding="utf-8") as svg_file:
            svg_file.write(svg_content)

        return output_path


def _generate_ticks(start: float, end: float, count: int) -> List[float]:
    if count <= 1 or math.isclose(start, end):
        return [start]

    span = end - start
    step = span / (count - 1)
    ticks = [start + step * i for i in range(count)]
    return ticks


def main() -> None:
    random.seed(RANDOM_SEED)
    name_path_input = "1R-20I-150C-2P"
    params = AHADiscreteParams(population_size=50, max_iterations=100)
    aha = AHADiscrete(name_path_input, params)
    best_solution, best_fitness, best_evaluation = aha.run()

    print("\nBest tardiness:", best_fitness)
    if best_evaluation is not None:
        print("Best evaluation detail:", best_evaluation)

    plot_path = aha.plot_best_history()
    print(f"Progress plot saved to {plot_path}")


if __name__ == "__main__":
    main()
