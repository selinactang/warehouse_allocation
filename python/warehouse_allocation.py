"""Recovered warehouse allocation project.

This script reconstructs the two optimization models described in:

    Selina Chen Tang and Changkun Zuo, "Warehouse Allocation Problem",
    DOI: 10.4108/eai.19-5-2023.2334330

The old MATLAB `part2.mat` data file is not available in the surviving project folder.
The algorithm is recovered from the paper and document screenshots, while numerical
inputs live in `data/provinces.csv` so they can be replaced if the original data is found.
"""

from __future__ import annotations

import csv
import heapq
import math
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_FILE = ROOT / "data" / "provinces.csv"
OUT_DIR = ROOT / "outputs"


PART1_PAPER_RESULT = [
    "Zhejiang",
    "Liaoning",
    "Henan",
    "Heilongjiang",
    "Jiangxi",
    "Gansu",
]

PART2_PAPER_FINAL_CAPITALS = [
    "Beijing",
    "Tianjin",
    "Shanghai",
    "Chongqing",
    "Nanning",
    "Changchun",
    "Jinan",
    "Taiyuan",
    "Chengdu",
    "Wuhan",
    "Hefei",
    "Nanjing",
    "Fuzhou",
    "Changsha",
    "Guiyang",
    "Guangzhou",
]


@dataclass(frozen=True)
class Province:
    province: str
    capital: str
    group: str
    part1_candidate: bool
    rent_value_from_paper: float | None
    population_index: float
    economic_index: float
    supply_capacity: float
    demand: float
    paper_final_capacity: float
    latitude: float
    longitude: float


def read_provinces(path: Path = DATA_FILE) -> list[Province]:
    rows: list[Province] = []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rent = row["rent_value_from_paper"].strip()
            rows.append(
                Province(
                    province=row["province"],
                    capital=row["capital"],
                    group=row["group"],
                    part1_candidate=row["part1_candidate"] == "1",
                    rent_value_from_paper=float(rent) if rent else None,
                    population_index=float(row["population_index"]),
                    economic_index=float(row["economic_index"]),
                    supply_capacity=float(row["supply_capacity"]),
                    demand=float(row["demand"]),
                    paper_final_capacity=float(row["paper_final_capacity"]),
                    latitude=float(row["latitude"]),
                    longitude=float(row["longitude"]),
                )
            )
    return rows


def solve_part1_from_paper_fragment(provinces: list[Province]) -> list[str]:
    """Recover the integer program shown in the paper's MATLAB screenshot.

    The original code minimizes the number of selected candidate provinces:

        f = ones(14,1)
        intcon = 1:14
        b = -386.9 * ones(1,14)
        A_row = -[31.66, 29.5, ..., 150]
        A = repmat(A_row, 14, 1)
        lb = zeros(14, 1)
        ub = ones(14, 1)
        intlinprog(f, intcon, A, b, [], [], lb, ub)

    Because every row of A is identical, this is equivalent to:

        minimize sum(x)
        subject to rent_values @ x >= 386.9
        x in {0, 1}

    With the candidate order below, this reproduces the paper's six selected provinces.
    """

    candidates = [p for p in provinces if p.part1_candidate]
    rents = [float(p.rent_value_from_paper) for p in candidates]

    # With only 14 binary variables, exhaustive search is simpler and more
    # reproducible than depending on a particular MILP package.
    best_mask: int | None = None
    best_count = len(candidates) + 1
    best_total = float("inf")
    for mask in range(1 << len(candidates)):
        selected_count = mask.bit_count()
        if selected_count > best_count:
            continue
        total = sum(rent for i, rent in enumerate(rents) if mask & (1 << i))
        if total >= 386.9 and (selected_count < best_count or total < best_total):
            best_mask = mask
            best_count = selected_count
            best_total = total

    if best_mask is None:
        raise RuntimeError("Part 1 search did not find a feasible solution.")

    return [p.province for i, p in enumerate(candidates) if best_mask & (1 << i)]


def haversine_km(a: Province, b: Province) -> float:
    radius = 6371.0
    lat1 = math.radians(a.latitude)
    lat2 = math.radians(b.latitude)
    dlat = lat2 - lat1
    dlon = math.radians(b.longitude - a.longitude)
    h = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * radius * math.asin(math.sqrt(h))


def distance_index(provinces: list[Province], low: float = 6.0, width: float = 4.0) -> list[list[float]]:
    n = len(provinces)
    dist = [[0.0 for _ in range(n)] for _ in range(n)]
    for i, origin in enumerate(provinces):
        for j, dest in enumerate(provinces):
            dist[i][j] = haversine_km(origin, dest)
    max_dist = max(max(row) for row in dist)
    return [[low + width * value / max_dist for value in row] for row in dist]


class MinCostFlow:
    def __init__(self, node_count: int) -> None:
        self.graph: list[list[list[float]]] = [[] for _ in range(node_count)]

    def add_edge(self, start: int, end: int, capacity: float, cost: float) -> None:
        forward = [end, len(self.graph[end]), capacity, cost, 0.0]
        backward = [start, len(self.graph[start]), 0.0, -cost, 0.0]
        self.graph[start].append(forward)
        self.graph[end].append(backward)

    def send(self, source: int, sink: int, required_flow: float) -> tuple[float, list[list[list[float]]]]:
        n = len(self.graph)
        potential = [0.0] * n
        total_cost = 0.0
        flow_sent = 0.0
        eps = 1e-9

        while flow_sent + eps < required_flow:
            dist = [float("inf")] * n
            parent_node = [-1] * n
            parent_edge = [-1] * n
            dist[source] = 0.0
            heap: list[tuple[float, int]] = [(0.0, source)]

            while heap:
                current_dist, u = heapq.heappop(heap)
                if current_dist > dist[u] + eps:
                    continue
                for edge_index, edge in enumerate(self.graph[u]):
                    v, _, capacity, cost, _ = edge
                    if capacity <= eps:
                        continue
                    reduced = cost + potential[u] - potential[int(v)]
                    candidate = current_dist + reduced
                    if candidate + eps < dist[int(v)]:
                        dist[int(v)] = candidate
                        parent_node[int(v)] = u
                        parent_edge[int(v)] = edge_index
                        heapq.heappush(heap, (candidate, int(v)))

            if math.isinf(dist[sink]):
                raise RuntimeError("Transportation LP is infeasible with the current supply and demand.")

            for i in range(n):
                if not math.isinf(dist[i]):
                    potential[i] += dist[i]

            add = required_flow - flow_sent
            v = sink
            while v != source:
                u = parent_node[v]
                edge = self.graph[u][parent_edge[v]]
                add = min(add, edge[2])
                v = u

            v = sink
            while v != source:
                u = parent_node[v]
                edge_index = parent_edge[v]
                edge = self.graph[u][edge_index]
                reverse_index = int(edge[1])
                edge[2] -= add
                edge[4] += add
                reverse = self.graph[v][reverse_index]
                reverse[2] += add
                reverse[4] -= add
                total_cost += add * edge[3]
                v = u

            flow_sent += add

        return total_cost, self.graph


def solve_transport_lp(
    provinces: list[Province],
    supply: list[float],
    demand: list[float],
    use_distance: bool = True,
) -> tuple[float, list[list[float]]]:
    """Solve the Part 2 transportation linear program."""

    n = len(provinces)
    p = [province.population_index for province in provinces]
    e = [province.economic_index for province in provinces]
    d = distance_index(provinces) if use_distance else [[1.0 for _ in range(n)] for _ in range(n)]

    cost = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cost[i][j] = p[j] * ((e[i] + e[j]) / 2.0) * d[i][j]

    total_demand = sum(demand)
    if sum(supply) + 1e-9 < total_demand:
        raise RuntimeError("Total supply is smaller than total demand.")

    source = 0
    origin_offset = 1
    destination_offset = origin_offset + n
    sink = destination_offset + n
    network = MinCostFlow(sink + 1)

    for i, amount in enumerate(supply):
        network.add_edge(source, origin_offset + i, amount, 0.0)
    origin_destination_edges: list[tuple[int, int, int]] = []
    for i in range(n):
        for j in range(n):
            edge_index = len(network.graph[origin_offset + i])
            network.add_edge(origin_offset + i, destination_offset + j, total_demand, cost[i][j])
            origin_destination_edges.append((i, j, edge_index))
    for j, amount in enumerate(demand):
        network.add_edge(destination_offset + j, sink, amount, 0.0)

    total_cost, graph = network.send(source, sink, total_demand)

    solution = [[0.0 for _ in range(n)] for _ in range(n)]
    for i, j, edge_index in origin_destination_edges:
        solution[i][j] = graph[origin_offset + i][edge_index][4]

    return total_cost, solution


def cut_down_warehouses(
    provinces: list[Province],
    minimum_active: int = 16,
) -> tuple[list[float], list[dict[str, object]]]:
    """Recovered version of the old MATLAB cut-down loop.

    This mirrors the screenshot:

        b_old = b;
        b0 = [];
        m0 = [];
        for z = 1:15
            total_costs = zeros(31,1);
            for i = 1:31
                b = b_old;
                b_temp = b(i,1)/(30-length(b0));
                ...
                weight = linprog(...);
            end
            [m,k] = min(total_costs);
            b_temp = b_old(k,1)/(30-length(b0));
            b_old(k,1) = 0;
            ...
            b0 = [b0, k];
            m0 = [m0, m];
        end

    The recovered Python version only tests currently active warehouses, which avoids
    repeatedly selecting an already-removed zero-capacity province.
    """

    supply = [p.supply_capacity for p in provinces]
    demand = [p.demand for p in provinces]
    if sum(supply) < sum(demand):
        scale = sum(demand) / sum(supply) * 1.02
        supply = [amount * scale for amount in supply]

    active = [amount > 0 for amount in supply]
    history: list[dict[str, object]] = []

    while sum(active) > minimum_active:
        best_cost = float("inf")
        best_index: int | None = None
        best_supply: list[float] | None = None

        for i, is_active in enumerate(active):
            if not is_active:
                continue
            candidate_supply = supply.copy()
            removed_amount = candidate_supply[i]
            candidate_supply[i] = 0.0
            receivers = [j for j, active_j in enumerate(active) if active_j and j != i]
            for j in receivers:
                candidate_supply[j] += removed_amount / len(receivers)

            try:
                total_cost, _ = solve_transport_lp(provinces, candidate_supply, demand)
            except RuntimeError:
                continue

            if total_cost < best_cost:
                best_cost = total_cost
                best_index = int(i)
                best_supply = candidate_supply

        if best_index is None or best_supply is None:
            raise RuntimeError("No feasible warehouse removal was found.")

        supply = best_supply
        active[best_index] = False
        history.append(
            {
                "step": len(history) + 1,
                "removed_index_1_based": best_index + 1,
                "removed_province": provinces[best_index].province,
                "removed_capital": provinces[best_index].capital,
                "total_cost": best_cost,
                "active_count": sum(active),
            }
        )

    return supply, history


def write_outputs(provinces: list[Province], solution: list[list[float]], history: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(exist_ok=True)

    with (OUT_DIR / "part2_solution.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["origin"] + [p.capital for p in provinces])
        for origin, row in zip(provinces, solution):
            writer.writerow([origin.capital] + [round(float(x), 6) for x in row])

    with (OUT_DIR / "removal_history.csv").open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["step", "removed_index_1_based", "removed_province", "removed_capital", "total_cost", "active_count"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(history)


def main() -> None:
    provinces = read_provinces()

    print("Part 1 recovered integer-program result")
    part1_selected = solve_part1_from_paper_fragment(provinces)
    print("  Literal screenshot-code solver selected:", ", ".join(part1_selected))
    print("  Paper text reports:                  ", ", ".join(PART1_PAPER_RESULT))
    if set(part1_selected) != set(PART1_PAPER_RESULT):
        print("  Note: these differ; both are preserved for historical recovery.")
    print()

    supply = [p.supply_capacity for p in provinces]
    demand = [p.demand for p in provinces]
    if sum(supply) < sum(demand):
        scale = sum(demand) / sum(supply) * 1.02
        supply = [amount * scale for amount in supply]
    base_cost, base_solution = solve_transport_lp(provinces, supply, demand)
    print("Part 2 recovered transportation LP")
    print(f"  Baseline active warehouses: {sum(amount > 0 for amount in supply)}")
    print(f"  Baseline reconstructed cost: {base_cost:,.2f}")
    print()

    final_supply, history = cut_down_warehouses(provinces, minimum_active=16)
    final_cost, final_solution = solve_transport_lp(provinces, final_supply, demand)
    active_final = [p.capital for p, amount in zip(provinces, final_supply) if amount > 1e-9]
    print("Cut-down reconstruction")
    print(f"  Removed warehouses: {len(history)}")
    print(f"  Final active warehouses: {len(active_final)}")
    print("  Final capitals:", ", ".join(active_final))
    print(f"  Final reconstructed cost: {final_cost:,.2f}")
    print()

    paper_caps = [p.capital for p in provinces if p.paper_final_capacity > 0]
    print("Final 16 warehouse capitals reported/marked in the paper")
    print("  " + ", ".join(PART2_PAPER_FINAL_CAPITALS))
    print()

    print("Nonzero capacity numbers visible in final paper image")
    print("  " + ", ".join(paper_caps))
    print()

    write_outputs(provinces, final_solution, history)
    print(f"Wrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
