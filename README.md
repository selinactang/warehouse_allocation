# Recovered Warehouse Allocation Project

This folder reconstructs the MATLAB project described in `Warehouse Allocation Problem`
and the EUDL PDF `10.4108/eai.19-5-2023.2334330`.

The original `part2.mat` file is not present in the surviving folder, so this recovery keeps
two things separate:

1. Values and logic that are directly visible in the PDF or Word document.
2. Reconstructed/default data that can be edited in `data/provinces.csv`.

## What Was Recovered

### Part 1: Minimum Number of Warehouses

The MATLAB screenshot in the paper is transcribed into:

- `matlab/part1_min_warehouses.m`
- `python/warehouse_allocation.py`

The Part 1 candidate provinces are the 14 "non" provinces. Important recovery note:
the literal MATLAB screenshot does **not** reproduce the six-province result when solved
as written. It selects the minimum number of binary variables satisfying the visible
rent-threshold constraint, while the paper text reports a different six-province answer.

For historical recovery, both are preserved:

- `solve_part1_from_paper_fragment` in Python and `part1_min_warehouses.m` transcribe the
  visible screenshot logic.
- `PART1_PAPER_RESULT` records the paper-reported answer.

The paper's reported selected warehouses are:

- Zhejiang
- Liaoning
- Henan
- Heilongjiang
- Jiangxi
- Gansu

The independent provinces listed in the paper are:

- Shanghai, Beijing, Fujian, Guizhou, Hunan, Hubei, Jilin, Shandong, Shanxi, Tianjin,
  Sichuan, Chongqing

The bidependent provinces listed in the paper are:

- Anhui, Jiangsu, Guangdong, Guangxi

### Part 2: Minimum Delivery Cost

The recovered model uses the formula from the paper:

```text
min sum_i sum_j X_ij * p_j * ((e_i + e_j) / 2) * d_ij
```

subject to supply and demand constraints:

```text
row sums of X <= supply capacity
column sums of X >= demand
X >= 0
```

The paper later refines the distance index to:

```text
d_ij = 6 + 4 * distance_ij / max(distance)
```

The warehouse cut-down loop was recovered from the Word document images:

- Try removing each currently active warehouse.
- Redistribute that warehouse's capacity equally to the remaining active warehouses.
- Solve the transportation LP for each candidate removal.
- Permanently remove the candidate with the lowest total cost.
- Continue until only 16 warehouses remain.

The paper's final 16 warehouse capitals are recorded in the scripts:

- Beijing, Tianjin, Shanghai, Chongqing, Nanning, Changchun, Jinan, Taiyuan, Chengdu,
  Wuhan, Hefei, Nanjing, Fuzhou, Changsha, Guiyang, Guangzhou.

The visible result image also contains nonzero capacity numbers for several locations;
those are stored as `paper_final_capacity` in `data/provinces.csv` where readable. The
image is not fully machine-readable, so these values are treated as evidence, not as a
complete original `part2.mat` replacement.

## Run the Python Version

From this folder. If `python3` resolves oddly on your machine, use the absolute Anaconda
Python shown by `which python3` in the parent project folder.

```bash
python3 python/warehouse_allocation.py
```

The Python script prints:

- Part 1 selected warehouses
- The recovered paper result
- The reconstructed Part 2 removal path
- Final active warehouses from the reconstructed dataset

It also writes:

- `outputs/part2_solution.csv`
- `outputs/removal_history.csv`

## Run the MATLAB Version

In MATLAB, set the current folder to `recovered_warehouse_project`, then run:

```matlab
addpath('matlab')
part1_min_warehouses
part2_delivery_lp
cut_down_warehouses
```

`intlinprog` and `linprog` require MATLAB Optimization Toolbox.

## How To Improve The Recovery

If you later find the old Excel table or `part2.mat`, replace the values in
`data/provinces.csv` or load the `.mat` file in the MATLAB scripts. The model code is
structured so the historical data can be swapped in without changing the algorithm.
