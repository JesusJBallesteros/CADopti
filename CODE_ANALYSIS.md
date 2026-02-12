# CADopti Code Analysis: Compute and Memory Optimization Opportunities

## Scope reviewed
- `CADopti.m`
- `TestPair_ref.m`
- `FindAssemblies_recursive_prepruned.m`
- `assemblies_across_bins.m`
- `assembly_assignment_matrix.m`

## High-impact findings

1. **Pairwise screening is the dominant runtime cost (`O(n_neurons^2 * n_bins * n_lags)`).**
   The nested loops over `(w1, w2)` and over all `BinSizes` in `CADopti.m` dominate elapsed time.

2. **`TestPair_ref.m` repeatedly allocates large temporary cell arrays and matrices.**
   In particular, per-call allocation of `Zaa`, `chunked`, `MargPr_t`, `Mx0`, `covABAB`, and `covABBA` is expensive for both CPU and memory.

3. **Memory growth is amplified by dynamic array/cell expansion.**
   Multiple places append inside loops (e.g., `p_values=[p_values; ...]`, cell array growth for assemblies), causing repeated reallocations.

4. **Some user-facing parameters are currently not enforced in core computation.**
   `bytelimit` and `No_th` are stored in output parameters but not actively used to gate agglomeration/pruning during search.

## Optimization recommendations (no change in algorithmic intent)

### A) Runtime improvements

- **Precompute and reuse per-bin metadata once.**
  In `CADopti.m`, compute reusable values like `(2*MaxLags(gg)+1)` and indexed bin lookups outside inner loops.

- **Replace repeated concatenation with preallocation + indexing.**
  Preallocate likely upper bounds for `p_values` and assembly containers, and track a write pointer.

- **Avoid repeated `find(BinSizes==value)` in hot loops.**
  Store the bin index directly with each assembly (e.g., `bin_idx`) to avoid linear search each time.

- **Convert some `parfor` payload outputs to reduction-safe containers.**
  Current logic appends to shared arrays inside `parfor`; refactoring to per-worker local buffers merged after loop can reduce overhead and improve compatibility.

### B) Memory efficiency improvements

- **Use narrower numeric types where mathematically safe.**
  `lag`, indices, and binary masks can often be stored as `int16/uint16/logical` rather than default `double`.

- **Avoid storing full intermediate cell-of-cell covariance structures when only sums are needed.**
  In `TestPair_ref.m`, if only accumulated covariance totals are used, accumulate directly into matrices (`varT`, `covX`) without retaining every `covABAB{i,j}`/`covABBA{i,j}` entry.

- **Minimize duplicate assembly copies.**
  During pruning, keep references/indices where possible rather than copying entire structures repeatedly.

- **Use sparse representation for assignment matrices when appropriate.**
  In `assembly_assignment_matrix.m`, `AAT` can be sparse while constructing/ordering for large populations.

### C) Algorithm-preserving structural improvements

- **Early reject candidate pairs before full test.**
  Cheap heuristics (e.g., low co-activity upper bounds based on marginals) can skip expensive variance/covariance computation when significance is impossible.

- **Apply `No_th` and `bytelimit` during agglomeration, not only as metadata.**
  Enforce thresholds while expanding higher-order assemblies to avoid unnecessary branch growth.

- **Cache transformed spike trains for repeated unit usage.**
  Units are tested many times against different assembly seeds; caching reusable transforms can reduce repeated work.

## Suggested implementation order

1. Preallocation and removal of dynamic growth (`p_values`, assembly buffers).
2. Eliminate unnecessary covariance cell storage in `TestPair_ref.m`.
3. Add active enforcement of `No_th` and `bytelimit` during candidate expansion.
4. Add optional sparse mode for visualization matrix construction.

## Validation strategy after optimization

- Verify exact assembly outputs (`elements`, `lag`, `pr`) against baseline on `Data.mat`.
- Profile with MATLAB `profile on/off` before and after each step.
- Track peak memory usage for representative neuron counts and bin configurations.

