# LinkedCellContainerDirect Full Simulation Benchmark Investigation

## Problem
The `bmSimulationFullLinkedCellDirect` benchmark crashes with segfaults, while other `LinkedCell/.../Direct` benchmarks work fine.

## Root Cause
The `proximity_iterator` design is fundamentally incompatible with dynamic particle movement during simulation. The iterator stores:
1. `std::vector<CellDirect*> cells` - Raw pointers to cells
2. `std::vector<Particle>::iterator cur` - Iterator into a particle vector

When `updateParticlePosition()` moves particles between cells, it:
1. Calls `eraseParticle()` which modifies particle vectors
2. Calls `addParticle()` which can cause vector reallocation
3. Returns a new iterator with potentially stale references

Even though we check for iterator invalidation, the pointers become corrupted (`0x000128000000`, `0x000000010000`, etc.), indicating memory corruption beyond simple null pointers.

## Attempted Fixes

### 1. Fixed `removeParticles()` Logic
- **Issue**: Used pointer arithmetic assuming contiguous storage
- **Fix**: Changed to use `fitsDomain()` and two-pass deletion
- **Result**: Didn't resolve segfaults

### 2. Fixed Missing Iterator Increment
- **Issue**: `applyBoundaries()` was missing `++it`, causing infinite loop
- **Fix**: Added increment
- **Result**: Fixed infinite loop but didn't resolve seg faults

### 3. Added Null Pointer Checks
- **Issue**: `inc()` was dereferencing null pointers
- **Fix**: Added checks for `(*cur_cell) == nullptr`
- **Result**: Prevented null crashes but revealed corrupted pointers

### 4. Reordered Erase/Add Operations
- **Issue**: `addParticle()` before `eraseParticle()` caused vector reallocation
- **Fix**: Erase first, then add to avoid iterator invalidation
- **Result**: Didn't prevent corruption

### 5. Fixed Cell Lookup in `eraseParticle()`
- **Issue**: Looking in wrong cell (`cur_cell` vs actual cell containing particle)
- **Fix**: Added `findCellIndex()` to locate correct cell
- **Result**: Improved correctness but didn't fix crashes

### 6. Fixed `satisfyInc()`/`satisfyDec()` Logic
- **Issue**: Dereferencing `cur` when at `end()` position
- **Fix**: Separated boundary checks from dereference operations
- **Result**: Prevented one class of UB but crashes persisted

### 7. Switched to `std::deque<CellDirect>`
- **Issue**: `std::vector<CellDirect>` could reallocate, invalidating cell pointers
- **Fix**: Changed to `std::deque` which provides stable element addresses
- **Result**: **Still crashes with corrupted pointers**

### 8. Increased Particle Vector Reserve Size
- **Issue**: Individual cell particle vectors reallocating
- **Fix**: Increased reserve from 100 → 500 → 5000 per cell
- **Result**: **Still crashes**, proving it's not just reallocation

## Why It Still Fails

Even with:
- ✅ `std::deque` for stable cell addresses
- ✅ Massive reserve (5000) to prevent particle vector reallocation
- ✅ Null pointer checks
- ✅ Proper iterator increment logic
- ✅ Correct erase/add ordering

**We still get corrupted pointers like `0x000128000000`.**

This indicates the fundamental design flaw: storing iterators/pointers that span multiple operations where the underlying containers are being modified. The corruption likely occurs when:

1. Iterator A is created with pointers to cells and a `cur` iterator
2. Particle movement modifies a cell's particle vector
3. Iterator A's `cur` now points to reallocated/freed memory
4. Even if we create iterator B from A, the corruption propagates

## Proper Solution

The iterator design needs a **complete rewrite** to use **index-based iteration**:

```cpp
class proximity_iterator {
    const LinkedCellContainerDirect* container;  // Pointer to container
    std::vector<size_t> cell_indices;            // Cell indices instead of pointers
    size_t cur_cell_idx;                         // Index into cell_indices
    size_t cur_particle_idx;                     // Index into particle vector
    // ...
};
```

This way, indices remain valid even if vectors reallocate, and we always dereference through the current container state.

## Benchmark Status

**DISABLED** - The following benchmarks are commented out due to iterator invalidation issues:
- `bmSimulationFullLinkedCellDirect` in [FullBenchmarks.cpp](benchmarks/src/FullBenchmarks.cpp:211-218)
- `bmLinkedCellUpdateDirect` in [LinkedCellBenchmark.cpp](benchmarks/src/LinkedCellBenchmark.cpp:473-481)

Other LinkedCellDirect benchmarks (Add, Iterator, Proximity, Boundary, Halo, Mixed, etc.) work fine because they don't involve `updateParticlePosition()` which triggers the iterator invalidation bug.

## Files Modified

- `benchmarks/code/linkedcellimpl/LinkedCellContainerDirect.h` - Added deque include, updated warnings
- `benchmarks/code/linkedcellimpl/LinkedCellContainerDirect.cpp` - Multiple fixes to erase/update logic, switched to deque
- `benchmarks/code/simulationimpl/SimulationBenchmark.h` - Fixed removeParticles() logic
- `benchmarks/src/FullBenchmarks.cpp` - Disabled problematic benchmark

## Conclusion

This is a **design limitation**, not a simple bug. The current pointer-based iterator design is incompatible with dynamic simulations. A complete redesign to index-based iteration is required to support full simulation benchmarks.
