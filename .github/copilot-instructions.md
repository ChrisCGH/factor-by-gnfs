# factor-by-gnfs Copilot Instructions

## Build and test

- The C++ sources and makefiles are in `gnfs/`; run native commands from that directory.
- Prefer the Clang build: `make -f makefile.clang all`. It builds into `gnfs/gbin/` and requires Clang, GMP, CppUnit, libuuid, and the C++ standard library. The Linux configuration uses C++20 and accepts `OPENMP=yes`.
- Run the native suite with `make -f makefile.clang test`.
- Build and run one C++ test with, for example:
  ```sh
  make -f makefile.clang LatticeSieverUnitTest
  ./gbin/LatticeSieverUnitTest
  ```
  Every target in `ALL_TESTS` follows this pattern. The unit-test targets link CppUnit even though the executables are invoked directly.
- Run the coverage target with `make -f makefile.clang coverage`; it rebuilds with coverage instrumentation and requires `gcovr` plus `llvm-cov gcov`.
- Root-level Python tools have standard-library `unittest` suites. Run all of them from the repository root with `python3 -m unittest discover -s . -p 'test_*.py'`, or one module with `python3 -m unittest test_trial_sieve`.

## Architecture

This is a General Number Field Sieve implementation, not a single executable:

1. `skewed` performs polynomial selection from `skewed.cfg`. Its output supplies the number, common root, and polynomials for a `sieve.cfg`.
2. `lsieve` is the lattice-sieving entry point and delegates the relation-generation work to `LatticeSiever`. It reads `sieve.cfg` in its current working directory and appends generated relations to `RELATION_FILE`. The older `sieve`/`Siever` path is a separate, non-lattice sieve implementation.
3. Post-processing converts relations into a solvable sparse matrix: `filter` deduplicates, removes singletons, prunes, and merges relations into relation sets; `buildMatrix` reads the relations, relation sets, relation-set primes, and sieve config to write the matrix; `bl` runs block Lanczos to produce dependencies.
4. `processDependencies` selects relations for a dependency. `processDependencies.pl` coordinates it with `calcroot`, which performs the algebraic and rational square-root step to recover a factor. `postProcessing.pl` is the legacy driver for this pipeline.

Core arithmetic is shared between stages: `VeryLong`/`VeryLongModular` handle arbitrary-precision arithmetic, `Polynomial` and `NumberField` represent GNFS polynomials and algebraic structures, and `SparseMatrix`/`blockLanczos` implement the GF(2) linear-algebra stage.

## Repository-specific conventions

- Keep changes compatible with the makefile's explicit per-executable source lists. Adding a source file to a component requires adding it to that target's `*_CPP_SRCS` list in `makefile.clang`; there is no automatic source discovery.
- `SieveConfig` and the polynomial-selection config parser are intentionally simple line parsers. Use the exact `KEY = value` format with one space after `=`; values are extracted from two characters after `=`, so `KEY=value` is not equivalent. Lines beginning with `#` are comments.
- Both inline polynomials and multiline `Polynomial` blocks are supported. In a multiline block, preserve the `DEGREE = n` header followed by the coefficients expected by `Polynomial::read_polynomial`.
- Validate generated or edited sieve configs with `python3 validate_sieve_cfg.py path/to/sieve.cfg` before using them. `trial_sieve.py`, `tune_sieve.py`, and `validate_sieve_cfg.py` operate from the repository root; `trial_sieve.py` deliberately gives each candidate a separate working directory because `lsieve` resolves `sieve.cfg` relative to its working directory.
- Relation files and intermediate matrix/dependency artifacts are pipeline inputs, not interchangeable text: preserve the relations/relation-sets/relation-set-primes tuple expected by `filter`, `buildMatrix`, `bl`, and `processDependencies`.
