# MadVoro regression tests

THUNDER-driven regression suite for MadVoro geometry and MPI correctness.

## Layout

- `cases/` — one directory per test with `REGRESSION_INFO` and `main.cpp`
- `lib/regression_checks.sh` — validation functions referenced by `CHECK_FUNCTION`
- `lib/voronoi_test_common.hpp` — shared helpers for Voronoi test cases
- `THUNDER/` — git submodule ([THUNDER](https://github.com/maormizrachi/THUNDER))

## Running locally

From the MadVoro repository root:

```bash
./regression_tests/run_all.sh --list-tests
./regression_tests/run_all.sh --test voronoi_mock_mesh_periodic --with-mpi --local
```

When embedded in RICH, the same tests are discovered automatically via
`regression_tests/config.json` as the `madvoro` subproject.
