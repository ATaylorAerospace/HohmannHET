# Contributing to HohmannHET

Thank you for your interest in HohmannHET! This document explains how to
contribute effectively.

## Tri-Language Parity Rule

Every physics function must be implemented identically in Python, C++, and
MATLAB. If you change a formula or add a feature in one language, you must
update all three. Cross-language numerical agreement must hold to `1e-6`.

## Development Setup

### Python
```bash
cd python
pip install -e ".[dev]"
pytest tests/ -v
```

### C++
```bash
cd cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cd build && ctest --output-on-failure
```

### MATLAB
```matlab
cd matlab/tests
results = runtests({'TestDynamics','TestPropulsion','TestOptimization'});
table(results)
```

## Pull Request Checklist

- [ ] All three language implementations are updated (if applicable)
- [ ] Tests pass in all languages
- [ ] New functions include docstrings / Doxygen / help blocks
- [ ] Physical equations are documented in LaTeX notation
- [ ] Constants match the shared reference values in the README

## Reporting Issues

Please include:
- Which language implementation is affected (Python / C++ / MATLAB / all)
- Minimal code to reproduce the issue
- Expected vs. actual output
- Your environment (OS, compiler/Python version, MATLAB release)

## License

By contributing, you agree that your contributions will be licensed under the
Apache License 2.0, consistent with the project license.
