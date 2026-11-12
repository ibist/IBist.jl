# Julia package instructions

This is the native Julia sister package to ibist-R and ibist-py.
Read the shared policy before adding methods or regenerating data:
https://github.com/ibist/ibist-R/blob/main/CROSS_LANGUAGE_DEVELOPMENT.md
In the shared workspace it is ../ibist-R/CROSS_LANGUAGE_DEVELOPMENT.md.
Inspect Git status and preserve unrelated changes.

- Use module/package name `IBist`, snake_case functions and arguments,
  CamelCase types, and four spaces for indentation.
- Keep source and documentation within 80 characters per line. Preserve
  unrelated wrapping; use compact calls and aligned continuation arguments.
- R defines the settled statistical specification. Julia has no R or
  Python runtime dependency. Reuse authoritative data and reference exports;
  never edit generated data manually or translate their cleaning rules.
- Preserve method order and repetitions, table orientation, defaults,
  interval sides, and full stored numerical precision. Document deliberate
  interface differences. Add independent checks alongside R comparisons.
- Use native Julia numerical code except for algorithms designated as
  shared C++ cores. Never maintain a separate handwritten Wang algorithm.
- Run focused affected checks during development. Do not stage, commit,
  register, or publish unless requested.
