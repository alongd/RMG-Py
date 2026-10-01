# Row 7 — grammar closure

- Verdict: `BLOCKED-BY-GAP`
- Verifier-definition SHA-256: `8d43a49a2a39040cb7884c99d2323797a84e9d4da3536de90548e7dfe51084db`
- Copied-literal-block SHA-256: `ab5105e678c52c6f74dd1f1fc293a1081f7076fc46e458b6a679ce0c570cc798`
- Product commit SHA: `f42bb28bbb605fac424f785a0cd95cf2babae98d`
- Run log: `/tmp/i028-verifier-logs/slow-suite.stdout.log`
- Run-log SHA-256: `c8809605b31ce353bb18f93ba86aec35fea55c1509a25b75256999a8fda0c549`
- Stderr log: `/tmp/i028-verifier-logs/slow-suite.stderr.log`
- Stderr-log SHA-256: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`

The registered literal hashes reproduce exactly, including closure counts
`(3, 12, 18, 12, 3)`, 48 states, 96 transitions, and transition SHA-256
`fce0f94f8b2562b7fcc6ba7d15485e8b1892240b1a087c6b15eeea686e8b13d1`.
This is not a row pass: the product has no counterparts for the registered
alkane roots and mapped sites, `RW00`–`RW03`, embedding enumeration,
`grammar_class`, `accepts_g8`, or the four-level closure BFS. The complete GAP
evidence is in `../i028_registered_verifier_mapping.md`.

The required slow command completed with `54 passed, 1 failed, 1 skipped`; its
one failure is the independent row-10 J_para rate disagreement, not a row-7
execution. No test-side grammar implementation was used as a product stand-in.
