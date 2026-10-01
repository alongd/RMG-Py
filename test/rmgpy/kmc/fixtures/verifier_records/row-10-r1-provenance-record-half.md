# Row 10 — R1 provenance/fixtures, record half

- Verdict: `FAIL` (also `BLOCKED-BY-GAP` for the unmapped registered elements)
- Verifier-definition SHA-256: `8d43a49a2a39040cb7884c99d2323797a84e9d4da3536de90548e7dfe51084db`
- Copied-literal-block SHA-256: `ab5105e678c52c6f74dd1f1fc293a1081f7076fc46e458b6a679ce0c570cc798`
- Product commit SHA: `f42bb28bbb605fac424f785a0cd95cf2babae98d`
- Run log: `/tmp/i028-verifier-logs/slow-suite.stdout.log`
- Run-log SHA-256: `c8809605b31ce353bb18f93ba86aec35fea55c1509a25b75256999a8fda0c549`
- Stderr log: `/tmp/i028-verifier-logs/slow-suite.stderr.log`
- Stderr-log SHA-256: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`

Mapped record-half checks passed for the real compiled J_para, J_ortho S7, and
J_ortho S6 create/dissociate pairs: full content IDs, artifact/provenance hashes,
formed/deleted single bonds, radical deltas, zero H transfer, exact graph
round-trips, S6/S7 reflection, persistent UUID bindings, reverse handles, and
injected hash/class/inverse/reflection/binding defects.

The row fails because the strict J_para fix exposes that the real compiled
J_para record does not match the archived J_para rate. Its compiled/archived
ratios at 600, 650, 700, 750, and 800 K are respectively `0.0107820852`,
`0.0152551817`, `0.0206536539`, `0.0269841241`, and `0.0342390731`.
`test_real_r1_tables_compile_ortho_for_every_transport_arm` therefore raises
`ValueError: junction provenance 'J_para' does not match an archived junction
rate`. The required slow command reproduced `54 passed, 1 failed, 1 skipped`.

The row is additionally gap-blocked for immutable parent/cage IDs in the
compiled record, the provenance-bound whole-component G8 validator and R1-3,
and the explicitly deferred R1-4 group-count half. Detailed evidence is in
`../i028_registered_verifier_mapping.md`.
