# Registered grammar-closure and R1 record-verifier mapping

This table maps the byte-authoritative elements in
`R-010_v8_2026-09-27_volatility-numerics-prereg.md` lines 1009–1371 and
the row-7/row-10 failure clauses restated in v13. `GAP` means the closest
production code does not implement the registered operation; the verifier must
not replace it with the illustrative v8 implementation.

## Row 7 — grammar closure

| Registered element | Product entry point or GAP | Evidence |
|---|---|---|
| `GRAMMAR_TERMINAL_SCHEMA` | **GAP** | `rmgpy/kmc/compiler.py` compiles RMG molecule/reaction graphs into `EventRecord` objects; it has no `C[h,a,r,e]`/`P[o,c]` grammar-terminal representation or R1-exemption labels. `rg -n 'G8|GRAMMAR_TERMINAL_SCHEMA|r1_label|P\[o,c\]' rmgpy/kmc test/rmgpy/kmc` returns no product match. |
| `GRAMMAR_PRODUCTIONS` (`Component`, `Tree`, `Benzene`, `Localisation`, `Boundary`) | **GAP** | The closest code is reaction generation in `EventCompiler.compile`; it enumerates RMG reaction-family outputs, not the registered graph productions or boundary markers. |
| `GRAMMAR_INITIAL_GRAPHS` (literal n-hexane, n-heptane, n-octane graphs and `S00`–`S03`) | **GAP** | `short_ps_molecule_catalogue` is a product catalogue, but does not expose these three literal roots with the four registered mapped sites. No `S00`–`S03` symbols occur in product code. |
| `RW00_attach_methyl` | **GAP** | No compiled product record is identified with this literal rewrite ID or with the registered `S00` embedding. The compiler's RMG-family records are chemically generated events, not this one-shot grammar attachment. |
| `RW01_attach_benzene` and `benzene_fragment` | **GAP** | No compiled product record is identified with this literal rewrite ID, `S01`, or the registered benzene fragment. |
| `RW02_R1_para` and `j_fragment` | **GAP** | The closest code is `R1_JUNCTION_LABELS['J_para']`, `_junction_ops`, and real compiled J_para records. Those records join two RMG reactants at `S9/P9`; they do not attach the registered J literal fragment at `S02` to the registered alkane roots. |
| `RW03_R1_ortho` and `q_fragment` | **GAP** | The closest code is `R1_JUNCTION_LABELS['J_ortho_S7'/'J_ortho_S6']`, `_ortho_junction_proxy`, and `_ortho_junction_records`. Those records compile S6/S7 cage capture; they do not attach the registered Q literal fragment at `S03`. |
| `compiled_rewrite_embeddings` | **GAP** | The compiler has reaction atom maps (`extract_atom_map`, `_canonical_atom_map`) and the executor has graph matching (`SparseState._match`), but neither enumerates embeddings of each registered grammar rewrite into each registered closure state in lexicographic order. |
| `render_grammar_graph` | **GAP** | `apply_record` and `SparseState.apply` execute real compiled records. Neither accepts `(base_graph, four-bit rewrite mask)` nor renders the registered literal fragments, so they are not the same construction. |
| `grammar_class` | **GAP** | Product component ledgers and volatility descriptors do not compute the registered tuple `(C,H,n_ar,critical groups,R1 labels,rewrite mask)` for grammar states. |
| `accepts_g8` | **GAP** | `rmgpy/kmc/volatility.py` has domain checks for volatility inputs, but no G8 recognizer, no boundary-marker validation, no exemption-label scoping, and no `E_GRAMMAR_CLOSURE` path. `rg -n 'accepts_g8|E_GRAMMAR_CLOSURE|G8' rmgpy/kmc test/rmgpy/kmc` returns no product match. |
| canonical graph representation used for closure deduplication | **GAP** | `_canonical_adjacency` canonicalises individual RMG molecules and `canonical_json_bytes` canonicalises artifact JSON, but neither canonicalises the registered coloured graph plus rewrite mask. |
| exhaustive four-level BFS, `seen_states`, and termination bound | **GAP** | `short_ps_molecule_catalogue` has an independently bounded catalogue traversal, but there is no product BFS over the four registered rewrites and no 16-mask termination proof. |
| `GRAMMAR_FIXTURES` hashes | **GAP** | Product artifact/input hashes bind compiler inputs, not the four registered grammar literals; none equals the registered terminal/production/initial/rewrite hashes. |
| `GRAMMAR_CLOSURE` (`depth_counts`, 48 states, depth transitions, `outside_G8=()`) | **GAP** | Without product counterparts for the registered roots, rewrites, embeddings, classes, and G8 recognizer, the production path cannot emit this table. |
| `GRAMMAR_TRANSITIONS` (96 rows; SHA-256 `fce0f94f8b2562b7fcc6ba7d15485e8b1892240b1a087c6b15eeea686e8b13d1`) | **GAP** | Product records have content-addressed IDs but no registered closure-transition table to hash. |
| `graph_colour` and `CLOSED_G8_GRAPHS` | **GAP** | No product colour renderer or complete G8 closure-state export exists. |
| Failure: outside-G8 product / fail-open handling | **GAP** | There is no product G8 call at which an outside-G8 product can be rejected fail-closed. |
| Failure: nontermination | **GAP** | There is no product traversal of the registered four-bit rewrite system whose termination can be checked. |
| Failure: transition/hash/table mismatch | **GAP** | There is no product transition table for comparison with the registered constants. |

**Row-7 premise verdict: `BLOCKED-BY-GAP`.** The real event compiler has no
counterpart for any of the four registered illustrative rewrites on the
registered alkane roots. Re-running the v8 block in a test would test only the
registered illustration, not the product.

## Row 10 — R1 provenance/fixtures, record half

| Registered element | Product entry point or GAP | Evidence |
|---|---|---|
| `EVENT_ARTIFACT` literal | **GAP as a literal; mapped production authority exists** | The literal bytes are an illustrative fixture and are not a production event artifact. Production artifacts are emitted by `EventCompiler.write_artifact`; their filename is SHA-256 of their canonical bytes and is checked by the verifier. |
| `EVENT_ARTIFACT_SHA256` / artifact hash check | `EventCompiler.write_artifact`, `canonical_json_bytes`; verifier checks `artifact_path.stem == sha256(artifact_path.read_bytes()).hexdigest()` | This is the product's content-addressed artifact authority. The registered literal hash remains a fixed oracle and is checked separately, not substituted for the production artifact. |
| `canonical_graph` | `_canonical_adjacency`, `_graph_adjacencies`, and canonical `reactant_graphs`/`product_graphs` fields | The encodings differ, so the verifier compares the registered literal hashes only to the copied literals and uses product adjacency strings for compiled-record graph checks. |
| `graph_sha256` | `EventRecord._compute_event_id` for semantic record content; artifact SHA for complete artifact bytes | Product records do not carry the v8 `motif_sha256`/`parent_sha256` fields. Exact product parent/product graphs are nevertheless content covered through `reactant_graphs` and `product_graphs`. |
| `literal_graph`; `J_H`, `J_BONDS`, `J_PRODUCT_GRAPH`, `J_PARENT_GRAPH`; `Q_H`, `Q_BONDS`, `Q_PRODUCT_GRAPH`, `Q_PARENT_GRAPH` | Fixed verifier oracles only | These are registered literals, not product entry points. Their hashes are copied verbatim and checked independently; they are not passed through the product as stand-ins. |
| `J_MOTIF`, `Q_MOTIF`, J/Q parent hashes | Fixed verifier constants; **GAP** for direct product motif-hash fields | The compiler records full `reactant_graphs` and `product_graphs`, but has no `motif_sha256`, `parent_sha256`, or mapped-motif atom-list field matching the illustrative J/Q records. |
| `FiredR1Record` schema | `EventRecord`; `inventory_class`, `junction_ops`, `reverse_of`, `reactant_graphs`, `product_graphs`, `provenance` | `validate_artifact` reconstructs each real record, validates its full ID, requires artifact provenance equality, and checks reciprocal reverse links and R1 operation fields. |
| `J_RECORD` / J_para fields | Real compiled J_para create/dissociate records from `EventCompiler.compile`; `_junction_ops` | Maps junction kind, external `P9`, attacked `S9`, formed single bond, para localisation, exact inverse rewrite, reverse handle, full graphs, and radical/bond operations. |
| J_ortho S7/S6 records and reflection | `_ortho_junction_records`, `_linked_ortho_pair`, `R1_JUNCTION_LABELS` | Product compiles distinct S7 and S6 create/dissociate pairs. S6 carries reflection `{S6:S7,S7:S6,S8:S10,S10:S8}`; S4 and S9 are absent from the map and therefore fixed. |
| `Q_RECORD` / `quinoid_free` illustrative record | **GAP** | Product inventory has `R1:quinoid_disproportionation` and volatility accepts a caller-supplied `quinoid_free` class, but no compiled record matches the literal Q record's attacker, parent/cage, motif hash, and exact inverse schema. |
| `full_event_id`; `J_ID`, `Q_ID` | `EventRecord._compute_event_id` and `EventRecord.validate` | Product IDs are `evt_` plus the full 64-hex SHA-256. The prefix means their total string length is 68, not the illustrative bare-hex length 64. |
| `TRUSTED_EVENT_RECORDS` | artifact `records` after `validate_artifact` | Consumers build an ID-indexed mapping from validated real records; no caller-supplied illustrative record is trusted. |
| `RUNTIME_BINDINGS` | `SparseState._preflight`, `_derive_r1_bindings`, `AppliedEvent.bindings`, `_open_captures`, `_find_reverse_capture` | Executor bindings tie record labels to persistent atom UUIDs. Reverse execution requires the captured binding and unchanged component-scope hash. |
| persistent atom UUIDs | `AtomRef.uuid`, `SparseState._allocate_uuid`, `assert_uuid_uniqueness` | UUIDs survive graph rewrites and strand moves; duplicate/reused UUIDs are rejected. |
| parent IDs and cage ID stored in the fired record | **GAP** | The compiled `EventRecord` has graph snapshots and reverse linkage but no fields for both parent IDs or a cage ID. The executor captures participant strands/components at runtime, but that is not immutable record content. |
| `graph_adjacency`, `graph_cycles` | RMG molecule graphs and product adjacency serialization; **GAP** for the exact illustrative helpers | Product code can inspect molecular topology, but the compiler record validator does not expose the registered whole-component cycle routine. |
| `inverse_r1` | `junction_ops[*].exact_inverse_rewrite`, reciprocal reverse records, `apply_record`, `SparseState.apply`, and `apply_compiled_graph_rewrite` | The real product inverse is executed and compared to real compiled parents; it is not reconstructed from the v8 helper. |
| `r1_source_groups` | **GAP / out of scope** | No source-exact `{4,9,10,88,89}` group counter is in the event compiler. `rmgpy/kmc/volatility.py` consumes group data but does not realise the registered source aggregation. |
| `add_cyclohexane` | **GAP** | The product has no whole-component G8 validation tied to event provenance, so an unrelated aliphatic ring cannot be rejected through the row-10 record path. |
| `check_r1_fixture` artifact/ID/record/provenance stages | `validate_artifact`, `EventRecord.validate`, and executor preflight/binding checks | Mapped stages fail closed on artifact, full-ID, provenance, reverse-link, operation, graph, UUID, and binding corruption. |
| `check_r1_fixture` motif/inverse stages | Real compiled graphs plus `apply_record`/reverse record | The verifier checks each real R1 pair's bond/radical/H delta and exact restoration. Direct equality to the illustrative J/Q motif hashes is not claimed. |
| `check_r1_fixture` whole-component G8/group stages | **GAP** | No product call combines trusted event provenance with whole-component G8 and source group counting. |
| R1-1 metadata/class swap → `E_R1_PROVENANCE` | `validate_artifact` rejects semantic mutation because `event_id` no longer matches content | Product raises `ValueError` rather than returning the registered symbolic code; the failure class is mapped to provenance rejection by the verifier. |
| R1-2 external `P9` and mapped attacked atom → accept | `validate_artifact` R1 field checks on real compiled records | Real J_para, J_ortho S7, and J_ortho S6 records must carry the exact attacker/attacked labels and formed bond. |
| R1-3 extra cyclohexane → `E_DOMAIN_ALIPHATIC_RING` | **GAP** | `volatility.py` can reject an aliphatic ring in a separate volatility-domain path, but it is not bound to the immutable fired event record or executor binding required by row 10. |
| R1-4 external conjugation and exact groups | **GAP / out of scope** | Group-count half belongs to the later evaporation build, and no compiler-side implementation exists. |
| no-H-transfer delta | `implicit_h_delta`, `formula_delta`, and bond/feature operations on each real record | The verifier requires zero implicit-H delta, zero H formula delta, and no hydrogen-setting operation. |
| formed/deleted bond and radical changes | `bond_ops`, `junction_ops.formed_bond`, `exact_inverse_rewrite`, `radical_delta` | The verifier compares the real create and dissociate records to the v13 table and executes them against their recorded graphs. |
| full 256-bit content ID | `EventRecord._compute_event_id`, `validate`; artifact `validate_artifact` | Semantic corruption changes the computed ID; truncated IDs are invalid. A control truncates a real ID and requires rejection. |
| provenance binding | record `provenance` equality to artifact provenance; `SparseState` runtime UUID capture | Compile-time provenance and runtime component binding both exist, but immutable parent/cage IDs remain the gap noted above. |
| Failure: hash mismatch | `EventRecord.validate`, artifact filename/content hash check | Injected semantic and artifact-byte defects must fail. |
| Failure: binding/UUID mismatch | `SparseState._find_reverse_capture`, `_binding_scope_hash`, `assert_uuid_uniqueness` | Existing executor tests cover stale/incorrect UUID captures; the row verifier records this mapped product path. |
| Failure: provenance/class/reverse mismatch | `validate_artifact` | Injected field, class, reverse-link, and inverse defects must fail. |
| Failure: inverse mismatch | `validate_artifact` plus execution of real forward/reverse records | Wrong inverse operation/bond order is rejected or fails restoration. |
| Failure: G8/valence/conjugation mismatch | **GAP** | Compiler validation checks record structure, not whole-component G8, valence, or unrelated conjugation under the registered authority rule. |
| Failure: fail-open declared J_para rate | `_channel_id` in `rmgpy/kmc/met.py`; confirmed product defect | Before this work, unmatched declared J_para fell through as J_para while J_ortho raised. The only allowed production change makes both declared junction kinds strict and adds an injected wrong-rate control. The real compiled J_para record then also rejects: over the registered 600–800 K grid its compiled/archived rate ratio is `0.0107821, 0.0152552, 0.0206537, 0.0269841, 0.0342391`, so the dispatch premise that every real compiled J_para still matches is false. |

**Row-10 record-half verdict:** `FAIL` because the real compiled J_para rate does
not match the archived J_para channel once fail-open handling is removed. The row
is additionally `BLOCKED-BY-GAP` for immutable parent/cage record fields,
whole-component G8 validation, and R1-3. R1-4/group counting is explicitly the
later evaporation half and is also recorded as a gap here.

## Search evidence

The premise audit used only the named product worktree and the registered
reports. The decisive searches were:

```text
rg -n 'G8|accepts_g8|GRAMMAR|RW00_attach|compiled_rewrite_embeddings|grammar_class|E_GRAMMAR_CLOSURE' rmgpy/kmc test/rmgpy/kmc
    no product implementation match

rg -n 'J_para|J_ortho|reverse_of|persistent_uuid|artifact|sha256|canonical|_channel_id|compiled' rmgpy/kmc test/rmgpy/kmc
    maps to compiler.py, event_record.py, state.py, met.py and existing product tests

rg -n 'E_R1_PROVENANCE|E_DOMAIN_ALIPHATIC_RING|source_groups|motif_sha256|parent_sha256|cage' rmgpy/kmc test/rmgpy/kmc
    no compiler-side registered authority path; only the separate volatility-domain ring rejection exists
```
