# Provenance

`modules/provenance.py` is the list of evidence statuses stored on scientific
results. Engine detection uses a separate, lowercase vocabulary in
`modules/engine_validation.py`.

## Evidence statuses

| Status | Meaning |
| --- | --- |
| `RETRIEVED` | The value came from a remote record or file that HelixScope did not calculate |
| `COMPUTED` | The value was calculated here from the input and the named method |
| `PREDICTED` | A model produced the value. RNA folds and AlphaFold files are predictions |
| `HEURISTIC` | A simplified local score. The CRISPR positional score is in this class |
| `EXPERIMENTAL` | Experimental coordinates or an experimental structure record |
| `ILLUSTRATIVE` | A schematic drawing. Not a measured structure |
| `UNAVAILABLE` | The method or source did not produce a value. This is not zero |
| `ERROR` | The operation failed. The failure is not replaced with a placeholder result |
| `PARTIAL` | Only part of the requested object was obtained. The gap stays visible |
| `UNMAPPED` | No corresponding feature was found. Absence is not a negative scientific claim |
| `STALE` | The displayed object no longer matches the current input |
| `RESOURCE_LIMIT` | A declared size, time or memory cap stopped the operation |

`NaN` is used when a ratio has no canonical bases. It is not coerced to zero.

## Engine statuses

| Status | Meaning |
| --- | --- |
| `not_installed` | The executable or module is absent. HelixScope does not imitate it |
| `detected` | A candidate binary was found. It has not been shown to run |
| `live_validated` | A local execution check succeeded for that tool |
| `remote_validated` | A remote service answered a declared check |
| `broken` | A check ran and failed |
| `unavailable` | No usable status could be established |
| `version_unavailable` | The file is present and the version string was not obtained |
| `invalid_executable` | A candidate existed and was rejected |

`detected` is not `live_validated`.

## Other labels you will see

These are not members of `EVIDENCE_STATUSES`, but they appear on specific
objects:

- `NO_HITS`: a BLAST search finished and reported no hits.
- `READY`: a local BLAST database or reference is ready to search. Ready is
  not a completed genome-wide search.
- `UNCERTAIN`: a structure-mapping label when the chain mapping is not settled.
- `TEST_ONLY`: an object that exists for a test or a test reference, not as a
  public genome.

## What a record keeps

Where the module supports it, the result includes the method name, the
software version `0.24.3-19`, and the time of a remote retrieval. Sequence
history in the browser is session state. It is not a laboratory notebook on
disk.
