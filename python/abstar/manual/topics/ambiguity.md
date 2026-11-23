How abstar reports gene calls when several alleles tie. The full record is decision 0008 (`abstar manual decisions 0008`).

## One primary call

For each of V, D, J and C (shown for V):

| Field | Contents |
|---|---|
| `v_call` | Exactly one allele, always. Searchable, countable, sortable. |
| `v_gene` | The gene of `v_call`. |
| `v_species` | The species of `v_call`. |
| `v_call_ties` | Every allele with the top score, including `v_call`. A native list in Parquet, comma-joined in TSV. |
| `v_species_ties` | Aligned entry by entry with `v_call_ties`. |

- Tied alleles share the top score, so `v_score` covers them all. There is no score list.
- There is no ambiguity flag. A call is ambiguous when `len(v_call_ties) > 1`.
- The name `v_calls` is deliberately avoided: it is one letter from `v_call`.

## Which alleles are searched

Alleles the database marks non-functional (IMGT pseudogene `P`, OGRDB `NF`) are not searched, so they never appear in calls or ties. `F`, `ORF` and `unknown` alleles always are. `--include-nonfunctional` (`abstar.run(..., include_nonfunctional=True)`) searches them too, e.g. for genomic-DNA repertoires or genotyping (decision 0038).

## Choosing the primary call

`v_call` is picked among the ties by these tie-breakers, in order:

1. Best alignment score.
2. Functional before ORF before pseudogene.
3. Population usage frequency, if the database provides an optional `usage` column in `alleles.tsv`.
4. A fixed natural sort order on allele name.

The same input and database always give the same call.

## Consequences for analysis

- `v_call` is always usable directly.
- Allele counts within tie sets lean toward common, functional alleles. The bias is small and runs in the least harmful direction.
- Rigorous analyses should use `v_call_ties`.
- Evaluation: a call is correct if `v_call` is in the generator's set of equally valid answers. Coverage of that set by `v_call_ties` is scored separately.
