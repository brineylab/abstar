# Transparent human TCR fixtures

These four repository-owned synthetic sequences contain one TRA, TRB, TRD and
TRG construction, using the committed human TCR database. They are assembly
controls, not patient reads or a comprehensive TCR validation set. No external
corpus, user database, annotation output or network service supplies their
expected biological answers.

The database manifest records IMGT GENE-DB, Homo sapiens, downloaded March 27,
2025, functional TR sequences (with the manifest's F+ORF+in-frame P display
option). All sequence and file hashes below use SHA-256. Allele hashes cover
uppercase FASTA sequence characters without whitespace, retaining dots for
IMGT-gapped records; source file hashes cover exact committed bytes. Raw allele
IDs retain species suffixes; public calls remove `__homo_sapiens`.

## Construction and biological evidence

Coordinates in JSON are zero-based, half-open. `germline_start/end` address a
named ungapped source; `query_start/end` address the assembled forward query.
Each V retains its prefix through the last complete codon. The listed synthetic
N payload follows V. TRB adds full TRBD1*01 and TRD adds full TRDD2*01, followed
by a literal A to join the J coding frame. TRA/TRG have no D segment. Each J
omits its last single nucleotide, then a 120-nt C slice `[3:123]` is appended.
This deliberately synthetic C splice tests paired V(D)J/C serialization; it
makes no claim about a natural splice or complete constant-region transcript.
All retained germline bases are unchanged. C slices avoid the leading ambiguous
base of packaged TRAC/TRDC, and all complete query codons are stop-free.

The conserved V Cys is the codon at gapped IMGT position 104 (`[309:312]`).
Its ungapped offset anchors the junction. The independently identified J
F-G-X-G motif supplies the terminal Phe codon. Translating the constructed
junction and full coding sequence establishes productive=True with no
productivity issues. JSON records the full nucleotide junction, its translation,
anchor offsets, coding translation and every concatenation boundary/payload.

V/J alignment evidence uses the longest exact reference match at the known
construction origin. An adjacent payload/C base identical to the omitted
terminal germline base cannot be distinguished by alignment: e.g. TRA V extends
one nucleotide into N1; TRA/TRB/TRD J extends one nucleotide into the C slice.
These sequence-supported boundaries are recorded separately from construction
boundaries. Junction coordinates are currently internal; public tests assert
the exact junction and its oriented sequence interval without adding fields to
the public schema.

Allowed call values are exact sorted comma-separated allele sets, calculated
from identity across each retained reference slice. They are not lists of
interchangeable preferred calls. TRBC1*01/*02 and TRGC1*01/*03/*04/*05 are
indistinguishable in the retained C slices and must remain tied.

| Locus | V | D | J | C | N1 / N2 | Junction AA | Sequence SHA-256 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TRA | TRAV1-1*01__homo_sapiens | null | TRAJ10*01__homo_sapiens | TRAC*01 | GCC / none | CAVRAILTGGGNKLTF | e5117e0feaf6be9878065256b7556ed59b35906c2ffcee9c5fa6bcf5e7657d77 |
| TRB | TRBV10-1*01__homo_sapiens | TRBD1*01__homo_sapiens | TRBJ1-1*01__homo_sapiens | TRBC1*01 | GCC / A | CASSEAGTGGMNTEAFF | 7e16e68bbb5b929f92b04850d2e27686d4df4caccebb3133f1e25aea0486bbdf |
| TRD | TRDV1*01__homo_sapiens | TRDD2*01__homo_sapiens | TRDJ1*01__homo_sapiens | TRDC*01 | GCC / A | CALGEAPSYNTDKLIF | 147f63494f6e14d6ef7d9348019ad5e63ae09d97f3ca37a68517c77e9a3ac59c |
| TRG | TRGV2*01__homo_sapiens | null | TRGJ1*01__homo_sapiens | TRGC1*01 | GCCAA / none | CATWDGAKNYYKKLF | 79bfde60b69b558ae23781f901de6ab434eb719a3da583462c6d4eeb67ac3c93 |

## Exact boundaries

- **synthetic_TRA**: TRAV1-1*01__homo_sapiens[0:273] -> query[0:273]; n1=GCC -> query[273:276]; TRAJ10*01__homo_sapiens[0:63] -> query[276:339]; TRAC*01[3:123] -> query[339:459].
- **synthetic_TRB**: TRBV10-1*01__homo_sapiens[0:285] -> query[0:285]; n1=GCC -> query[285:288]; TRBD1*01__homo_sapiens[0:12] -> query[288:300]; n2=A -> query[300:301]; TRBJ1-1*01__homo_sapiens[0:47] -> query[301:348]; TRBC1*01[3:123] -> query[348:468].
- **synthetic_TRD**: TRDV1*01__homo_sapiens[0:285] -> query[0:285]; n1=GCC -> query[285:288]; TRDD2*01__homo_sapiens[0:9] -> query[288:297]; n2=A -> query[297:298]; TRDJ1*01__homo_sapiens[0:50] -> query[298:348]; TRDC*01[3:123] -> query[348:468].
- **synthetic_TRG**: TRGV2*01__homo_sapiens[0:300] -> query[0:300]; n1=GCCAA -> query[300:305]; TRGJ1*01__homo_sapiens[0:49] -> query[305:354]; TRGC1*01[3:123] -> query[354:474].

## Allele hashes

| Source ID | Ungapped SHA-256 | IMGT-gapped SHA-256 |
| --- | --- | --- |
| TRAV1-1*01__homo_sapiens | fff121a104864c98536ced43e4faa6503e8366b0b3528c5e81515b0665904c79 | fe41f0edfda9fe447e33af20c9ca8f879fe4480a051d3d77e85732078afe6679 |
| TRAJ10*01__homo_sapiens | 13ce117bd44b2895139c7ba9a0e2988af17aa7cfa3a543298ee3292f3fa8e340 | 13ce117bd44b2895139c7ba9a0e2988af17aa7cfa3a543298ee3292f3fa8e340 |
| TRAC*01 | a5f3276de9175dfcccfb4fe20e4eceeae0756bbc160deededb3dfda6bd2ca537 | 8c64470f70080ab7b8a3e606d622b0588bb5128d8fc5e12698fc493c89a53b58 |
| TRBV10-1*01__homo_sapiens | d4b6a2075b7efab90407ea0c01e36c6c8ed5d0aa5f61a494f94bb096bf2f9fd8 | d6819c1bfa2091a4f53b1fe2b8da812be5945f3fe7f305a7f7cd44a236e794be |
| TRBD1*01__homo_sapiens | 4843114c918f35a391be00552081b29cb0fe5eaef1e1c939b828ec2a347103bb | 4843114c918f35a391be00552081b29cb0fe5eaef1e1c939b828ec2a347103bb |
| TRBJ1-1*01__homo_sapiens | 905c2293af3815a9e41c9435aa0e2fa622eaa4b5d2ce268f27aa16ef0fed4f8f | 905c2293af3815a9e41c9435aa0e2fa622eaa4b5d2ce268f27aa16ef0fed4f8f |
| TRBC1*01 | f8b28bd6a1a532f75dcf9e8f5d1f1c7f52fa898194f31bd9848702debe807953 | b5b76a1e37c17d6214e2ae39134b685086e7be743e00e46ca33f12dbebf66a2d |
| TRDV1*01__homo_sapiens | 85c06eb0d91e5e921631647bc3e2d7882d7a14761fbd0ee5fdf327aa6bfcd23f | 42d11d7a0ac321518066c10cc5054e87d09f10a64ccea27024d8f3d93e307ebb |
| TRDD2*01__homo_sapiens | b58cfd9d646187a48d3c3bc67312499caa0894607c6cb79926b490a005187cfa | b58cfd9d646187a48d3c3bc67312499caa0894607c6cb79926b490a005187cfa |
| TRDJ1*01__homo_sapiens | a63cdeaa2fe0fc6a1800fdb235b3f7e0b01942745c13349de2b21bc002f6e7fa | a63cdeaa2fe0fc6a1800fdb235b3f7e0b01942745c13349de2b21bc002f6e7fa |
| TRDC*01 | eb55f101521b268771d2a3c5e968d40831b9813ed41646584ed2cfbcd965401f | a5a4c0914c0d3decabb53abd5cc2c8296eea1d780f4633411ba8e7b6dbdb54d9 |
| TRGV2*01__homo_sapiens | 5ed12b3eb7324824447278c67bfd9e2387b77397207b2cf995a2bf001096baf7 | 46aa1ead3e4a86e9c84297cb4471a0df997def57feebfa07a256234fa1aa810d |
| TRGJ1*01__homo_sapiens | 3464cc7f4f77d4b24d795179241d9c4f6dab03ed34a00e447af7bdd71f59456b | 3464cc7f4f77d4b24d795179241d9c4f6dab03ed34a00e447af7bdd71f59456b |
| TRGC1*01 | 5f6f55f5e50ba163cc40c2ab418b92088083e88ac1378327a3b843251b6f4add | 43b2c5acf5caab11c49aa7ef4160a40e4616a95c554426f170de035dded82d0e |

## Source file hashes

| File | SHA-256 |
| --- | --- |
| abstar/germline_dbs/tcr/human/manifest.txt | ce4db0ad056fb84dac636537bea6d98ab435de6ac83cfff555d06ff4d3833b7b |
| abstar/germline_dbs/tcr/human/ungapped/v.fasta | 9073979463d502a6b6a5e00d582828ceaab8c6ea07b1008680d3ef1dc9a2b9c5 |
| abstar/germline_dbs/tcr/human/ungapped/d.fasta | 59f4b4793552e339c95af6304823fef561e9e5e5d037ed0ebb1077f67cbe7fd2 |
| abstar/germline_dbs/tcr/human/ungapped/j.fasta | bb003eea50b8e8da24a55e50d8ae061b8959300620a6f829ba664b3d62888a72 |
| abstar/germline_dbs/tcr/human/ungapped/c.fasta | c2c7bbee0879c6014cbbc5374a87000cf3d97d1445835f8d5172dd6d8188dcf2 |
| abstar/germline_dbs/tcr/human/imgt_gapped/v.fasta | 182d9e5279ab5259e374be49f263355980be487bd174593f19872db73da738ec |
| abstar/germline_dbs/tcr/human/imgt_gapped/d.fasta | 59f4b4793552e339c95af6304823fef561e9e5e5d037ed0ebb1077f67cbe7fd2 |
| abstar/germline_dbs/tcr/human/imgt_gapped/j.fasta | bb003eea50b8e8da24a55e50d8ae061b8959300620a6f829ba664b3d62888a72 |
| abstar/germline_dbs/tcr/human/imgt_gapped/c.fasta | 79f445c5d0e3855781005dd561cc280318eb39f6281909db747618209fe82dd9 |

Fixture `sequences.fasta` file SHA-256: `99373274ce2bc48538c83980ce9ab5de19a819a572028a418fb79fd1722e1a36`.

## Reproduction

From the repository root, this standard-library-only program reconstructs all
four FASTA records from the JSON recipe and authenticates the result. The
fixture integrity test additionally authenticates every source file and allele,
checks gapped/ungapped correspondence, proves the anchors and coding frame, and
verifies exact allowed allele ties. It does not invoke abstar or an aligner.

```python
import hashlib
import json
from pathlib import Path

root = Path(".")
data = root / "abstar/tests/data/tcr"
metadata = json.loads((data / "cases.json").read_text())
db = root / "abstar/germline_dbs/tcr/human/ungapped"
references = {}
for segment in ("v", "d", "j", "c"):
    references[segment] = {
        chunk.splitlines()[0].split()[0]: "".join(chunk.splitlines()[1:])
        for chunk in (db / f"{segment}.fasta").read_text().split(">")[1:]
    }
fasta = ""
for case in metadata["cases"]:
    sequence = ""
    for piece in case["construction"]:
        assert len(sequence) == piece["query_start"]
        if "payload" in piece:
            sequence += piece["payload"]
        else:
            source = references[piece["segment"]][piece["allele"]]
            sequence += source[piece["germline_start"]:piece["germline_end"]]
        assert len(sequence) == piece["query_end"]
    assert hashlib.sha256(sequence.encode()).hexdigest() == case["sequence_sha256"]
    fasta += f">{case['sequence_id']}\n{sequence}\n"
assert fasta == (data / "sequences.fasta").read_text()
```

## Packaged BCR database smoke constructions

`test_pipeline.py::test_run_with_each_packaged_bcr_database` defines one input
per database. It reads packaged sources directly, chooses the lexicographically
first V ID and first J ID with compatible locus and species suffix, then uses
`V[:len(V)//3*3] + payload + J[:-1]`. The committed test table records the exact
source IDs, payload and resulting sequence hash. The test authenticates named
V/J membership in both source representations and their gapped/ungapped
identity. GCC supplies an alanine codon; GCCAA adds the two bases required by
the selected mouse J reading frame. Source-based checks establish a stop-free
coding frame and in-frame Cys-to-Trp junction. Expected annotation assertions
are deliberately limited to annotated status, preserved ID, logical database,
expected species/locus, and V/J calls belonging to that database. They do not
claim cross-species biological goldens or cross-species productivity accuracy.

| Database | V source | J source | Payload | Constructed junction AA |
| --- | --- | --- | --- | --- |
| human | IGHV1-18*01__homo_sapiens | IGHJ1*01__homo_sapiens | GCC | CARAAEYFQHW |
| macaque | IGHV1-105*01 | IGHJ1-1*01 | GCC | CARAAEYFEFW |
| c57bl6 | IGHV0-24BS*00__mus_musculus | IGHJ0-32C2*00__mus_musculus | GCCAA | CARANYWYFDVW |
| balbc | IGHV0-22XF*00__mus_musculus | IGHJ0-G76U*00__mus_musculus | GCCAA | CARANYWYFDVW |
| human+c57bl6 | IGHV0-24BS*00__mus_musculus | IGHJ0-32C2*00__mus_musculus | GCCAA | CARANYWYFDVW |

All new API tests redirect HOME to pytest's temporary directory so developer
`~/.abstar` databases cannot shadow these packaged sources.
