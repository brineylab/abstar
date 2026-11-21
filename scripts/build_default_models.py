"""Build the built-in human recombination models from digitized published distributions.

Usage: python scripts/build_default_models.py

Inputs are the digitized tables in germline-sources/models/digitized/ (see
germline-sources/models/SOURCES.md for every citation and extraction method) and the
built-in human germline database, which decides which genes can be generated. Output:
python/abstar/models/<name>/model.toml. Standard library only.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIG = ROOT / "germline-sources" / "models" / "digitized"
OUT = ROOT / "python" / "abstar" / "models"
DB = ROOT / "python" / "abstar" / "germline_dbs" / "human" / "alleles.tsv"

IGOR = (
    "Marcou Q, Mora T, Walczak AM (2018) High-throughput immune repertoire analysis with IGoR. "
    "Nat Commun 9:561. doi:10.1038/s41467-018-02832-w. CC BY 4.0 "
    "(http://creativecommons.org/licenses/by/4.0/); values digitized from the published vector "
    "figure; changes: renormalized and restricted to genes in the abstar human database"
)
ELHANATI = (
    "Elhanati Y, Sethna Z, Marcou Q, Callan CG Jr, Mora T, Walczak AM (2015) Inferring processes "
    "underlying B-cell repertoire diversity. Phil Trans R Soc B 370:20140243. "
    "doi:10.1098/rstb.2014.0243. CC BY 4.0; values digitized from the electronic supplementary "
    "material (vector figures); changes: averaged over individuals A and B, renormalized"
)
MIKOCZIOVA = (
    "Mikocziova I, Peres A, Gidoni M, Greiff V, Yaari G, Sollid LM (2021) Germline polymorphisms "
    "and alternative splicing of human immunoglobulin light chain genes. iScience 24:103192. "
    "doi:10.1016/j.isci.2021.103192. CC BY 4.0; box-plot medians digitized from the high-resolution "
    "raster figure (about +/-0.003); changes: renormalized, restricted to database genes"
)
TOWNSEND = (
    "Townsend CL, Laffy JMJ, Wu YB, et al. (2016) Significant differences in physicochemical "
    "properties of human immunoglobulin kappa and lambda CDR3 regions. Front Immunol 7:388. "
    "doi:10.3389/fimmu.2016.00388. CC BY 4.0; mean N additions read from the Results text"
)
TAIL = (
    "Engineered table, not from a published distribution: uniform 20-60 nt per junction, mixed in "
    "by the bnab preset's long_cdr3_weight. Rationale: Briney BS, Willis JR, Crowe JE Jr (2012) "
    "PLoS ONE 7:e36750 (CC BY), Fig. 3A, shows that long HCDR3s come mainly from long N additions"
)

# Mikocziova et al. label identical duplicate genes "E" (either copy). Map to the database names.
LIGHT_E = {
    "IGKV1E-12": ["IGKV1-12", "IGKV1D-12"], "IGKV1E-13": ["IGKV1-13", "IGKV1D-13"],
    "IGKV6E-21": ["IGKV6-21", "IGKV6D-21"], "IGKV2E-28": ["IGKV2-28", "IGKV2D-28"],
    "IGKV1E-33": ["IGKV1-33", "IGKV1D-33"], "IGKV1E-37": ["IGKV1-37", "IGKV1D-37"],
    "IGKV1E-39": ["IGKV1-39", "IGKV1D-39"], "IGKV2E-40": ["IGKV2-40", "IGKV2D-40"],
}
# Box plots at zero (no box drawn: "Unknown" or deleted in most donors) get this floor so the
# gene is still generated rarely. Genes the paper calls likely pseudogenes get no floor.
FLOOR = 0.001
PSEUDO = {"IGLJ6", "IGLV3-12"}


def load(name):
    return json.loads((DIG / f"{name}.json").read_text())


def generatable():
    """locus -> segment -> genes with a generatable allele (V and J need an anchor)."""
    out = {}
    with open(DB) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["species"] == "human" and (r["segment"] == "D" or r["anchor"]):
                out.setdefault(r["locus"], {}).setdefault(r["segment"], set()).add(r["gene"])
    return out


def norm(d):
    s = sum(d.values())
    return {k: v / s for k, v in d.items()}


def restrict(usage, genes, notes, what):
    """Keep genes in the database; record dropped genes and the source mass they carried."""
    usage = norm(usage)
    drop = sorted(g for g in usage if g not in genes)
    if drop:
        notes.append(f"{what}: dropped (not in database) {', '.join(drop)}; "
                     f"{sum(usage[g] for g in drop):.4f} of source mass")
    return norm({g: p for g, p in usage.items() if g in genes})


def dist(d):
    """{int-like: p} -> (values, probs), dropping zeros, renormalized."""
    d = norm({int(k): v for k, v in d.items() if v > 0})
    ks = sorted(d)
    return ks, [d[k] for k in ks]


def avg(*ds):
    keys = set().union(*ds)
    return {k: sum(d.get(k, 0.0) for d in ds) / len(ds) for k in keys}


def geometric(mean, eps=1e-4):
    p = 1 / (1 + mean)
    out, n = {}, 0
    while True:
        v = p * (1 - p) ** n
        if v < eps:
            return out
        out[n] = v
        n += 1


def n_markov():
    e = load("igh_elhanati_s9_s10")
    rows = {"start": [(a + b) / 2 for a, b in zip(e["start_A"], e["start_B"])]}
    for x in "ACGT":
        rows[x] = [(a + b) / 2 for a, b in zip(e["trans_A"][x], e["trans_B"][x])]
    return {k: [v / sum(r) for v in r] for k, r in rows.items()}


def fmt(x):
    return f"{x:.8f}".rstrip("0").rstrip(".") if x else "0"


def q(s):
    return json.dumps(s)


def write(m):
    L = [f"# Generated by scripts/build_default_models.py; edit the sources, not this file.",
         f"# Citations and extraction methods: germline-sources/models/SOURCES.md"]
    L += [f"# {n}" for n in m["notes"]]
    L += ["format_version = 1", f"name = {q(m['name'])}", 'species = "human"', f"locus = {q(m['locus'])}",
          f"description = {q(m['description'])}", "d_inversion = 0.0", f"shm_model = {q(m['shm'])}", ""]
    for cit, used in m["sources"]:
        L += ["[[source]]", f"citation = {q(cit)}", f"used_for = {q(used)}", ""]

    def table(header, d):
        L.append(f"[{header}]")
        L.extend(f"{q(k)} = {fmt(v)}" for k, v in sorted(d.items()))
        L.append("")

    table("v_usage", m["v"])
    table("j_usage", m["j"])
    for j, row in sorted(m.get("d_given_j", {}).items()):
        table(f"d_given_j.{q(j)}", row)

    def dtab(header, d):
        ks, ps = dist(d)
        L.extend([f"[{header}]", f"values = [{', '.join(map(str, ks))}]",
                  f"probs = [{', '.join(fmt(p) for p in ps)}]", ""])

    for key, d in m["trim"].items():
        dtab(f"trim.{key}.default", d)
    for key, d in m["ins"].items():
        dtab(f"insertions.{key}", d)
    for key, d in m.get("tail", {}).items():
        dtab(f"insertions_long_tail.{key}", d)
    L.append("[n_markov]")
    L.extend(f"{k} = [{', '.join(fmt(v) for v in r)}]" for k, r in n_markov().items())
    path = OUT / m["name"] / "model.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L) + "\n")
    print(f"{path.relative_to(ROOT)}")
    for n in m["notes"]:
        print(f"  {n}")


N_COMP = (ELHANATI, "n_markov: start = Fig. S9 'VD sense' mononucleotide frequencies; rows = Fig. S10 "
          "dinucleotide (VD sense) conditional frequencies")


def igh(db):
    s13, s4, notes = load("igh_igor_s13"), load("igh_igor_s4"), []
    d_given_j = {j: restrict(row, db["IGH"]["D"], notes if j == "IGHJ1" else [], "D")
                 for j, row in s4.items()}
    tail = {n: 1 / 41 for n in range(20, 61)}
    return dict(
        name="human-igh", locus="IGH", shm="motif", notes=notes,
        description="Human IGH; gene usage, P(D|J), trims and N lengths from IGoR (Marcou et al. 2018) "
                    "digitized figures; N composition from Elhanati et al. 2015.",
        v=restrict(s13["V"], db["IGH"]["V"], notes, "V"), j=restrict(s13["J"], db["IGH"]["J"], notes, "J"),
        d_given_j=d_given_j,
        trim={"v_3": s13["del_v3"], "d_5": s13["del_d5"], "d_3": s13["del_d3"], "j_5": s13["del_j5"]},
        ins={"vd": s13["ins_vd"], "dj": s13["ins_dj"]}, tail={"vd": tail, "dj": tail},
        sources=[
            (IGOR, "v_usage, j_usage: Supplementary Fig. 13 IGHV/IGHJ panels, IGoR series"),
            (IGOR, "d_given_j: Supplementary Fig. 4 (P(D|J) per IGHD gene)"),
            (IGOR, "trim v_3, d_5, d_3, j_5 and insertions vd, dj: Supplementary Fig. 13, IGoR series "
                   "(deletions -5..15, insertions 0..40 as plotted)"),
            N_COMP,
            (TAIL, "insertions_long_tail vd, dj"),
        ],
    )


def trb(db):
    f2, s3, s7, notes = load("tcr_igor_fig2"), load("trb_igor_s3"), load("trb_igor_s7_s14"), []

    def gene_usage(prefix):
        u = {}
        for allele, v in s3.items():
            if allele.startswith(prefix):
                g = allele.split("*")[0]
                u[g] = u.get(g, 0.0) + (sum(v["dna"]) / 2 + sum(v["rna"]) / 5) / 2
        return u

    p2 = s7["trb_p_d2_given_j"]
    return dict(
        name="human-trb", locus="TRB", shm="none", notes=notes,
        description="Human TRB; gene usage, P(D|J), trims and N lengths from IGoR (Marcou et al. 2018) "
                    "digitized figures; N composition borrowed from IGH (Elhanati et al. 2015).",
        v=restrict(gene_usage("TRBV"), db["TRB"]["V"], notes, "V"),
        j=restrict(gene_usage("TRBJ"), db["TRB"]["J"], notes, "J"),
        d_given_j={j: {"TRBD1": 1 - p2[j], "TRBD2": p2[j]} for j in p2},
        trim={"v_3": avg(f2["trb_v3_dna"], f2["trb_v3_rna"]), "d_5": s7["trb_d5"], "d_3": s7["trb_d3"],
              "j_5": avg(f2["trb_j5_dna"], f2["trb_j5_rna"])},
        ins={"vd": avg(f2["trb_vd_dna"], f2["trb_vd_rna"]), "dj": avg(f2["trb_dj_dna"], f2["trb_dj_rna"])},
        sources=[
            (IGOR, "v_usage, j_usage: Supplementary Fig. 3 (allele usage summed to genes; mean of the 2 DNA "
                   "individuals and mean of the 5 mRNA individuals, averaged)"),
            (IGOR, "d_given_j: Supplementary Fig. 14, IGoR P(D2|J); P(D1|J) = 1 - P(D2|J)"),
            (IGOR, "trim v_3, j_5 and insertions vd, dj: Fig. 2a-c, mean of the TRB DNA and TRB RNA curves"),
            (IGOR, "trim d_5, d_3: Supplementary Fig. 7, 'true model' curves (the data-inferred TRB model)"),
            (N_COMP[0], "n_markov (surrogate: IGH TdT composition; no TRB N composition in the sources used)"),
        ],
    )


def tra(db):
    f2, notes = load("tcr_igor_fig2"), []
    top_v = {"TRAV19": 0.08, "TRAV13-1": 0.13, "TRAV27": 0.06, "TRAV13-2": 0.05, "TRAV4": 0.04}
    top_j = {g: 0.04 for g in ["TRAJ45", "TRAJ49", "TRAJ37", "TRAJ53", "TRAJ42", "TRAJ48"]}

    def fill(top, genes, what):
        rest = sorted(genes - set(top))
        notes.append(f"{what}: {len(top)} genes from Supplementary Fig. 5 (mass {sum(top.values()):.2f}); "
                     f"fallback: remaining mass uniform over the other {len(rest)} database genes")
        r = (1 - sum(top.values())) / len(rest)
        return {**top, **{g: r for g in rest}}

    return dict(
        name="human-tra", locus="TRA", shm="none", notes=notes,
        description="Human TRA; trims and N lengths from IGoR (Marcou et al. 2018) digitized figures; "
                    "gene usage partly from IGoR (top genes) with a uniform fallback for the rest.",
        v=fill(top_v, db["TRA"]["V"], "V"), j=fill(top_j, db["TRA"]["J"], "J"),
        trim={"v_3": f2["tra_v3_rna"], "j_5": f2["tra_j5_rna"]}, ins={"vj": f2["tra_vj_rna"]},
        sources=[
            (IGOR, "v_usage, j_usage (partial): Supplementary Fig. 5 panel labels P(V), P(J) for the 6 most "
                   "used TRAV alleles (TRAV13-1*01 and *02 summed) and TRAJ alleles, read from the figure text"),
            ("FALLBACK: uniform over the remaining generatable genes in the abstar human database",
             "v_usage, j_usage (all genes not in Supplementary Fig. 5)"),
            (IGOR, "trim v_3, j_5 and insertions vj: Fig. 2a-c, TRA RNA curves"),
            (N_COMP[0], "n_markov (surrogate: IGH TdT composition; no TRA N composition in the sources used)"),
        ],
    )


def light(db, locus):
    mk, s13, notes = load("light_mikocziova_fig3_fig4"), load("igh_igor_s13"), []
    fig = "Fig. 4" if locus == "IGK" else "Fig. 3"
    n_mean = 3.046 if locus == "IGK" else 3.321

    def usage(seg):
        u = {}
        for g, v in mk[f"{locus}{seg}"].items():
            if g in PSEUDO:
                continue
            v = FLOOR if v is None or v < FLOOR else v
            names = [x for x in LIGHT_E.get(g, [g]) if x in db[locus][seg]] or LIGHT_E.get(g, [g])[:1]
            for x in names:
                u[x] = u.get(x, 0.0) + v / len(names)
        return restrict(u, db[locus][seg], notes, seg)

    return dict(
        name=f"human-{locus.lower()}", locus=locus, shm="motif", notes=notes,
        description=f"Human {locus}; gene usage from Mikocziova et al. 2021 {fig} (naive repertoires, "
                    "100 donors); N length geometric with the Townsend et al. 2016 mean; trims and N "
                    "composition are IGH surrogates.",
        v=usage("V"), j=usage("J"),
        trim={"v_3": s13["del_v3"], "j_5": s13["del_j5"]}, ins={"vj": geometric(n_mean)},
        sources=[
            (MIKOCZIOVA, f"v_usage, j_usage: {fig} box-plot medians; duplicate 'E' genes split equally over "
                         f"the database copies; boxes at zero floored at {FLOOR}; likely pseudogenes "
                         f"({', '.join(sorted(PSEUDO))}) omitted"),
            (TOWNSEND, f"insertions vj (FALLBACK shape): geometric with mean {n_mean} nt "
                       f"({'kappa' if locus == 'IGK' else 'lambda'} mean N addition), truncated at p < 1e-4"),
            (IGOR, "trim v_3, j_5 (FALLBACK surrogate: IGH V 3' and J 5' deletion profiles from "
                   "Supplementary Fig. 13; no open published light-chain trim distribution found)"),
            (N_COMP[0], "n_markov (surrogate: IGH TdT composition)"),
        ],
    )


def main():
    db = generatable()
    for m in [igh(db), light(db, "IGK"), light(db, "IGL"), tra(db), trb(db)]:
        write(m)


if __name__ == "__main__":
    main()
