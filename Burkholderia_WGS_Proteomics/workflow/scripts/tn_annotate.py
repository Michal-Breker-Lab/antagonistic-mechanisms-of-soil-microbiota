"""CDS lookup shared by the WGS and Sanger insertion callers."""
from collections import defaultdict

def parse_gff(path):
    feats = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            tag = a.get("locus_tag")
            if not tag:
                continue
            feats[f[0]].append({
                "locus_tag": tag,
                "product": a.get("product", "").replace("%2C", ","),
                "gene": a.get("gene", ""),
                "start": int(f[3]), "end": int(f[4]), "strand": f[6],
            })
    for c in feats:
        feats[c].sort(key=lambda d: d["start"])
    return feats


def annotate(feats, contig, pos):
    genes = feats.get(contig, [])
    hit = next((g for g in genes if g["start"] <= pos <= g["end"]), None)
    if hit:
        length = hit["end"] - hit["start"] + 1
        into = (pos - hit["start"] + 1) if hit["strand"] == "+" else (hit["end"] - pos + 1)
        pct = 100.0 * into / length
        codon = (into + 2) // 3
        eff = ("5' 5% of CDS - null expected" if pct <= 5 else
               "3' 10% of CDS - may retain function" if pct >= 90 else
               "mid-CDS - null expected")
        return {"locus_tag": hit["locus_tag"], "gene": hit["gene"],
                "product": hit["product"], "cds_start": hit["start"],
                "cds_end": hit["end"], "cds_strand": hit["strand"],
                "cds_len_bp": length, "bp_into_cds": into,
                "pct_into_cds": round(pct, 1), "codon": codon,
                "aa_pos_of": (length // 3) - 1, "effect": eff,
                "context": "in CDS"}
    prev = [g for g in genes if g["end"] < pos]
    nxt = [g for g in genes if g["start"] > pos]
    L = prev[-1] if prev else None
    R = nxt[0] if nxt else None
    bits = []
    if L:
        bits.append(f"{L['locus_tag']} ends {pos - L['end']} bp upstream ({L['strand']})")
    if R:
        bits.append(f"{R['locus_tag']} starts {R['start'] - pos} bp downstream ({R['strand']})")
    prom = [g for g in (R,) if g and g["strand"] == "+" and g["start"] - pos <= 200]
    prom += [g for g in (L,) if g and g["strand"] == "-" and pos - g["end"] <= 200]
    return {"locus_tag": "", "gene": "", "product": "", "cds_start": "",
            "cds_end": "", "cds_strand": "", "cds_len_bp": "", "bp_into_cds": "",
            "pct_into_cds": "", "codon": "", "aa_pos_of": "",
            "effect": ("intergenic, within 200 bp 5' of "
                       + ",".join(g["locus_tag"] for g in prom)
                       + " - possible promoter disruption") if prom else "intergenic",
            "context": "; ".join(bits)}
