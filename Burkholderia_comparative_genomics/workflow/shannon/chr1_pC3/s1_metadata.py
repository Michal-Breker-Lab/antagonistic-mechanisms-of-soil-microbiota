#!/usr/bin/env python3
"""Stage 1 - acquire assembly + replicon metadata for all complete Burkholderia
sensu lato genomes from NCBI Datasets v2alpha.

Two traps this script exists to avoid:

1. GCA/GCF double counting. A naive taxon query returns both the GenBank (GCA)
   and RefSeq (GCF) record for the same physical genome. For Burkholderia the
   raw count is 1283 = 612 refseq + 650 genbank; the distinct-genome count is
   roughly half that. We keep GCF when a pair exists, GCA when it does not.

2. Replicon data is not in dataset_report. Per-sequence molecule types and
   lengths come from the separate /sequence_reports endpoint, one call per
   assembly.
"""
import json
import os
import sys
import time
import urllib.error
import urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
GENERA = ["Burkholderia", "Paraburkholderia", "Caballeronia",
          "Trinickia", "Mycetohabitans", "Robbsia"]
OUT = os.environ.get("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
MD = f"{OUT}/metadata"
NCBI_KEY = os.environ.get("NCBI_API_KEY")          # optional, raises rate limit
DELAY = 0.12 if NCBI_KEY else 0.40                  # 10/s with key, ~3/s without


def get(url, tries=5):
    if NCBI_KEY:
        url += ("&" if "?" in url else "?") + "api_key=" + NCBI_KEY
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return json.load(fh)
        except (urllib.error.URLError, urllib.error.HTTPError,
                TimeoutError, json.JSONDecodeError) as e:
            if a == tries - 1:
                raise
            time.sleep(2 ** a)
    raise RuntimeError("unreachable")


def fetch_reports():
    """All complete-genome dataset reports across the six genera."""
    seen, out = set(), []
    for g in GENERA:
        tok, n = None, 0
        while True:
            u = (f"{API}/genome/taxon/{g}/dataset_report"
                 f"?filters.assembly_level=complete_genome"
                 f"&page_size=1000&returned_content=COMPLETE")
            if tok:
                u += f"&page_token={tok}"
            d = get(u)
            for r in d.get("reports", []):
                acc = r["accession"]
                if acc in seen:      # genera overlap via reclassified taxa
                    continue
                seen.add(acc)
                r["_query_genus"] = g
                out.append(r)
                n += 1
            tok = d.get("next_page_token")
            if not tok:
                break
            time.sleep(DELAY)
        print(f"  {g:18s} {n:5d} records", flush=True)
    return out


def dedupe(reports):
    """Keep GCF over its paired GCA. Returns (kept, dropped_count)."""
    by_acc = {r["accession"]: r for r in reports}
    drop = set()
    for r in reports:
        acc = r["accession"]
        pair = (r.get("assembly_info", {}).get("paired_assembly", {}) or {}).get("accession")
        if not pair or pair not in by_acc:
            continue
        # drop whichever of the pair is the GenBank (GCA) record
        loser = acc if acc.startswith("GCA_") else pair
        drop.add(loser)
    kept = [r for r in reports if r["accession"] not in drop]
    return kept, len(drop)


def biosample_fields(bs):
    """Flatten BioSample: named fields plus every attribute key/value."""
    out = {}
    if not bs:
        return out
    for k in ("accession", "host", "isolation_source", "geo_loc_name",
              "collection_date", "collected_by", "lat_lon", "strain",
              "package", "publication_date"):
        v = bs.get(k)
        if v:
            out["bs_" + k] = str(v).replace("\t", " ").replace("\n", " ")
    for a in bs.get("attributes", []) or []:
        name, val = a.get("name"), a.get("value")
        if name and val and ("bs_" + name) not in out:
            out["attr_" + name] = str(val).replace("\t", " ").replace("\n", " ")
    return out


def flatten(r):
    ai = r.get("assembly_info", {}) or {}
    st = r.get("assembly_stats", {}) or {}
    cm = r.get("checkm_info", {}) or {}
    an = r.get("average_nucleotide_identity", {}) or {}
    best = an.get("best_ani_match", {}) or {}
    org = r.get("organism", {}) or {}
    row = {
        "accession": r["accession"],
        "paired_accession": ai.get("paired_assembly", {}).get("accession", ""),
        "source_database": r.get("source_database", ""),
        "query_genus": r.get("_query_genus", ""),
        "organism_name": org.get("organism_name", ""),
        "tax_id": org.get("tax_id", ""),
        "strain": (org.get("infraspecific_names", {}) or {}).get("strain", ""),
        "assembly_level": ai.get("assembly_level", ""),
        "assembly_name": ai.get("assembly_name", ""),
        "release_date": ai.get("release_date", ""),
        "submitter": ai.get("submitter", ""),
        "refseq_category": ai.get("refseq_category", ""),
        "sequencing_tech": ai.get("sequencing_tech", ""),
        "assembly_method": ai.get("assembly_method", ""),
        "bioproject": ai.get("bioproject_accession", ""),
        "ncbi_n_chromosomes": st.get("total_number_of_chromosomes", ""),
        "total_length": st.get("total_sequence_length", ""),
        "n_contigs": st.get("number_of_contigs", ""),
        "contig_n50": st.get("contig_n50", ""),
        "gc_percent": st.get("gc_percent", ""),
        "genome_coverage": st.get("genome_coverage", ""),
        "checkm_completeness": cm.get("completeness", ""),
        "checkm_contamination": cm.get("contamination", ""),
        "ani_taxonomy_check": an.get("taxonomy_check_status", ""),
        "ani_match_status": an.get("match_status", ""),
        "ani_best_organism": best.get("organism_name", ""),
        "ani_best_value": best.get("ani", ""),
        "ani_best_category": best.get("category", ""),
    }
    row.update(biosample_fields(ai.get("biosample")))
    return row


def write_tsv(rows, path):
    cols, seen = [], set()
    for r in rows:                       # stable union of keys, first-seen order
        for k in r:
            if k not in seen:
                seen.add(k)
                cols.append(k)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    return cols


def fetch_replicons(accs):
    rows, miss = [], []
    for i, acc in enumerate(accs, 1):
        try:
            d = get(f"{API}/genome/accession/{acc}/sequence_reports?page_size=1000")
        except Exception as e:                       # noqa: BLE001
            miss.append((acc, repr(e)))
            continue
        reps = d.get("reports", []) or []
        if not reps:
            miss.append((acc, "empty"))
        for s in reps:
            rows.append({
                "accession": acc,
                "sequence_name": s.get("sequence_name", ""),
                "chr_name": s.get("chr_name", ""),
                "genbank_accession": s.get("genbank_accession", ""),
                "refseq_accession": s.get("refseq_accession", ""),
                "length": s.get("length", ""),
                "gc_percent": s.get("gc_percent", ""),
                "molecule_type": s.get("assigned_molecule_location_type", ""),
                "role": s.get("role", ""),
                "assembly_unit": s.get("assembly_unit", ""),
            })
        if i % 100 == 0:
            print(f"  replicons {i}/{len(accs)}", flush=True)
        time.sleep(DELAY)
    return rows, miss


def main():
    os.makedirs(MD, exist_ok=True)
    print("=== fetching dataset reports ===", flush=True)
    reports = fetch_reports()
    print(f"raw records (accession-unique): {len(reports)}")

    kept, dropped = dedupe(reports)
    print(f"dropped as GCA duplicates of a GCF: {dropped}")
    print(f"DISTINCT GENOMES: {len(kept)}")
    n_gcf = sum(1 for r in kept if r["accession"].startswith("GCF_"))
    print(f"  of which RefSeq (GCF): {n_gcf}   GenBank-only (GCA): {len(kept)-n_gcf}")

    with open(f"{MD}/dataset_reports.json", "w") as fh:
        json.dump(kept, fh)

    rows = [flatten(r) for r in kept]
    cols = write_tsv(rows, f"{MD}/assemblies.tsv")
    print(f"wrote assemblies.tsv  {len(rows)} rows x {len(cols)} cols")

    accs = [r["accession"] for r in kept]
    with open(f"{MD}/accessions.txt", "w") as fh:
        fh.write("\n".join(accs) + "\n")

    print("=== fetching sequence reports (replicons) ===", flush=True)
    rrows, miss = fetch_replicons(accs)
    write_tsv(rrows, f"{MD}/replicons.tsv")
    print(f"wrote replicons.tsv   {len(rrows)} rows for "
          f"{len(set(r['accession'] for r in rrows))} assemblies")
    if miss:
        with open(f"{MD}/replicons_missing.txt", "w") as fh:
            for a, e in miss:
                fh.write(f"{a}\t{e}\n")
        print(f"WARNING: {len(miss)} assemblies returned no sequence report")

    # quick sanity summary
    host = sum(1 for r in rows if r.get("bs_host"))
    iso = sum(1 for r in rows if r.get("bs_isolation_source"))
    print(f"\nBioSample coverage: host={host}/{len(rows)}  "
          f"isolation_source={iso}/{len(rows)}")
    print("STAGE1_DONE")


if __name__ == "__main__":
    sys.exit(main())
