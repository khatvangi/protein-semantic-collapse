#!/usr/bin/env python3
"""
Fetch domain sequences from InterPro API.

Strategy:
1. Query InterPro for each Pfam family of interest
2. Get all proteins with that domain + coordinates
3. Match against our RBP/non-RBP protein sets
4. Extract domain sequences

This is faster than querying each protein individually.
"""

import requests
import json
import time
from pathlib import Path
import sys

INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"

def load_protein_ids(fasta_path):
    """extract protein IDs from fasta file, returns set (uppercase)"""
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].split()[0]
                if "|" in parts:
                    prot_id = parts.split("|")[1].upper()
                else:
                    prot_id = parts.upper()
                ids.add(prot_id)
    return ids

def load_sequences(fasta_path):
    """load sequences into dict by ID (uppercase keys)"""
    seqs = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                parts = line[1:].split()[0]
                if "|" in parts:
                    current_id = parts.split("|")[1].upper()
                else:
                    current_id = parts.upper()
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            seqs[current_id] = "".join(current_seq)

    return seqs

def load_pfam_families(path):
    """load Pfam families from file, returns list of (pfam_id, name)"""
    families = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                families.append((parts[0], parts[1]))
    return families

def fetch_pfam_proteins(pfam_id, reviewed_only=True, max_retries=3):
    """
    Fetch all proteins with a given Pfam domain from InterPro API.
    Returns list of dicts: {accession, start, end, protein_length, source_db}.

    Uses pagination to handle large families.
    """
    all_matches = []

    # endpoint: /protein/reviewed/ for SwissProt only, /protein/UniProt/ for all
    if reviewed_only:
        base_url = f"{INTERPRO_API}/protein/reviewed/entry/pfam/{pfam_id}/"
    else:
        base_url = f"{INTERPRO_API}/protein/UniProt/entry/pfam/{pfam_id}/"
    params = "?format=json&page_size=200"

    url = base_url + params
    page = 1

    while url:
        for attempt in range(max_retries):
            try:
                resp = requests.get(url, timeout=120)
                if resp.status_code == 200:
                    data = resp.json()

                    for result in data.get("results", []):
                        meta = result.get("metadata", {})
                        acc = meta.get("accession", "").upper()
                        length = meta.get("length", 0)
                        source = meta.get("source_database", "")

                        # get domain coordinates from entries
                        for entry in result.get("entries", []):
                            if entry.get("accession") == pfam_id:
                                for loc in entry.get("entry_protein_locations", []):
                                    for frag in loc.get("fragments", []):
                                        start = frag.get("start")
                                        end = frag.get("end")
                                        if start and end:
                                            all_matches.append({
                                                "accession": acc,
                                                "start": int(start),
                                                "end": int(end),
                                                "protein_length": length,
                                                "source_db": source
                                            })

                    # pagination
                    url = data.get("next")
                    page += 1
                    break

                elif resp.status_code == 429:
                    time.sleep(10 * (attempt + 1))
                elif resp.status_code == 404:
                    # family not found
                    return all_matches
                else:
                    print(f"  HTTP {resp.status_code} for {pfam_id}", file=sys.stderr)
                    time.sleep(3)

            except requests.exceptions.Timeout:
                print(f"  Timeout for {pfam_id}, retrying...", file=sys.stderr)
                time.sleep(5)
            except Exception as e:
                print(f"  Error fetching {pfam_id}: {e}", file=sys.stderr)
                time.sleep(3)
        else:
            # all retries failed
            print(f"  Failed to fetch {pfam_id} after {max_retries} retries", file=sys.stderr)
            break

        # rate limit delay between pages
        time.sleep(0.5)

    return all_matches

def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")
    data_dir = base_dir / "data"
    seq_dir = base_dir / "sequences"
    out_dir = base_dir / "data"
    out_dir.mkdir(exist_ok=True)

    # load protein sets
    print("Loading protein IDs...")
    rbp_ids = load_protein_ids(seq_dir / "rbp_sequences.fasta")
    nonrbp_ids = load_protein_ids(seq_dir / "nonrbp_sequences.fasta")
    print(f"  RBPs: {len(rbp_ids)}")
    print(f"  non-RBPs: {len(nonrbp_ids)}")

    # load sequences
    print("Loading sequences...")
    rbp_seqs = load_sequences(seq_dir / "rbp_sequences.fasta")
    nonrbp_seqs = load_sequences(seq_dir / "nonrbp_sequences.fasta")
    all_seqs = {**rbp_seqs, **nonrbp_seqs}
    print(f"  Total: {len(all_seqs)} sequences")

    # load RBD Pfam families
    rbd_families = load_pfam_families(data_dir / "rbd_pfam_families.txt")
    print(f"\nLoaded {len(rbd_families)} RBD Pfam families")

    # common non-RBD families for comparison
    # selected for: well-characterized, common, diverse functions
    nonrbd_families = [
        ("PF00069", "Pkinase"),          # protein kinase
        ("PF00017", "SH2"),               # SH2 domain
        ("PF00018", "SH3_1"),             # SH3 domain
        ("PF00046", "Homeodomain"),       # homeodomain
        ("PF00010", "HLH"),               # helix-loop-helix
        ("PF00170", "bZIP_1"),            # bZIP domain
        ("PF00071", "Ras"),               # Ras GTPase
        ("PF00096", "zf-C2H2"),           # C2H2 zinc finger (DNA-binding)
        ("PF00400", "WD40"),              # WD40 repeat
        ("PF00560", "LRR_1"),             # leucine-rich repeat
        ("PF00089", "Trypsin"),           # serine protease
        ("PF07679", "I-set"),             # Ig domain
        ("PF00041", "fn3"),               # fibronectin type III
        ("PF00112", "Peptidase_C1"),      # papain-like protease
        ("PF00102", "Y_phosphatase"),     # tyrosine phosphatase
        ("PF00106", "adh_short"),         # short-chain dehydrogenase
        ("PF00085", "Thioredoxin"),       # thioredoxin
        ("PF00171", "Aldedh"),            # aldehyde dehydrogenase
        ("PF00134", "Cyclin_N"),          # cyclin N-term
        ("PF00023", "Ank"),               # ankyrin repeat
    ]
    print(f"Defined {len(nonrbd_families)} non-RBD Pfam families for comparison")

    # output file
    domain_seqs_file = out_dir / "domain_sequences.jsonl"

    # =========================================================
    # Fetch RBD domains
    # =========================================================
    print("\n" + "=" * 60)
    print("Fetching RBD domain matches from InterPro")
    print("=" * 60)

    rbd_domains = []
    for i, (pfam_id, name) in enumerate(rbd_families):
        print(f"  [{i+1}/{len(rbd_families)}] {pfam_id} ({name})...", end=" ", flush=True)

        # fetch from InterPro (reviewed/SwissProt only)
        matches = fetch_pfam_proteins(pfam_id, reviewed_only=True)

        count_in_set = 0
        for m in matches:
            acc = m["accession"]
            if acc in all_seqs:
                seq = all_seqs[acc]
                # validate coordinates
                if m["start"] >= 1 and m["end"] <= len(seq):
                    dom_seq = seq[m["start"]-1:m["end"]]
                    if len(dom_seq) >= 10:  # min domain length
                        is_rbp = acc in rbp_ids
                        rbd_domains.append({
                            "accession": acc,
                            "pfam": pfam_id,
                            "pfam_name": name,
                            "start": m["start"],
                            "end": m["end"],
                            "domain_length": len(dom_seq),
                            "is_rbp": is_rbp,
                            "is_rbd": True,
                            "sequence": dom_seq
                        })
                        count_in_set += 1

        print(f"{len(matches)} InterPro, {count_in_set} in our set")
        time.sleep(1)  # rate limit between families

    print(f"\nTotal RBD domains: {len(rbd_domains)}")
    print(f"  from RBPs: {sum(1 for d in rbd_domains if d['is_rbp'])}")
    print(f"  from non-RBPs: {sum(1 for d in rbd_domains if not d['is_rbp'])}")

    # =========================================================
    # Fetch non-RBD domains
    # =========================================================
    print("\n" + "=" * 60)
    print("Fetching non-RBD domain matches from InterPro")
    print("=" * 60)

    nonrbd_domains = []
    for i, (pfam_id, name) in enumerate(nonrbd_families):
        print(f"  [{i+1}/{len(nonrbd_families)}] {pfam_id} ({name})...", end=" ", flush=True)

        matches = fetch_pfam_proteins(pfam_id, reviewed_only=True)

        count_in_set = 0
        for m in matches:
            acc = m["accession"]
            if acc in all_seqs:
                seq = all_seqs[acc]
                if m["start"] >= 1 and m["end"] <= len(seq):
                    dom_seq = seq[m["start"]-1:m["end"]]
                    if len(dom_seq) >= 10:
                        is_rbp = acc in rbp_ids
                        nonrbd_domains.append({
                            "accession": acc,
                            "pfam": pfam_id,
                            "pfam_name": name,
                            "start": m["start"],
                            "end": m["end"],
                            "domain_length": len(dom_seq),
                            "is_rbp": is_rbp,
                            "is_rbd": False,
                            "sequence": dom_seq
                        })
                        count_in_set += 1

        print(f"{len(matches)} InterPro, {count_in_set} in our set")
        time.sleep(1)

    print(f"\nTotal non-RBD domains: {len(nonrbd_domains)}")
    print(f"  from RBPs: {sum(1 for d in nonrbd_domains if d['is_rbp'])}")
    print(f"  from non-RBPs: {sum(1 for d in nonrbd_domains if not d['is_rbp'])}")

    # =========================================================
    # Save results
    # =========================================================
    all_domains = rbd_domains + nonrbd_domains
    with open(domain_seqs_file, "w") as f:
        for d in all_domains:
            f.write(json.dumps(d) + "\n")

    print(f"\n" + "=" * 60)
    print("FINAL SUMMARY")
    print("=" * 60)
    print(f"RBD domains: {len(rbd_domains)}")
    print(f"non-RBD domains: {len(nonrbd_domains)}")
    print(f"Total: {len(all_domains)}")
    print(f"\nSaved to: {domain_seqs_file}")

if __name__ == "__main__":
    main()
