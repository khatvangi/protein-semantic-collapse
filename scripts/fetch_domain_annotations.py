#!/usr/bin/env python3
"""
Fetch domain annotations from UniProt for all proteins.
Extracts Pfam domain coordinates to enable domain-level analysis.

Uses UniProt REST API with batching for efficiency.
Saves progress for resumability.
"""

import requests
import json
import time
from pathlib import Path
import sys

BASE_URL = "https://rest.uniprot.org/uniprotkb"

def load_ids(fasta_path):
    """extract protein IDs from fasta file"""
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # format: >sp|P12345|NAME or >P12345
                parts = line[1:].split()[0]
                if "|" in parts:
                    prot_id = parts.split("|")[1]
                else:
                    prot_id = parts
                ids.append(prot_id)
    return ids

def load_sequences(fasta_path):
    """load sequences into dict"""
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
                    current_id = parts.split("|")[1]
                else:
                    current_id = parts
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = "".join(current_seq)

    return seqs

def fetch_domains_batch(ids, max_retries=3):
    """
    Fetch domain annotations for a batch of IDs.
    Returns dict: {id: [(pfam_id, name, start, end), ...]}
    """
    results = {}

    # UniProt query for multiple IDs
    query = " OR ".join(f"accession:{id}" for id in ids)
    params = {
        "query": query,
        "fields": "accession,ft_domain",
        "format": "json",
        "size": len(ids)
    }

    for attempt in range(max_retries):
        try:
            resp = requests.get(BASE_URL + "/search", params=params, timeout=60)
            if resp.status_code == 200:
                data = resp.json()

                for entry in data.get("results", []):
                    acc = entry.get("primaryAccession", "")
                    domains = []

                    # extract domain features
                    for feature in entry.get("features", []):
                        if feature.get("type") == "Domain":
                            loc = feature.get("location", {})
                            start = loc.get("start", {}).get("value")
                            end = loc.get("end", {}).get("value")
                            desc = feature.get("description", "")

                            # try to extract Pfam ID from cross-references
                            pfam_id = ""
                            for xref in feature.get("evidences", []):
                                if xref.get("source") == "Pfam":
                                    pfam_id = xref.get("id", "")
                                    break

                            if start and end:
                                domains.append({
                                    "pfam": pfam_id,
                                    "name": desc,
                                    "start": int(start),
                                    "end": int(end)
                                })

                    results[acc] = domains

                return results

            elif resp.status_code == 429:
                # rate limited
                time.sleep(5 * (attempt + 1))
            else:
                print(f"  HTTP {resp.status_code} for batch", file=sys.stderr)
                time.sleep(2)

        except Exception as e:
            print(f"  Error: {e}", file=sys.stderr)
            time.sleep(2)

    return results

def main():
    base_dir = Path("/storage/kiran-stuff/protein-semantic-collapse")
    seq_dir = base_dir / "sequences"
    out_dir = base_dir / "data"
    out_dir.mkdir(exist_ok=True)

    # output files
    domains_file = out_dir / "domain_annotations.jsonl"
    progress_file = out_dir / "domain_fetch_progress.txt"

    # load all protein IDs
    print("Loading protein IDs...")
    rbp_ids = load_ids(seq_dir / "rbp_sequences.fasta")
    nonrbp_ids = load_ids(seq_dir / "nonrbp_sequences.fasta")

    all_ids = [(id, "RBP") for id in rbp_ids] + [(id, "non-RBP") for id in nonrbp_ids]
    print(f"  Total: {len(all_ids)} proteins")

    # check progress
    done_ids = set()
    if progress_file.exists():
        with open(progress_file) as f:
            done_ids = set(line.strip() for line in f)
        print(f"  Resuming: {len(done_ids)} already done")

    # filter to remaining
    remaining = [(id, label) for id, label in all_ids if id not in done_ids]
    print(f"  Remaining: {len(remaining)} to fetch")

    if not remaining:
        print("All done!")
        return

    # batch fetch
    batch_size = 25  # UniProt recommends small batches
    total_batches = (len(remaining) + batch_size - 1) // batch_size

    start_time = time.time()
    fetched = 0

    with open(domains_file, "a") as out_f, open(progress_file, "a") as prog_f:
        for i in range(0, len(remaining), batch_size):
            batch = remaining[i:i+batch_size]
            batch_ids = [id for id, _ in batch]
            batch_labels = {id: label for id, label in batch}

            results = fetch_domains_batch(batch_ids)

            # write results
            for prot_id in batch_ids:
                domains = results.get(prot_id, [])
                record = {
                    "id": prot_id,
                    "label": batch_labels[prot_id],
                    "domains": domains
                }
                out_f.write(json.dumps(record) + "\n")
                prog_f.write(prot_id + "\n")

            fetched += len(batch)

            # progress
            if (i // batch_size + 1) % 100 == 0 or i + batch_size >= len(remaining):
                elapsed = time.time() - start_time
                rate = fetched / elapsed if elapsed > 0 else 0
                eta = (len(remaining) - fetched) / rate / 60 if rate > 0 else 0
                batch_num = i // batch_size + 1

                print(f"  [{fetched}/{len(remaining)}] {rate:.1f}/s | ETA: {eta:.1f}min")
                out_f.flush()
                prog_f.flush()

            # small delay to respect rate limits
            time.sleep(0.3)

    print(f"\nDone! Saved to {domains_file}")

if __name__ == "__main__":
    main()
