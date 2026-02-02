#!/usr/bin/env python3
"""
Fetch protein sequences from UniProt for RBP vocabulary analysis.
Uses direct accession fetch which is more reliable.
"""

import requests
import time
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# config
MAX_WORKERS = 10  # parallel fetches
MAX_RETRIES = 3
DELAY_ON_ERROR = 1.0

PROJECT_DIR = Path("/storage/kiran-stuff/protein-semantic-collapse")
OUTPUT_DIR = PROJECT_DIR / "sequences"
OUTPUT_DIR.mkdir(exist_ok=True)

# force unbuffered output
sys.stdout.reconfigure(line_buffering=True)


def load_ids(filepath):
    """Load UniProt IDs from file."""
    with open(filepath) as f:
        return [line.strip() for line in f if line.strip()]


def fetch_one(accession, retries=MAX_RETRIES):
    """Fetch FASTA for a single UniProt accession."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"

    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return accession, response.text
            elif response.status_code == 404:
                return accession, None  # ID doesn't exist
            elif response.status_code == 429:
                wait = int(response.headers.get("Retry-After", 10))
                time.sleep(wait)
            else:
                time.sleep(DELAY_ON_ERROR)
        except requests.exceptions.RequestException:
            time.sleep(DELAY_ON_ERROR)

    return accession, None


def get_existing_ids(fasta_file):
    """Get set of IDs already in fasta file."""
    existing = set()
    if fasta_file.exists():
        with open(fasta_file) as f:
            for line in f:
                if line.startswith(">"):
                    # >sp|P14866|NAME or >tr|A0A123|NAME
                    parts = line.split("|")
                    if len(parts) >= 2:
                        existing.add(parts[1])
    return existing


def fetch_all(ids, output_file, label):
    """Fetch all sequences with parallel workers."""
    print(f"\n{'='*60}")
    print(f"Fetching {label}: {len(ids)} IDs")
    print(f"{'='*60}")

    # resume support
    existing = get_existing_ids(output_file)
    ids_to_fetch = [id for id in ids if id not in existing]

    if existing:
        print(f"Resuming: {len(existing)} already fetched, {len(ids_to_fetch)} remaining")

    if not ids_to_fetch:
        print("All IDs already fetched!")
        return len(existing), []

    fetched = len(existing)
    failed = []
    start_time = time.time()

    with open(output_file, "a") as out, \
         ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:

        futures = {executor.submit(fetch_one, id): id for id in ids_to_fetch}
        done = 0

        for future in as_completed(futures):
            accession, fasta = future.result()
            done += 1

            if fasta:
                out.write(fasta)
                fetched += 1
            else:
                failed.append(accession)

            # progress every 100
            if done % 100 == 0 or done == len(ids_to_fetch):
                elapsed = time.time() - start_time
                rate = done / elapsed if elapsed > 0 else 0
                eta = (len(ids_to_fetch) - done) / rate if rate > 0 else 0
                print(f"  [{done}/{len(ids_to_fetch)}] {fetched} ok, {len(failed)} failed | "
                      f"{rate:.1f}/s | ETA: {eta/60:.1f}min")

            # flush periodically
            if done % 500 == 0:
                out.flush()

    print(f"\nDone: {fetched} sequences saved")

    if failed:
        failed_file = output_file.with_suffix(".failed.txt")
        with open(failed_file, "w") as f:
            f.write("\n".join(failed))
        print(f"Warning: {len(failed)} IDs failed (obsolete/merged), saved to {failed_file}")

    return fetched, failed


def main():
    start = time.time()

    # load IDs
    rbp_ids = load_ids("/tmp/rbp_uniprot_ids.txt")
    nonrbp_ids = load_ids("/tmp/nonrbp_uniprot_ids.txt")

    print(f"Loaded {len(rbp_ids)} RBP IDs")
    print(f"Loaded {len(nonrbp_ids)} non-RBP IDs")

    # fetch RBPs first (smaller set)
    rbp_file = OUTPUT_DIR / "rbp_sequences.fasta"
    rbp_fetched, rbp_failed = fetch_all(rbp_ids, rbp_file, "RBPs")

    # fetch non-RBPs
    nonrbp_file = OUTPUT_DIR / "nonrbp_sequences.fasta"
    nonrbp_fetched, nonrbp_failed = fetch_all(nonrbp_ids, nonrbp_file, "non-RBPs")

    # summary
    elapsed = time.time() - start
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"RBPs: {rbp_fetched} fetched, {len(rbp_failed)} failed")
    print(f"non-RBPs: {nonrbp_fetched} fetched, {len(nonrbp_failed)} failed")
    print(f"Total time: {elapsed/60:.1f} minutes")
    print(f"Files saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
