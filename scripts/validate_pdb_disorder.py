#!/usr/bin/env python3
"""
Structural validation: Check if D/E-rich regions correspond to
disordered regions in PDB structures.

Uses RCSB PDB API to fetch structure quality metrics.
"""

import json
import requests
from pathlib import Path

# families to validate
HIGH_DE_FAMILIES = ["PF00658", "PF00012", "PF00076"]  # PABP, HSP70, RRM_1
LOW_DE_FAMILIES = ["PF00189", "PF01479", "PF01280"]   # Ribosomal

def search_pdb_by_pfam(pfam_id, limit=5):
    """Search PDB for structures containing a Pfam domain."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "exact_match",
                "value": pfam_id
            }
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "return_all_hits": False,
            "pager": {"start": 0, "rows": limit}
        }
    }

    try:
        resp = requests.post(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            json=query,
            timeout=30
        )
        if resp.status_code == 200:
            data = resp.json()
            return [hit["identifier"] for hit in data.get("result_set", [])]
        else:
            return []
    except Exception as e:
        print(f"  Error searching {pfam_id}: {e}")
        return []


def get_structure_quality(pdb_id):
    """Get resolution and missing residue info for a structure."""
    try:
        # get entry info
        resp = requests.get(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}",
            timeout=30
        )
        if resp.status_code != 200:
            return None

        data = resp.json()

        # resolution
        resolution = None
        if "rcsb_entry_info" in data:
            resolution = data["rcsb_entry_info"].get("resolution_combined", [None])[0]

        # get polymer entity info for missing residues
        resp2 = requests.get(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1",
            timeout=30
        )

        n_missing = 0
        total_residues = 0
        if resp2.status_code == 200:
            poly_data = resp2.json()
            if "entity_poly" in poly_data:
                seq = poly_data["entity_poly"].get("pdbx_seq_one_letter_code_can", "")
                total_residues = len(seq) if seq else 0

            # check for unobserved residues
            if "rcsb_polymer_entity_feature" in poly_data:
                for feat in poly_data["rcsb_polymer_entity_feature"]:
                    if feat.get("type") == "unobserved":
                        n_missing += 1

        return {
            "pdb_id": pdb_id,
            "resolution": resolution,
            "total_residues": total_residues,
            "has_missing": n_missing > 0
        }

    except Exception as e:
        print(f"  Error getting quality for {pdb_id}: {e}")
        return None


def get_unobserved_regions(pdb_id):
    """Get detailed unobserved (disordered) regions from PDB."""
    try:
        resp = requests.get(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id}/A",
            timeout=30
        )
        if resp.status_code != 200:
            return []

        data = resp.json()
        unobserved = []

        if "rcsb_polymer_instance_feature" in data:
            for feat in data["rcsb_polymer_instance_feature"]:
                if "UNOBSERVED" in feat.get("type", "").upper():
                    for region in feat.get("feature_positions", []):
                        start = region.get("beg_seq_id", 0)
                        end = region.get("end_seq_id", 0)
                        if start and end:
                            unobserved.append((start, end, end - start + 1))

        return unobserved

    except Exception as e:
        return []


def main():
    print("=" * 70)
    print("STRUCTURAL VALIDATION: D/E REGIONS vs DISORDER")
    print("=" * 70)

    results = {"high_de": [], "low_de": []}

    # high D/E families
    print("\n[HIGH D/E FAMILIES]")
    for pfam in HIGH_DE_FAMILIES:
        print(f"\n{pfam}:")
        pdb_ids = search_pdb_by_pfam(pfam, limit=3)
        print(f"  Found {len(pdb_ids)} structures: {pdb_ids}")

        for pdb_id in pdb_ids:
            quality = get_structure_quality(pdb_id)
            if quality:
                unobserved = get_unobserved_regions(pdb_id)
                quality["unobserved_regions"] = unobserved
                quality["n_unobserved_regions"] = len(unobserved)
                quality["total_unobserved"] = sum(r[2] for r in unobserved)
                results["high_de"].append(quality)

                res_str = f"{quality['resolution']:.2f}Å" if quality['resolution'] else "N/A"
                print(f"    {pdb_id}: res={res_str}, unobserved_regions={len(unobserved)}")
                if unobserved:
                    for start, end, length in unobserved[:3]:
                        print(f"      - residues {start}-{end} ({length} aa)")

    # low D/E families (controls)
    print("\n[LOW D/E FAMILIES - CONTROLS]")
    for pfam in LOW_DE_FAMILIES:
        print(f"\n{pfam}:")
        pdb_ids = search_pdb_by_pfam(pfam, limit=3)
        print(f"  Found {len(pdb_ids)} structures: {pdb_ids}")

        for pdb_id in pdb_ids:
            quality = get_structure_quality(pdb_id)
            if quality:
                unobserved = get_unobserved_regions(pdb_id)
                quality["unobserved_regions"] = unobserved
                quality["n_unobserved_regions"] = len(unobserved)
                quality["total_unobserved"] = sum(r[2] for r in unobserved)
                results["low_de"].append(quality)

                res_str = f"{quality['resolution']:.2f}Å" if quality['resolution'] else "N/A"
                print(f"    {pdb_id}: res={res_str}, unobserved_regions={len(unobserved)}")

    # summary comparison
    print("\n" + "=" * 70)
    print("SUMMARY COMPARISON")
    print("=" * 70)

    if results["high_de"]:
        high_unobs = [r["n_unobserved_regions"] for r in results["high_de"]]
        high_total = [r["total_unobserved"] for r in results["high_de"]]
        print(f"\nHigh D/E families ({len(results['high_de'])} structures):")
        print(f"  Mean unobserved regions: {sum(high_unobs)/len(high_unobs):.1f}")
        print(f"  Mean total unobserved residues: {sum(high_total)/len(high_total):.1f}")

    if results["low_de"]:
        low_unobs = [r["n_unobserved_regions"] for r in results["low_de"]]
        low_total = [r["total_unobserved"] for r in results["low_de"]]
        print(f"\nLow D/E families ({len(results['low_de'])} structures):")
        print(f"  Mean unobserved regions: {sum(low_unobs)/len(low_unobs):.1f}")
        print(f"  Mean total unobserved residues: {sum(low_total)/len(low_total):.1f}")

    # interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    if results["high_de"] and results["low_de"]:
        high_mean = sum(high_unobs) / len(high_unobs)
        low_mean = sum(low_unobs) / len(low_unobs)

        if high_mean > low_mean * 1.5:
            print("\n[+] HIGH D/E families have MORE unobserved (disordered) regions")
            print("    → VALIDATES hypothesis: D/E-rich regions are structurally disordered")
        elif high_mean < low_mean * 0.67:
            print("\n[-] HIGH D/E families have FEWER unobserved regions")
            print("    → Does NOT support hypothesis")
        else:
            print("\n[~] No clear difference in disorder between groups")
            print("    → Inconclusive - may need larger sample or finer analysis")


if __name__ == "__main__":
    main()
