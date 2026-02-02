from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

# RAA2: Hydrophobic -> A, Polar/Charged -> S
RAA2_MAP_CANONICAL = {
    'A': 'A', 'M': 'A', 'W': 'A', 'L': 'A', 'Y': 'A', 'C': 'A', 'F': 'A', 'I': 'A', 'V': 'A',
    'P': 'S', 'G': 'S', 'H': 'S', 'T': 'S', 'S': 'S', 'D': 'S', 'E': 'S', 'K': 'S', 'N': 'S', 'Q': 'S', 'R': 'S'
}


def reduce_to_raa2(seq: str) -> str:
    return ''.join(RAA2_MAP_CANONICAL.get(aa, 'X') for aa in seq)


def process_fasta(input_path: Path, output_path: Path):
    all_records = []
    for record in SeqIO.parse(input_path, "fasta"):
        reduced_seq = reduce_to_raa2(str(record.seq))
        new_record = SeqRecord(Seq(reduced_seq), id=record.id, description="RAA2-reduced")
        all_records.append(new_record)

    with open(output_path, "w") as f:
        SeqIO.write(all_records, f, "fasta")
    print(f"✅ Wrote {len(all_records)} RAA2-reduced sequences to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Reduce FASTA sequences to RAA2 (canonical A/S) format for PLM fine-tuning")
    parser.add_argument("--input", type=str, required=True, help="Input UniRef50 FASTA file")
    parser.add_argument("--output", type=str, default="raa2_uniref50.fasta", help="Output RAA2 FASTA file")
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()

    process_fasta(input_path, output_path)


if __name__ == "__main__":
    main()

