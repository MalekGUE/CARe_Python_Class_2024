from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, molecular_weight
import csv


class FastaParser:
    def __init__(self, fasta_file, output_csv):
        self.fasta_file = fasta_file
        self.output_csv = output_csv

    def parse_file(self):
        with open(self.fasta_file, "r") as file:
            sequences = SeqIO.parse(file, "fasta")

            with open(self.output_csv, "w", newline="") as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(["Sequence Name", "Description", "Length", "GC Content (%)", "Molecular Weight"])

                for seq_record in sequences:
                    seq_name = seq_record.id
                    seq_description = seq_record.description
                    seq_length = len(seq_record.seq)
                    gc_content = gc_fraction(seq_record.seq) * 100
                    mol_weight = molecular_weight(seq_record.seq)

                    print(f"{seq_name}, {seq_description}, {seq_length}, {gc_content:.2f}, {mol_weight:.2f}")

                    csv_writer.writerow(
                        [seq_name, seq_description, seq_length, f"{gc_content:.2f}", f"{mol_weight:.2f}"])


if __name__ == "__main__":
    fasta_file = "sequences_beginner.fasta"
    output_csv = "ex_2_beginner_submit.csv"
    parser = FastaParser(fasta_file, output_csv)
    parser.parse_file()
    print(f"Data written to {output_csv}")

