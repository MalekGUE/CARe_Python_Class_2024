from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import molecular_weight, gc_fraction
import csv


# Base class for sequences
class Sequence:
    def __init__(self, name, description, data):
        self.name = name
        self.description = description
        self.data = data

    def calculate_length(self):
        return len(self.data)

    def gc_content(self):
        return gc_fraction(self.data)

    def transcription(self):
        raise NotImplementedError("This method should be implemented by subclasses.")

    def translation(self):
        raise NotImplementedError("This method should be implemented by subclasses.")

    def molecular_weight(self):
        raise NotImplementedError("This method should be implemented by subclasses.")


# DNA sequence class
class DnaSequence(Sequence):
    def transcription(self):
        return str(Seq(self.data).transcribe())

    def translation(self):
        return str(Seq(self.data).translate())


# RNA sequence class
class RnaSequence(Sequence):
    def transcription(self):
        return self.data  # RNA is already transcribed

    def translation(self):
        return str(Seq(self.data).translate())


# Protein sequence class
class ProteinSequence(Sequence):
    def translation(self):
        return self.data  # Protein is already translated

    def molecular_weight(self):
        return molecular_weight(self.data)


# FastaParser class
class FastaParser:
    def __init__(self, fasta_file, output_csv):
        self.fasta_file = fasta_file
        self.output_csv = output_csv

    def parse_file(self):
        sequences = []

        # Read and parse the sequences in the FASTA file
        with open(self.fasta_file, "r") as file:
            for seq_record in SeqIO.parse(file, "fasta"):
                seq_name = seq_record.id
                seq_description = seq_record.description
                seq_data = str(seq_record.seq)

                if "U" in seq_data:  # RNA sequence
                    sequence = RnaSequence(seq_name, seq_description, seq_data)
                elif "T" in seq_data:  # DNA sequence
                    sequence = DnaSequence(seq_name, seq_description, seq_data)
                else:  # Protein sequence
                    sequence = ProteinSequence(seq_name, seq_description, seq_data)

                sequences.append(sequence)

        # Write the results to the CSV file
        with open(self.output_csv, "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow([
                "Sequence Name", "Description", "Sequence Data", "Length",
                "GC Content (%)", "RNA Transcript", "Protein Translation", "Molecular Weight"
            ])

            for seq in sequences:
                seq_name = seq.name
                seq_description = seq.description
                seq_data = seq.data
                seq_length = seq.calculate_length()
                gc_content = seq.gc_content()
                rna_transcript = seq.transcription()
                protein_translation = seq.translation() if isinstance(seq, (DnaSequence, RnaSequence)) else ""
                mol_weight = seq.molecular_weight() if isinstance(seq, ProteinSequence) else None

                # Print sequence details
                print(f"{seq_name}: {seq_description}")
                print(f"Sequence: {seq_data}")
                print(f"Length: {seq_length}, GC Content: {gc_content:.2f}%")
                print(f"RNA Transcript: {rna_transcript}")
                print(f"Protein Translation: {protein_translation}")
                if mol_weight is not None:
                    print(f"Molecular Weight: {mol_weight:.2f}")
                else:
                    print(f"Molecular Weight: N/A")
                print("-" * 50)

                # Write to CSV
                csv_writer.writerow([
                    seq_name, seq_description, seq_data, seq_length, f"{gc_content:.2f}",
                    rna_transcript, protein_translation, f"{mol_weight:.2f}" if mol_weight is not None else "N/A"
                ])


if __name__ == "__main__":
    fasta_file = "sequences_advanced.fasta"  # Ensure this file exists in the directory
    output_csv = "ex_2_advanced_submission.csv"
    parser = FastaParser(fasta_file, output_csv)
    parser.parse_file()
    print(f"Data written to {output_csv}")

