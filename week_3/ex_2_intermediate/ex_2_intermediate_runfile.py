from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, molecular_weight
import csv


# Base class
class Sequence:
    def __init__(self, name, data):
        self.name = name
        self.data = data

    def get_name(self):
        return self.name

    def get_data(self):
        return self.data



class DnaSequence(Sequence):
    def __init__(self, name, data):
        super().__init__(name, data)

    def gc_content(self):

        return gc_fraction(self.data) * 100



class RnaSequence(Sequence):
    def __init__(self, name, data):
        super().__init__(name, data)

    def has_start_codon(self):

        return self.data.startswith("AUG")



class FastaParser:
    def __init__(self, fasta_file, output_csv):
        self.fasta_file = fasta_file
        self.output_csv = output_csv

    def parse_file(self):
        sequences = []


        with open(self.fasta_file, "r") as file:
            for seq_record in SeqIO.parse(file, "fasta"):
                seq_name = seq_record.id
                seq_data = seq_record.seq


                if 'U' in seq_data:
                    sequence = RnaSequence(seq_name, seq_data)
                else:
                    sequence = DnaSequence(seq_name, seq_data)

                sequences.append(sequence)


        with open(self.output_csv, "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["Sequence Name", "Description", "Length", "GC Content (%)", "Molecular Weight", "Type",
                                 "Start Codon Found"])

            for seq in sequences:
                seq_name = seq.get_name()
                seq_data = seq.get_data()
                seq_length = len(seq_data)
                seq_type = "RNA" if isinstance(seq, RnaSequence) else "DNA"
                start_codon = "Yes" if isinstance(seq, RnaSequence) and seq.has_start_codon() else "No"


                if isinstance(seq, DnaSequence):
                    gc_content = seq.gc_content()
                    mol_weight = molecular_weight(seq_data)
                else:

                    mol_weight = molecular_weight(seq_data.replace("U", "T"))
                    gc_content = ""


                print(
                    f"{seq_name}, {seq_type}, Length: {seq_length}, GC Content: {gc_content}, Molecular Weight: {mol_weight}, Start Codon: {start_codon}")

                # Write to CSV
                csv_writer.writerow(
                    [seq_name, "", seq_length, f"{gc_content:.2f}" if gc_content else "", f"{mol_weight:.2f}", seq_type,
                     start_codon])


if __name__ == "__main__":
    fasta_file = "sequences_intermediate.fasta"
    output_csv = "ex_2_intermediate_submit.csv"
    parser = FastaParser(fasta_file, output_csv)
    parser.parse_file()
    print(f"Data written to {output_csv}")
