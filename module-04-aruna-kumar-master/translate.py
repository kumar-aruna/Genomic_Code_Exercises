import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO





def main():
    for record in SeqIO.parse(infile, "fasta"):
        if re.match("^\d{1}\D*$", record.id):

            print(" ")

    def find_first_orf(rna):
        # Takes an RNA sequence and returns the first complete ORF as a Seq object

        orf = re.search('AUG([AUGC]{3})+?(UAA|UAG|UGA)', str(rna)).group()

        # Regex return string so convert string to Sequence
        org_seq = Seq(orf)

        return org_seq

    def translate_first_orf(dna):
        """takes a DNA sequence, transcribes it into RNA,
        finds the first ORF, translates said ORF into a protein, and returns that protein"""

        # Transcribe a DNA sequence into RNA and return the RNA sequence as a new Seq object

        try:
            rna = dna.transcribe()
            print("RNA:", rna)

            orf = find_first_orf(rna)
            print("orf: ", orf)

            protein = orf.translate()
            print(protein)
        except AttributeError:
            orf = ""
            return Seq(orf)


    dna = Seq("ATGGCCATTGTAATGGGCCCCGCTGAAAGGGTGCCCGATAG")
    translate_first_orf(dna)


if __name__ == "__main__":
    infile = str(input("Enter file name: "))
    main()

# <FASTA record ID>:\t<translated first ORF>