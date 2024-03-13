import random


"""
This Function Identifies the ORF
"""

# Global Variables
QUALITY_SCORE_CUTOFF = 7.25
MIN_ORF_LENGTH = 60  # 60 Nucleotides coding for amino acids as min ORF length
base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3} # letter to number
# matching for index
MOTIF_MATRIX = [[.5,.5,.5,.5,0,0,0,0,0,2,-99,-99,.5], #A
[0,0,0,0,0,0,0,0,0,-99,2,-99,0],  # T
[0,0,0,0,0,0,0,0,0,-99,-99,-99,0],  # C
[.5,.5,.5,.5,0,0,0,0,0,.5,-99,2,0]]  # G


# def identifyORFs(): #main function: performs the goals of the program.
# Finding the most likely coding region
def identify_ORFs(input_file, output_file):
    pass


# def scanSeq(): #searches a single sequence for possible orfs. Returns
# start positions, lengths and sequences of ORFs
def scanSeq(sequence):
    pass


# def scoreMotif(): #scores a sequence the same length as the motif[13
# bases]. Returns motif score.
def scoreMotif():
    pass




l = 20 #read length
k = 4
bases = ["A","T","C","G"]
def generateReads(l,k): #prints k reads of length l
    for i in range(k):
        seq = ""
        for j in range(l):
            seq += bases[random.randint(0,3)]
        print(">sequence " + str(i+1))
        print(seq)

generateReads(l,k)
