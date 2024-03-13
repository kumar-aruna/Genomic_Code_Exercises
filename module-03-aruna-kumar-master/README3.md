
# Identifying the Sequence

In this program, we are refactoring the code of identify sequence file and pytesting the given file that contains Nucleic acid or Aminoacid or lower case nucleic acid or rna sequences.


**Created by: Aruna Kumar**

**Created Date: Feb-13-2023**


```
``$python3 identify_sequence_refactored.py sequence1.txt sequence2.txt sequence3.txt

```


## 1. Finding a Nucleotide sequence/Aminoacid sequence/rna sequence/lowe case sequence

```
sequence1.txt
Output:
The length of the sequence is : 243536 
Filename is sequence1.txt and sequence type is nucleic acid


Process finished with exit code 0
```


```
sequence2.txt
Output:
The length of the sequence is : 86 
Filename is sequence2.txt and sequence type is lower nucleic acid
Process finished with exit code 0

```
```
sequence3.txt

sequence3 Output
The length of the sequence is : 169 
Filename is sequence3.txt and sequence type is amino acid

Process finished with exit code 0
```

###Testing  behavior of identify sequence refactored file

```
pytest 
```
 ##output
 
```
 assert identify_sequence("ZZXy43") == "No sequence", "Invalid sequence"
E       AssertionError: Invalid sequence
E       assert '' == 'No sequence'
E         - No sequence

tests/test_identify_sequence_refactored.py:32: AssertionError
=========================================================================== short test summary info ============================================================================
FAILED tests/test_identify_sequence_refactored.py::test_nonsequence - AssertionError: Invalid sequence
========================================================================= 1 failed, 4 passed in 0.07s ==========================================================================
```




