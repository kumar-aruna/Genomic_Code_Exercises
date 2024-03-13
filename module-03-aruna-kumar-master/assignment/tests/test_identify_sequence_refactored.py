#!user/bin/env python3
"""Test behavior of identify_sequence_refactored.py"""

from ..identify_sequence_refactored import identify_sequence # import the function to test (make sure you rename your script to "indentify_sequence_refactored.py")


def test_dna_sequence():
    """Identify a DNA sequence"""
    assert identify_sequence("ATGATGA") == "nucleic acid", "expect ATGATGA identifies as nucleic acid"


def test_dna_lowercase_sequence():
    """Identify a DNA sequence in lowercase"""
    assert identify_sequence("atgcat") == "lower nucleic acid", "expect atgcat identifies as lowercase_nucleic_acid"


def test_rna_sequence():
    """Identify an RNA sequence"""
    assert identify_sequence("AUGCAU") == "rna sequence", "expect AUGCAU identifies as rna sequence"


def test_aminoacid_sequence():
    """Identify an amino acid sequence"""
    assert identify_sequence("MSNKA") == "amino acid", "expect MSNKA identifies as amino acid"

# Note: we expect this test to fail right now, but perhaps we can refactor the code again to make it pass in the future!
# For now, let's write the test to have that reminder.


def test_nonsequence():
    """Identify invalid sequence"""
    assert identify_sequence("ZZXy43") == "No sequence", "Invalid sequence"
