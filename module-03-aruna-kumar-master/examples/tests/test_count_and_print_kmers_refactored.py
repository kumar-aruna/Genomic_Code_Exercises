#!user/bin/env python3
"""Test behavior of count_and_print_kmers_refactored.py"""

from count_and_print_kmers_refactored import count_kmers # import the function to test

def test_simple_sequence():
    """Count 6-mers in a simple 7 nt sequence."""
    assert count_kmers("ATGTATG", 6) == {'ATGTAT':1, 'TGTATG':1}, "expect a two sequence dictionary"

def test_hanging_kmer():
    """Count 4-mers in a 3 nt sequence - expect no results."""
    assert count_kmers("ATG", 4) == {}, "expect an empty dictionary, sequence length < kmer length"

def test_should_fail():
	"""A test that should fail as an example."""
	assert count_kmers("ATG", 3) == {"TTG":1}, "expect the test to fail"