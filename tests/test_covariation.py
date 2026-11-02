import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from covariation import pairs_from_ss, mi

def test_pairs_from_ss():
    assert pairs_from_ss("<<..>>") == [(0, 5), (1, 4)]

def test_mi_perfect_and_independent():
    a = list("AAAAUUUUGGGGCCCC"); b = list("UUUUAAAACCCCGGGG")
    assert mi(a, b) > 1.9
    c = list("ACGUACGUACGUACGU"); d = list("AAAACCCCGGGGUUUU")
    assert mi(c, d) < 0.01
