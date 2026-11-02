import random, sys
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from structure import dinuc_shuffle

def test_dinuc_shuffle_preserves_counts_and_ends():
    s = "GCCGAGUAGGAUCGAGGGUACAGAAUCGGAUUCCAGGCUAGGAUCC"
    rng = random.Random(1)
    for _ in range(20):
        t = dinuc_shuffle(s, rng)
        assert Counter(zip(t, t[1:])) == Counter(zip(s, s[1:]))
        assert t[0] == s[0] and t[-1] == s[-1]
