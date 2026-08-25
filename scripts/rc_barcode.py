import sys


comp = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N",
}


for line in sys.stdin:
    seq = line.strip()

    rc_seq = "".join(
        comp.get(base, base)
        for base in reversed(seq)
    )

    print(rc_seq)
