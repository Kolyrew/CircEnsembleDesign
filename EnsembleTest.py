import subprocess

"""
[IRES] = TTTATAATTTCTTCTTCCAGAA
[Spicer]: UC
[Kozak]: GCCACCAUG
[Junction]: UUAAGCUA
"""

prefix = "TTTATAATTTCTTCTTCCAGAAUCGCCACCAUG"  # IRES=22
suffix = "UUAAGCUA"
iresEND = 22

prefix = "TTTATAATTTCTTCTTCCAGAA"  # IRES=22
suffix = ""
iresEND = 22


beam = 500
iters = 10
lr = 0.03
epsilon = 0.2
L = 2.5


def exec(props, input, cwd="./"):
    proc = subprocess.run(
        props, input=input + "\n", cwd=cwd,
        capture_output=True, text=True, check=True)
    out_lines = proc.stdout.splitlines()
    err_lines = proc.stderr.splitlines()
    return "\n".join(out_lines + err_lines).split("\n")


def find(lines, Q):
    for l in lines:
        if Q in l:
            return l.split(":")[1].strip()
    return None


def colored_ires(seq, start, color=31):
    return f"\033[1;{color}m{seq[0:start]}\033[m{seq[start+1:]}"


def count_cross_pairs(s: str, k: int):
    stack = []
    cross = 0
    for i, c in enumerate(s):
        if c == '(':
            stack.append(i)
        elif c == ')':
            j = stack.pop()
            cross += (j <= k) != (i <= k)
    return cross


def EnsembleWrapper(seq, init, iresEND):
    params = ["./bin/EnsembleDesign"]
    params += list(map(str, [beam, iters, lr, epsilon, 1, init]))
    params += ["--prefix", prefix, "--suffix",
               suffix, "--iresEND", str(iresEND)]
    res = exec(params, seq, "./")
    mRNA = find(res, "Final mRNA")

    res = exec(["RNAfold", "--noPS"], mRNA, "./")
    struct = res[1].split(" ")[0]

    crp = count_cross_pairs(struct, 22)

    print(f" - mRNA:  {colored_ires(mRNA, 22)}")
    print(f"          {colored_ires(res[1], 22)}")
    print(f"   CROSS: {crp}")

    return mRNA, struct


def wrapper(seq):
    res = exec([
        "./bin/LinearDesign_2D", f"{L:.2f}", "0", "codon_usage_freq_table_human.csv"
    ], seq, "./tools/LinearDesign/")
    INIT_MRNA_SOLUTION = find(res, "mRNA sequence")

    print(f" - SEQ:   {seq}")
    print(f" - Init:  {INIT_MRNA_SOLUTION}")
    mRNA, struct = EnsembleWrapper(seq, INIT_MRNA_SOLUTION, 0)
    mRNA, struct = EnsembleWrapper(seq, INIT_MRNA_SOLUTION, 22)


with open('data/uniprot.fasta', 'r') as fa:
    for example in ("".join(fa.readlines())).split("\n>"):
        if not example.strip():
            continue
        name, prot = example.split("\n")
        seq = prot[0:30]
        print("-" * 60)
        print("Name: ", name)
        wrapper(seq)
