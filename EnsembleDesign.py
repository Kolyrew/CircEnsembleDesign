import os
import re
import sys
import subprocess
import threading
import argparse
import math
import shutil
from queue import Queue
from multiprocessing import Value, Lock

prog_path = "./bin/EnsembleDesign"


def read_fasta(file_path):
    records = []
    with open(file_path, "r") as file:
        seq_id = None
        seq_lines = []
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    records.append((seq_id, "".join(seq_lines)))
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)

        if seq_id is not None:
            records.append((seq_id, "".join(seq_lines)))
    return records

def get_mfe_solutoin(seq, ires=""):
    try:
        if ires:
            # For now, use fallback when IRES is present
            codon_map = {
                'A': 'GCU', 'R': 'CGU', 'N': 'AAU', 'D': 'GAU', 'C': 'UGU',
                'Q': 'CAA', 'E': 'GAA', 'G': 'GGU', 'H': 'CAU', 'I': 'AUU',
                'L': 'CUA', 'K': 'AAA', 'M': 'AUG', 'F': 'UUU', 'P': 'CCU',
                'S': 'UCU', 'T': 'ACU', 'W': 'UGG', 'Y': 'UAU', 'V': 'GUU'
            }
            rna = []
            for aa in seq.strip():
                rna.append(codon_map.get(aa, 'AUG'))
            return ''.join(rna)
        else:
            full_command = f"cd ./tools/LinearDesign && echo {seq} | ./lineardesign"
            result = subprocess.run(full_command, shell=True, capture_output=True, text=True)
            mfe_solution = result.stdout.strip()
            return mfe_solution.split("\n")[-3].strip().split(" ")[-1]
    except Exception:
        # Fallback: generate a simple codon-mapped RNA (first-choice codons)
        codon_map = {
            'A': 'GCU', 'R': 'CGU', 'N': 'AAU', 'D': 'GAU', 'C': 'UGU',
            'Q': 'CAA', 'E': 'GAA', 'G': 'GGU', 'H': 'CAU', 'I': 'AUU',
            'L': 'CUA', 'K': 'AAA', 'M': 'AUG', 'F': 'UUU', 'P': 'CCU',
            'S': 'UCU', 'T': 'ACU', 'W': 'UGG', 'Y': 'UAU', 'V': 'GUU'
        }
        rna = []
        for aa in seq.strip():
            rna.append(codon_map.get(aa, 'AUG'))
        return ''.join(rna)

def eval_partition(seq):
    full_command = f"cd ./tools/LinearPartition && echo {seq} | ./linearpartition -V -d0 -b0"
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True)
    partition_result = result.stderr.strip()
    return float(partition_result.split(' ')[-2].strip())


def count_cross_pairs_probabilistic(seq, ires_len, prob_threshold=0.01):
    """Count cross-pairs between IRES and all rest (ORF + suffix) using P(i,j) > threshold. 1-based: IRES = [1, ires_len]."""
    if ires_len == 0:
        return 0, [], 0.0
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        full_command = f"cd ./tools/LinearPartition && echo {seq} | ./linearpartition -V -c {prob_threshold} -r {tmp_path}"
        subprocess.run(full_command, shell=True, capture_output=True, text=True)
        cross_pairs = []
        sum_prob = 0.0
        with open(tmp_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) != 3:
                    continue
                i, j, prob = int(parts[0]), int(parts[1]), float(parts[2])
                i_in_ires = (1 <= i <= ires_len)
                j_in_ires = (1 <= j <= ires_len)
                if i_in_ires != j_in_ires:
                    cross_pairs.append((i, j, prob))
                    sum_prob += prob
        return len(cross_pairs), cross_pairs, sum_prob
    except Exception:
        return 0, [], 0.0
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

def process_run_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    seq = re.search(r'Final mRNA sequence: (.*)', lines[-1]).group(1)
    efe = eval_partition(seq)
    return seq, efe

def execute_mrna_design(protein, output_dir, run_number, args, progress_tracker, lock):

    mfe_solutoin = get_mfe_solutoin(protein, args.ires)
    command = f"echo {protein} | {prog_path} {args.beam_size} {args.num_iters} {args.lr} {args.epsilon} {run_number} {mfe_solutoin} {args.ires} {args.ires_orf_lambda} {args.cross_pair_prob_threshold} {args.max_cross_pairs}"

    output_path = os.path.join(output_dir, f"run_{run_number}.txt")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as output_file:
        subprocess.run(command, shell=True, stdout=output_file, stderr=subprocess.STDOUT)

    with lock:
        progress_tracker.value += 1
        sys.stderr.write(f"Completed: {progress_tracker.value} runs\n")
        sys.stderr.flush()

def worker(task_queue, progress_tracker, lock):
    while True:
        task = task_queue.get()
        if task is None:
            break
        execute_mrna_design(*task, progress_tracker=progress_tracker, lock=lock)
        task_queue.task_done()

def run_mrna_design(args):
    records = read_fasta(args.fasta)
    num_runs = args.num_runs
    num_threads = args.num_threads

    task_queue = Queue()
    progress_tracker = Value('i', 0)
    lock = Lock()

    # Initialize progress tracker
    total_runs = len(records) * num_runs
    sys.stderr.write(f"{len(records)} input sequences.\n{num_runs} runs for each sequence.\nTotal {total_runs} runs started.\n")

    # Create and start threads
    threads = []
    for _ in range(num_threads):
        t = threading.Thread(target=worker, args=(task_queue, progress_tracker, lock))
        t.start()
        threads.append(t)

    # Clean log folders
    for seq_id, _ in records:
        log_dir = os.path.join(args.output_dir, seq_id)
        if os.path.exists(log_dir):
            shutil.rmtree(log_dir)
            
    # Enqueue tasks
    for run_number in range(1, num_runs + 1):
        for seq_id, protein in records:
            log_dir = os.path.join(args.output_dir, seq_id)
            task_queue.put((protein, log_dir, run_number, args))

    # Block until all tasks are done
    task_queue.join()

    # Stop workers
    for _ in range(num_threads):
        task_queue.put(None)
    for t in threads:
        t.join()

    sys.stderr.write(f"All {total_runs} runs completed.\n")

    ires_len = len(args.ires)
    max_cp = args.max_cross_pairs

    width = int(math.log10(num_runs)) + 1;
    for seq_id, protein in records:
        log_dir = os.path.join(args.output_dir, seq_id)
        results = []
        for run_file in sorted(os.listdir(log_dir)):
            if not run_file.endswith(".txt"):
                continue
            run_id = int(run_file[:-4].split("_")[-1])
            seq, efe = process_run_file(os.path.join(log_dir, run_file))
            results.append((run_id, seq, efe))

        prob_thresh = args.cross_pair_prob_threshold
        logs = []
        candidates = []
        for run_id, seq, efe in sorted(results, key=lambda x: x[0]):
            cp_count = 0
            cp_list = []
            sum_prob = 0.0
            if ires_len > 0:
                cp_count, cp_list, sum_prob = count_cross_pairs_probabilistic(seq, ires_len, prob_thresh)
            status = "OK" if cp_count <= max_cp else "REJECTED"
            logs.append(f"[{run_id:>{width}}] {seq} | EFE: {efe} kcal/mol | cross-pairs(P>{prob_thresh}): {cp_count} sum_P={sum_prob:.4f} ({status})")
            candidates.append((run_id, seq, efe, cp_count, cp_list))

        sorted_by_efe = sorted(candidates, key=lambda x: x[2])
                                                                                                                                                                                
        best_seq, best_efe, best_cp = None, 0, 0
        for run_id, seq, efe, cp_count, cp_list in sorted_by_efe:
            if cp_count <= max_cp:
                best_seq, best_efe, best_cp = seq, efe, cp_count
                break

        if best_seq is None and sorted_by_efe:
            fallback = min(sorted_by_efe, key=lambda x: (x[3], x[2]))
            best_seq, best_efe, best_cp = fallback[1], fallback[2], fallback[3]
            logs.append(f"[WARNING] No sequence with <= {max_cp} cross-pairs found, using best available ({best_cp} cross-pairs)")

        logs.append(f"[Best] {best_seq} | EFE: {best_efe} kcal/mol | cross-pairs: {best_cp}")

        cp_info = f"|Cross-pairs: {best_cp}" if ires_len > 0 else ""
        print(f">{seq_id}|Ensemble Free Energy: {best_efe} kcal/mol{cp_info}\n{best_seq}")

        result_file = os.path.join(log_dir, "results.txt")
        with open(result_file, 'w') as f:
            f.write("\n".join(logs))

def main():
    parser = argparse.ArgumentParser(description='Run EnsembleDesign on protein fasta file.')
    parser.add_argument('--fasta', type=str, default='./examples.fasta', help='Path to the input protein fasta file (default: ./examples.fasta)')
    parser.add_argument('--output_dir', type=str, default='./outputs', help='Directory to save output files (default: ./outputs)')
    parser.add_argument('--beam_size', type=int, default=200, help='Beam size for beam pruning (default: 200)')
    parser.add_argument('--lr', type=float, default=0.03, help='Learning rate for projected gradient decent (default: 0.03)')
    parser.add_argument('--epsilon', type=float, default=0.5, help='The epsilon paramter of soft-MFE initialization (default: 0.5)')
    parser.add_argument('--num_iters', type=int, default=30, help='Number of optimization steps (default: 30)')
    parser.add_argument('--num_runs', type=int, default=20, help='Number of runs per sample file (default: 20)')
    parser.add_argument('--num_threads', type=int, default=8, help='Number of threads in the thread pool (default: 16)')
    parser.add_argument('--ires', type=str, default='', help='IRES sequence to prepend to the mRNA (default: empty)')
    parser.add_argument('--ires_orf_lambda', type=float, default=1.0, help='Penalty factor for IRES-ORF base pairing (default: 1.0, no penalty; use < 1.0 like 0.1 to discourage IRES-ORF pairs)')
    parser.add_argument('--max_cross_pairs', type=int, default=5, help='Max allowed IRES-ORF cross-pairs (default: 5; sequences with more are rejected)')
    parser.add_argument('--cross_pair_prob_threshold', type=float, default=0.01, help='Base pair probability threshold for cross-pair detection (default: 0.01)')

    args = parser.parse_args()
    run_mrna_design(args)

if __name__ == "__main__":
    main()