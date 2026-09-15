import argparse
import os
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input FASTA file(s)", nargs="+", type=str)
    parser.add_argument("output", help="Output TXT file", type=str)

    args = parser.parse_args()

    step00.check_suffixes(args.input, {".fasta"})

    for input_file in tqdm.tqdm(args.input):
        total_length = 0
        small_length = 0
        with open(input_file, "r") as file:
            for line in tqdm.tqdm(file.readlines(), position=1, leave=False):
                line = line.strip()

                if line.startswith(">"):
                    continue

                total_length += sum(1 for char in line if char.isalpha())
                small_length += sum(1 for char in line if char.islower())

        print(os.path.basename(input_file), ":", f"{small_length} / {total_length}", "=", f"{small_length / total_length * 100:.2f}")
