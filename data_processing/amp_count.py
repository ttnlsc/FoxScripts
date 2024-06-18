import os
import sys
import subprocess
import argparse

def main(args):
    input_dir = args.input_dir
    log_dir = args.logs
    os.makedirs(log_dir, exist_ok=True)

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    if not os.listdir(input_dir):
        print("Empty input directory")
        sys.exit(1)
    else:
        for infile in os.listdir(input_dir):
            if infile.endswith("bam"):
                name = infile.split('.')[0]
                infile_path = os.path.join(input_dir, infile)
                outfile = os.path.join(output_dir, f"{name}.txt")
                logfile = os.path.join(log_dir, f"{name}.log")

                count_command = [
                    "bedtools", "multicov", "-bams", infile_path, "-bed", args.bed,
                ]

                with open(outfile, "w") as of, open(logfile, "w") as lf:
                    subprocess.run(count_command, stdout=of, stderr=lf, text=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count reads in BAM files using bedtools multicov")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--input_dir", type=str, required=True, help="Absolute path to input folder")
    parser.add_argument("--output_dir", type=str, required=True, help="Absolute path to output folder")
    parser.add_argument("--bed", type=str, required=True, help="Bed file")
    parser.add_argument("--logs", type=str, required=True, help="Absolute path to logs folder")
    args = parser.parse_args()

    main(args)
