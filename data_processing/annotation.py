import os
import sys
import subprocess
import argparse

def annotate_variants(input_dir, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # List to store base names of input files
    names = []

    # Check if input directory is empty
    if not os.listdir(input_dir):
        print("Empty input directory")
        sys.exit(1)
    else:
        # Collect names of files to be annotated
        for infile in os.listdir(input_dir):
            if infile.endswith(".vt.vcf.gz"):
                name = infile.replace(".vt.vcf.gz", "")
                names.append(name)

    print("Files to be annotated:", names)

    # List of annotators to be used
    annotators = [
        # "brca1_func_assay",
        "cancer_genome_interpreter",
        "cancer_hotspots",
        # "civic",
        # "civic_gene",
        "clingen",
        "clingen_allele_registry",
        "clinpred",
        "clinvar",
        "clinvar_acmg",
        "cosmic",
        "dbsnp",
        "fathmm",
        "gnomad3",
        "metalr",
        "metasvm",
        "mutation_assessor",
        "mutationtaster",
        "mutpred1",
        "mutpred_indel",
        "oncokb",
        "omim",
        "pharmgkb",
        "phylop",
        "polyphen2",
        "provean",
        "pseudogene",
        "sift",
        "spliceai",
        "uniprot",
        "uniprot_domain"
    ]

    # Annotate each file
    for name in names:
        print(f"Annotating {name}...")
        cravat_dir = os.path.join(output_dir, "annotated")
        os.makedirs(cravat_dir, exist_ok=True)

        cravat_command = [
            "oc", "run",
            os.path.join(input_dir, f"{name}.vt.vcf.gz"),
            "-a", *annotators,
            "-d", cravat_dir,
            "-l", "hg38",
            "-t", "excel"
        ]

        # Run CRAVAT annotation command
        subprocess.run(cravat_command, check=True, text=True)

def main():
    parser = argparse.ArgumentParser(description="Annotation of germline variants")
    parser.add_argument("--input_dir", type=str, required=True, help="Absolute path to input folder")
    parser.add_argument("--output_dir", type=str, default="./annotations", help="Absolute path to output folder")
    args = parser.parse_args()

    annotate_variants(args.input_dir, args.output_dir)

if __name__ == "__main__":
    main()
