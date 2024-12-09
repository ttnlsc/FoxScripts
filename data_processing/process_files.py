import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

AA_NAME_DICT = {'Gly': 'G', 'Leu': 'L', 'Tyr': 'Y',
                'Ser': 'S', 'Glu': 'E', 'Gln': 'Q',
                'Asp': 'D', 'Asn': 'N', 'Phe': 'F',
                'Ala': 'A', 'Lys': 'K', 'Arg': 'R',
                'His': 'H', 'Cys': 'C', 'Val': 'V',
                'Pro': 'P', 'Trp': 'W', 'Ile': 'I',
                'Met': 'M', 'Thr': 'T'}


def read_metrics_to_dataframe(filepath: str):
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('## METRICS CLASS'):
                break
        else:
            raise ValueError(f"File {filepath} does not contain '## METRICS CLASS'")

        try:
            header_line = next(file).strip()
            headers = header_line.split('\t')
            data_line = next(file).strip()
            values = data_line.split('\t')
        except StopIteration:
            raise ValueError(f"File {filepath} does not contain enough lines after '## METRICS CLASS'")
        if len(headers) > len(values):
            values.extend([None] * (len(headers) - len(values)))
    df = pd.DataFrame([values], columns=headers)

    return df


def collect_pcr_metrics(folder: str, outfile_name: str):
    metrics = ['Source.Name\tPF_READS\tPF_BASES\tMEAN_COVER\tUniformity']
    for file in os.listdir(folder):
        if file.endswith('.txt'):
            name = file.split('_')[0]
            path = os.path.join(folder, file)
            df = read_metrics_to_dataframe(path)
            reads = int(df['PF_READS'].iloc[0])
            bases = int(df['PF_BASES'].iloc[0])
            mean_cover = round(float(df['MEAN_TARGET_COVERAGE'].iloc[0]), 1)
            if mean_cover <= 100:
                uniformity = float(df['PCT_TARGET_BASES_10X'].iloc[0])
            elif 100 < mean_cover <= 150:
                uniformity = float(df['PCT_TARGET_BASES_20X'].iloc[0])
            elif 150 < mean_cover <= 200:
                uniformity = float(df['PCT_TARGET_BASES_30X'].iloc[0])
            elif 200 < mean_cover <= 250:
                uniformity = float(df['PCT_TARGET_BASES_40X'].iloc[0])
            elif 250 < mean_cover <= 500:
                uniformity = float(df['PCT_TARGET_BASES_50X'].iloc[0])
            elif 500 < mean_cover <= 1250:
                uniformity = float(df['PCT_TARGET_BASES_100X'].iloc[0])
            elif 1250 < mean_cover <= 2500:
                uniformity = float(df['PCT_TARGET_BASES_250X'].iloc[0])
            elif 2500 < mean_cover <= 5000:
                uniformity = float(df['PCT_TARGET_BASES_500X'].iloc[0])
            elif 5000 < mean_cover <= 12500:
                uniformity = float(df['PCT_TARGET_BASES_1000X'].iloc[0])
            elif 12500 < mean_cover:
                uniformity = float(df['PCT_TARGET_BASES_2500X'].iloc[0])
            uniformity_pct = round((uniformity * 100), 2)

            metrics.append(f'{name}\t{reads}\t{bases}\t{mean_cover}\t{uniformity_pct}')

    with open(f"{outfile_name}.txt", "w") as outfile:
        outfile.write("\n".join(metrics))


def draw_cover_plot(folder: str):
    for file in os.listdir(folder):
        name = file.split('.')[0]
        path = os.path.join(folder, file)
        try:
            df = pd.read_csv(path, sep='\t')
        except Exception as e:
            print(f"Failed to read {file}: {e}")
            continue
        plt.figure()
        plt.axes()
        sns.lineplot(data=df['read_count'])
        plt.yscale('log')
        plt.xlim(0)
        plt.ylim(1)
        plt.savefig(f'{name}.png', dpi=300, bbox_inches='tight')
        plt.close()


def parse_varscan_output(folder: str):
    # Ensure the "parsed" directory exists
    os.makedirs(os.path.join(folder, "parsed"), exist_ok=True)
    path_to_parsed = os.path.join(folder, "parsed")

    for input_file in os.listdir(folder):
        if input_file.endswith(".vcf"):
            name = input_file.replace(".vcf", "")  # Extract the file name without extension

            # Initialize a list for the header and new data lines
            new_line = ["Chrom\tPosition\tNote\tREF\tALT\tGT\tGQ\tSDP\tDP\tRD\tAD\tFREQ\tPVAL\tRBQ\tABQ\tRDF\tRDR\tADF\tADR"]

            with open(os.path.join(folder, input_file), mode='r') as file:
                for line in file:
                    line = line.strip()

                    if line.startswith('#'):
                        continue  # Skip header lines

                    # Process each data line
                    line = line.split('\t')
                    info = line[-1].split(':')  # Split the INFO column
                    current_line = line[0:5]  # Extract the first 5 fields
                    current_line.extend(info)  # Append additional information fields
                    new_line.append('\t'.join(current_line))  # Add the processed line to the result list

            # Write the processed data to a new file
            with open(os.path.join(path_to_parsed, f'{name}.txt'), mode='w') as outfile:
                outfile.write('\n'.join(new_line))  # Join the lines with newline character


def filter_cravat_xls(folder: str, total_reads: int = 50, vaf: float = 0.2):
    os.makedirs(os.path.join(folder, "filtered"), exist_ok=True)
    path_to_filtered = os.path.join(folder, "filtered")
    for file in os.listdir(folder):
        if file.endswith("vt.vcf.gz.xlsx"):
            name = file.split('.')[0]
            path = os.path.join(folder, file)
            df = pd.read_excel(path, sheet_name="Variant", header=1)
            columns_to_keep = ['UID', 'Chrom', 'Position', 'Ref Base', 'Alt Base', 'Coding', 'Gene', 'Transcript',
                               'Sequence Ontology', 'Exon Number', 'cDNA change',
                               'Protein Change', 'All Mappings', 'Sample Count', 'Samples', 'ID', 'Clinical Significance',
                               'Disease Names', 'Review Status',
                               'ClinVar ID', 'Phred', 'VCF filter', 'Zygosity', 'Alternate reads', 'Total reads',
                               'Variant AF', 'rsID', 'Global AF']
            df = df[columns_to_keep]
            filtered_df = df.query(
                "`Total reads` > @total_reads and `Variant AF` >= @vaf and (`Global AF`.isna() or `Global AF` < 0.01)"
            )

            dest_file = open(os.path.join(path_to_filtered, f'{name}.xlsx'), 'wb')
            filtered_df.to_excel(dest_file, index=False, engine='openpyxl')
            print("Атлищна! отфильтрованные данные в ", dest_file)


def filter_varscan_xls(folder: str, total_reads: int = 50):
    # Define paths
    path_to_xlsx = os.path.join(folder, "annotated_varscan")
    path_to_txt = os.path.join(folder, "varscan", "parsed")
    path_to_filtered = os.path.join(path_to_xlsx, "filtered")

    # Create the filtered directory if it doesn't exist
    os.makedirs(path_to_filtered, exist_ok=True)

    # Iterate through files in the annotated_varscan folder
    for file in os.listdir(path_to_xlsx):
        if file.endswith(".vcf.xlsx"):
            name = file.split('.vcf')[0]
            columns_to_keep = ['UID', 'Chrom', 'Position', 'Ref Base', 'Alt Base', 'Coding', 'Gene', 'Transcript',
                               'Sequence Ontology', 'Exon Number', 'cDNA change',
                               'Protein Change', 'All Mappings', 'Sample Count', 'Samples', 'ID',
                               'Clinical Significance', 'Disease Names', 'Review Status',
                               'ClinVar ID', 'Total reads', 'rsID', 'Global AF']

            # Load the Excel file and CSV file
            df_1 = pd.read_excel(os.path.join(path_to_xlsx, file), sheet_name="Variant", header=1)
            df1 = df1[columns_to_keep]
            df_2 = pd.read_csv(os.path.join(path_to_txt, f'{name}.txt'), sep='\t', header=0)
            if 'indel' in name and 'Position' in df_2.columns:
                df_2['Position'] = df_2['Position'] + 1

            # Merge the dataframes
            df = df_1.merge(df_2, how='left')  # Use 'how' to specify the type of merge
            df['FREQ'] = df['FREQ'].str.rstrip('%').astype(float)
            filtered_df = df.query(
                "`FREQ` >= 1 and `Total reads` > @total_reads and (`Global AF`.isna() or `Global AF` < 0.01)"
            )

            # Define destination file path
            dest_file_path = os.path.join(path_to_filtered, f'{name}.xlsx')

            # Write the merged dataframe to Excel
            filtered_df.to_excel(dest_file_path, index=False, engine='openpyxl')
            print(f"Очень хорошо! Объединенные данные сохранены в {dest_file_path}")


def convert_three_to_one(aa_string: str) -> str:
    def replace_match(match):
        aa_three_ref = match.group(1)
        position = match.group(2)
        aa_three_alt = match.group(3)

        aa_one_ref = AA_NAME_DICT.get(aa_three_ref.capitalize(), aa_three_ref)
        aa_one_alt = AA_NAME_DICT.get(aa_three_alt.capitalize(), aa_three_alt)

        return f"{aa_one_ref}{position}{aa_one_alt}"

    result = re.sub(r"([A-Za-z]{3})(\d+)([A-Za-z]{3})", replace_match, aa_string)
    return result