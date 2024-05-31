import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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


def parse_varscan_output(input_file: str, output_file: str = None) -> None:
    if output_file is None:
        name = os.path.basename(input_file).split('.')[0]
        output_file = name + '.tsv'

    current_line = []
    new_line = ["CHROM\tPOS\tNote\tREF\tALT\tGT\tGQ\tSDP\tDP\tRD\tAD\tFREQ\tPVAL\tRBQ\tABQ\tRDF\tRDR\tADF\tADR"]

    with open(input_file, mode='r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                continue
            else:
                line = line.split('\t')
                info = line[-1].split(':')
                current_line.extend(line[0:5])
                current_line.extend(info)
                new_line.append('\t'.join(current_line))
                current_line = []

    with open(output_file, mode='w') as outfile:
        outfile.write('\n'.join(new_line))


def filter_cravat_xls(folder: str, total_reads: int = 50):
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
                               'Variant AF', 'rsID', 'Global AF',
                               'African AF', 'Ashkenazi Jewish AF', 'East Asian AF', 'Finnish AF', 'Latino AF',
                               'Non-Fin Eur AF', 'Other AF', 'South Asian AF']
            df = df[columns_to_keep]
            filtered_df = df.query(f"`Total reads` > {total_reads}")
            sorted_df = filtered_df[filtered_df['Global AF'].isna() | (filtered_df['Global AF'] < 0.01)]

            dest_file = open(os.path.join(path_to_filtered, f'{name}.xlsx'), 'wb')
            sorted_df.to_excel(dest_file, index=False, engine='openpyxl')
            print("Атлищна! отфильтрованные данные в ", dest_file)