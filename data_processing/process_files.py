import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def collect_pcr_metrics(folder: str, outfile_name: str):
    metrics = ['Source.Name\tPF_READS\tPF_BASES\tMEAN_COVER\tUniformity']
    for file in os.listdir(folder):
        if file.endswith('.txt'):
            path = os.path.join(folder, file)
            name = file.split('_')[0]
            df = pd.read_csv(path, sep='\t', header=5)
            row = df.loc[0]
            reads = int(row['PF_READS'])
            bases = int(row['PF_BASES'])
            mean_cover = round(row['MEAN_TARGET_COVERAGE'], 1)
            if mean_cover <= 100:
                uniformity = row['PCT_TARGET_BASES_10X']
            elif 100 < mean_cover <= 150:
                uniformity = row['PCT_TARGET_BASES_20X']
            elif 150 < mean_cover <= 200:
                uniformity = row['PCT_TARGET_BASES_30X']
            elif 200 < mean_cover <= 250:
                uniformity = row['PCT_TARGET_BASES_40X']
            elif 250 < mean_cover <= 500:
                uniformity = row['PCT_TARGET_BASES_50X']
            elif 500 < mean_cover <= 1250:
                uniformity = row['PCT_TARGET_BASES_100X']
            elif 1250 < mean_cover <= 2500:
                uniformity = row['PCT_TARGET_BASES_250X']
            elif 2500 < mean_cover <= 5000:
                uniformity = row['PCT_TARGET_BASES_500X']
            elif 5000 < mean_cover <= 12500:
                uniformity = row['PCT_TARGET_BASES_1000X']
            elif 12500 < mean_cover:
                uniformity = row['PCT_TARGET_BASES_2500X']
            uniformity_pct = round((uniformity * 100), 2)

            metrics.append(f'{name}\t{reads}\t{bases}\t{mean_cover}\t{uniformity_pct}')

        with open(f"{outfile_name}.txt", "w") as outfile:
            outfile.write("\n".join(metrics))


def draw_cover_plot(folder: str):
    for file in os.listdir(folder):
        name = file.name.split('.')[0]
        df = pd.read_csv(file, sep='\t')
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
    new_line = []

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
            name = file.name.split('.')[0]
            df = pd.read_excel(file, sheet_name="Variant", header=1)
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

            dest_file = open(path_to_filtered / f'{name}.xlsx', 'wb')
            sorted_df.to_excel(dest_file, index=False, engine='openpyxl')
            print("Атлищна! отфильтрованные данные в ", dest_file)