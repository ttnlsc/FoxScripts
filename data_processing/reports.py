import os

clinical_significance = {
    'Pathogenic': ' как патогенный клинически значимый вариант',
    'Pathogenic/Likely pathogenic': ' как патогенный клинически значимый вариант',
    'Likely Pathogenic': ' как вероятно патогенный клинически значимый вариант',
    'Likely pathogenic': ' как вероятно патогенный клинически значимый вариант',
    'Uncertain significance': ' как вариант с неизвестным клиническим значением',
    'Conflicting interpretations of pathogenicity. risk factor': ', имеет конфликтующие интерпретации патогенности',
    'Conflicting interpretations of pathogenicity|risk factor': ', имеет конфликтующие интерпретации патогенности',
    'Conflicting interpretations of pathogenicity': ', имеет конфликтующие интерпретации патогенности',
    'not provided': ', критерии не предоставлены'
}

acmg_significance = {
    '5': 'патогенным вариантом',
    '4': 'вероятно патогенным вариантом',
    '3': 'вариантом с неизвестным клиническим значением'
}

sequence_ontology = {
    'missense_variant': 'замене аминокислоты',
    'frameshift_truncation': 'сдвигу рамки считывания и преждевременной остановке синтеза белка',
    'frameshift_elongation': 'сдвигу рамки считывания',
    'stop_gained': 'преждевременной остановке синтеза белка',
    'inframe_insertion': 'вставке без сдвига рамки считывания',
    'intron_variant': 'замене нуклеотида в интроне', 'splice_site_variant': 'изменению сайта сплайсинга',
    'synonymous_variant': 'изменению сайта сплайсинга',
    'start_lost': 'потере старт-кодона',
    'inframe_deletion': 'делеции без сдвига рамки считывания',
    '5_prime_UTR_variant': 'к нарушению сплайсинга'
}

zigosity = {'het': 'гетерозиготном', 'hom': 'гомозиготном'}


def read_data(filename):
    params = {}
    with open(filename, 'r') as file:
        content = file.read()

        splited_data = content.splitlines()

        header = splited_data[0].split(';')
        params['result'] = header.index('result')
        params['number'] = header.index('number')
        params['material'] = header.index('material')
        params['reads'] = header.index('PF_READS')
        # bases = header.index('PF_BASES')
        params['cover'] = header.index('MEAN_COVER')
        params['uniformity'] = header.index('Uniformity')
        params['numffpe'] = header.index('num_ffpe')
        params['gene'] = header.index('Gene')
        params['transcript'] = header.index('Transcript')
        params['chrom'] = header.index('Chrom')
        params['pos'] = header.index('Position')
        params['ref'] = header.index('Ref Base')
        params['alt'] = header.index('Alt Base')
        params['exon'] = header.index('Exon Number')
        params['coding'] = header.index('cDNA change')
        params['protein'] = header.index('Protein Change')
        params['effect'] = header.index('Sequence Ontology')
        params['clinvar'] = header.index('Clinical Significance')
        params['acmg'] = header.index('ACMG')
        params['amp'] = header.index('AMP')
        params['depth'] = header.index('Total reads')
        params['freq'] = header.index('Global AF')
        params['ab'] = header.index('Variant AF')
        params['rs'] = header.index('rsID')
        params['zyg'] = header.index('Zygosity')
        params['cosmic'] = header.index('ID')
        # params['panel'] = header.index('panel')

        for line in splited_data[1:]:
            line = line.split(';')
            if line[params['material']] == 'кровь':
                write_abc_blood_conclusion(params, line)
            elif line[params['material']] == 'блок':
                write_abc_ffpe_conclusion(params, line)
            else:
                raise ValueError(f"Wrong material {line[params['material']]}")


def write_abc_blood_conclusion(params, line):
    name = line[params['number']].split('/')[0]
    output_file = name + '.txt'

    conclusion = [
        'Выполнен анализ наличия мутаций в кодирующих и прилегающих областях генов BRCA1(ENST00000357654.9) и BRCA2(ENST00000380152.8).',
        'Исследование проводилось с использованием набора реагентов "Соло-тест АВС плюс" (РУ № РЗН 2023/20034) методом высокопроизводительного секвенирования по ТУ 21.20.23-005-91709359-2022.',
        'Гены, входящие в панель: BRCA1, BRCA2, ATM, CHEK2, PALB2.',
        'Диагностические характеристики запуска:',
        f"Количество прочтений на образец - {line[params['reads']]}",
        f"Среднее значение покрытия образца - {line[params['cover']]}",
        f"Равномерность покрытия - {line[params['uniformity']]}%"
    ]

    if line[params['result']] == 'wt':
        conclusion.append(
            "При исследовании ДНК, выделенной из лейкоцитов периферической крови, патогенных клинически значимых мутаций в генах BRCA1/BRCA2 не выявлено.")
    elif line[params['result']] == 'mut':
        if line[params['gene']] in ["BRCA1", "BRCA2"]:
            mut_description = ["При исследовании ДНК, выделенной из лейкоцитов периферической крови,"]
            if line[params['exon']] == '':
                mut_description.append(
                    f"в гене {line[params['gene']]}({line[params['transcript']]}) выявлена герминальная мутация {line[params['coding']]}")
                mut_description.append(
                    f"(hg38 {line[params['chrom']]}:{line[params['pos']]}-{line[params['ref']]}-{line[params['alt']]}, {line[params['protein']]}, {line[params['rs']]})")
                mut_description.append(
                    f"в {zigosity[line[params['zyg']]]} состоянии, приводящая к {sequence_ontology[line[params['effect']]]}.")
            else:
                mut_description.append(
                    f"в {line[params['exon']]} экзоне гена {line[params['gene']]}({line[params['transcript']]}) выявлена герминальная мутация {line[params['coding']]}")
                mut_description.append(
                    f"(hg38 {line[params['chrom']]}:{line[params['pos']]}-{line[params['ref']]}-{line[params['alt']]}, {line[params['protein']]}, {line[params['rs']]})")
                mut_description.append(
                    f"в {zigosity[line[params['zyg']]]} состоянии, приводящая к {sequence_ontology[line[params['effect']]]}.")
            mut_description = ' '.join(mut_description)
            conclusion.append(mut_description)
        else:
            conclusion.append(
                "При исследовании ДНК, выделенной из лейкоцитов периферической крови, патогенных клинически значимых мутаций в генах BRCA1/BRCA2 не выявлено.")
            mut_description = ["Примечание: "]
            if line[params['exon']] == '':
                mut_description.append(
                    f"в гене {line[params['gene']]}({line[params['transcript']]}) выявлена герминальная мутация {line[params['coding']]}")
                mut_description.append(
                    f"(hg38 {line[params['chrom']]}:{line[params['pos']]}-{line[params['ref']]}-{line[params['alt']]}, {line[params['protein']]}, {line[params['rs']]})")
                mut_description.append(
                    f"в {zigosity[line[params['zyg']]]} состоянии, приводящая к {sequence_ontology[line[params['effect']]]}.")
            else:
                mut_description.append(
                    f"в {line[params['exon']]} экзоне гена {line[params['gene']]}({line[params['transcript']]}) выявлена герминальная мутация {line[params['coding']]}")
                mut_description.append(
                    f"(hg38 {line[params['chrom']]}:{line[params['pos']]}-{line[params['ref']]}-{line[params['alt']]}, {line[params['protein']]}, {line[params['rs']]})")
                mut_description.append(
                    f"в {zigosity[line[params['zyg']]]} состоянии, приводящая к {sequence_ontology[line[params['effect']]]}.")
            mut_description = ' '.join(mut_description)
            conclusion.append(mut_description)
        if line[params['freq']] == '':
            conclusion.append("Вариант отсутствует в базе данных популяционных частот gnomAD_genomes v3.")
        else:
            conclusion.append(
                f"Частота минорного аллеля (gnomAD_genomes v3) {(float(line[params['freq']]) * 100): .4g}%.")
        conclusion.append(f"Глубина прочтения варианта {line[params['depth']]}x.")
        if line[params['clinvar']] != '' and line[params['rs']] != '':
            conclusion.append(
                f"Мутация {line[params['coding']]} в гене {line[params['gene']]} зарегистрирована в международных базах данных dbSNP, ClinVar{clinical_significance[line[params['clinvar']]]}.")
        elif line[params['clinvar']] != '' and line[params['rs']] == '':
            conclusion.append(
                f"Мутация {line[params['coding']]} в гене {line[params['gene']]} зарегистрирована в международной базе данных ClinVar{clinical_significance[line[params['clinvar']]]}.")
        elif line[params['clinvar']] == '' and line[params['rs']] != '':
            conclusion.append(
                f"Мутация {line[params['coding']]} в гене {line[params['gene']]} зарегистрирована в международной базе данных dbSNP.")
        else:
            conclusion.append(
                f"Мутация {line[params['coding']]} в гене {line[params['gene']]} не зарегистрирована в международных базах данных dbSNP, ClinVar.")
        conclusion.append(
            f"Согласно классификации ACMG мутация {line[params['coding']]} в гене {line[params['gene']]} является {acmg_significance[line[params['acmg']]]}, класс мутаций {line[params['acmg']]}.")
        if line[params['gene']] == "CHEK2" and line[params['coding']] == 'c.470T>C':
            conclusion.append(
                "Согласно клиническим рекомендациям NCCN, изменение тактики ведения для пациентов носителей данного варианта не требуется.")
        else:
            conclusion.append("Рекомендована консультация врача-генетика.")
    else:
        raise ValueError(f"Wrong result {line[params['result']]}")

    with open(output_file, mode='w') as outfile:
        outfile.write('\n'.join(conclusion))


def write_abc_ffpe_conclusion(params, line):
    name = line[params['number']].split('/')[0]
    output_file = name + '.txt'

    conclusion = [
        'Выполнен анализ наличия мутаций в кодирующих и прилегающих областях генов BRCA1(ENST00000357654.9) и BRCA2(ENST00000380152.8).',
        'Исследование проводилось с использованием набора реагентов "Соло-тест АВС плюс" (РУ № РЗН 2023/20034) методом высокопроизводительного секвенирования по ТУ 21.20.23-005-91709359-2022.',
        'Гены, входящие в панель: BRCA1, BRCA2, ATM, CHEK2, PALB2, PIK3CA.',
        f"Исследованный материал: блок № {line[params['numffpe']]}",
        'Диагностические характеристики запуска:',
        f"Количество прочтений на образец - {line[params['reads']]}",
        f"Среднее значение покрытия образца - {line[params['cover']]}",
        f"Равномерность покрытия - {line[params['uniformity']]}%"
    ]
    if line[params['result']] == 'wt':
        conclusion.append(
            "При исследовании ДНК, выделенной из фиксированной формалином и залитой в парафин ткани, патогенных клинически значимых мутаций в генах BRCA1/BRCA2 не выявлено.")
    elif line[params['result']] == 'mut':
        if line[params['gene']] in ["BRCA1", "BRCA2"]:
            mut_description = [
                "При исследовании ДНК, выделенной из фиксированной формалином и залитой в парафин ткани,"]
            if line[params['exon']] == '':
                mut_description.append(
                    f"в гене {line[params['gene']]}({line[params['transcript']]}) выявлена мутация {line[params['coding']]}")
                mut_description.append(
                    f"(hg38 {line[params['chrom']]}:{line[params['pos']]}-{line[params['ref']]}-{line[params['alt']]}, {line[params['protein']]}, {line[params['rs']]}),")
                mut_description.append(
                    f"аллельный баланс {(float(line[params['ab']]) * 100): .2g}%, приводящая к {sequence_ontology[line[params['effect']]]}.")
            else:
                mut_description.append(
                    f"в {line[params['exon']]} экзоне гена {line[params['gene']]}({line[params['transcript']]}) выявлена мутация {line[params['coding']]}")
                mut_description.append(
                    f"(hg38 {line[params['chrom']]}:{line[params['pos']]}-{line[params['ref']]}-{line[params['alt']]}, {line[params['protein']]}, {line[params['rs']]}),")
                mut_description.append(
                    f"аллельный баланс {(float(line[params['ab']]) * 100): .2g}%, приводящая к {sequence_ontology[line[params['effect']]]}.")
            mut_description = ' '.join(mut_description)
            conclusion.append(mut_description)
        else:
            conclusion.append(
                "При исследовании ДНК, выделенной из фиксированной формалином и залитой в парафин ткани, патогенных клинически значимых мутаций в генах BRCA1/BRCA2 не выявлено.")
            mut_description = ["Примечание:"]
            if line[params['exon']] == '':
                mut_description.append(
                    f"в гене {line[params['gene']]}({line[params['transcript']]}) выявлена мутация {line[params['coding']]}")
                mut_description.append(
                    f"(hg38 {line[params['chrom']]}:{line[params['pos']]}-{line[params['ref']]}-{line[params['alt']]}, {line[params['protein']]}, {line[params['rs']]}),")
                mut_description.append(
                    f"аллельный баланс {(float(line[params['ab']]) * 100): .2g}%, приводящая к {sequence_ontology[line[params['effect']]]}.")
            else:
                mut_description.append(
                    f"в {line[params['exon']]} экзоне гена {line[params['gene']]}({line[params['transcript']]}) выявлена мутация {line[params['coding']]}")
                mut_description.append(
                    f"(hg38 {line[params['chrom']]}:{line[params['pos']]}-{line[params['ref']]}-{line[params['alt']]}, {line[params['protein']]}, {line[params['rs']]}),")
                mut_description.append(
                    f"аллельный баланс {(float(line[params['ab']]) * 100): .2g}%, приводящая к {sequence_ontology[line[params['effect']]]}.")
            mut_description = ' '.join(mut_description)
            conclusion.append(mut_description)
        if line[params['freq']] == '':
            conclusion.append("Вариант отсутствует в базе данных популяционных частот gnomAD_genomes v3.")
        else:
            conclusion.append(
                f"Частота минорного аллеля (gnomAD_genomes v3) {(float(line[params['freq']]) * 100): .4g}%.")
        conclusion.append(f"Глубина прочтения варианта {line[params['depth']]}x.")
        if line[params['gene']] not in ["PIK3CA", "BRAF"]:
            if line[params['clinvar']] != '' and line[params['rs']] != '':
                conclusion.append(
                    f"Мутация {line[params['coding']]} в гене {line[params['gene']]} зарегистрирована в международных базах данных dbSNP, ClinVar{clinical_significance[line[params['clinvar']]]}.")
            elif line[params['clinvar']] != '' and line[params['rs']] == '':
                conclusion.append(
                    f"Мутация {line[params['coding']]} в гене {line[params['gene']]} зарегистрирована в международной базе данных ClinVar{clinical_significance[line[params['clinvar']]]}.")
            elif line[params['clinvar']] == '' and line[params['rs']] != '':
                conclusion.append(
                    f"Мутация {line[params['coding']]} в гене {line[params['gene']]} зарегистрирована в международной базе данных dbSNP.")
            else:
                conclusion.append(
                    f"Мутация {line[params['coding']]} в гене {line[params['gene']]} не зарегистрирована в международных базах данных dbSNP, ClinVar.")
            conclusion.append(
                f"Согласно классификации ACMG мутация {line[params['coding']]} в гене {line[params['gene']]} является {acmg_significance[line[params['acmg']]]}, класс мутаций {line[params['acmg']]}.")
            if line[params['gene']] == "CHEK2" and line[params['coding']] == 'c.470T>C':
                conclusion.append(
                    "Согласно клиническим рекомендациям NCCN, изменение тактики ведения для пациентов носителей данного варианта не требуется.")
            else:
                conclusion.append(f"Класс по AMP - {line[params['amp']]}")
                conclusion.append(
                    "С целью установления характера мутации: герминальный или соматический, рекомендовано выполнить исследование на наличие данной мутации в образце крови.")
                conclusion.append("Рекомендована консультация врача-генетика.")
        else:
            if line[params['cosmic']] != '':
                conclusion.append(
                    f"Мутация {line[params['coding']]} в гене {line[params['gene']]} зарегистрирована в международной базе данных COSMIC ({line[params['cosmic']]}).")
            conclusion.append(f"Класс по AMP - {line[params['amp']]}.")
    else:
        raise ValueError(f"Wrong result {line[params['result']]}")

    with open(output_file, mode='w') as outfile:
        outfile.write('\n'.join(conclusion))
