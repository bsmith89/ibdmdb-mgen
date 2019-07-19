# {{{1 Preamble

import pandas as pd
import os
from lib.snake import curl_recipe, alias_recipe, hardlink_recipe
import sqlite3
import sys

limit_numpy_procs = \
        """
        export MKL_NUM_THREADS={threads}
        export OPENBLAS_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        export VECLIB_MAXIMUM_THREADS={threads}

        """

# Utility wildcard constrains
noperiod='[^.]+'
wildcard_constraints:
    r='r1|r2|r3',
    group=noperiod,
    library=noperiod,

# {{{1 Configuration

MAX_THREADS = 99
if 'MAX_THREADS' in config:
    MAX_THREADS = config['MAX_THREADS']

# {{{2 Include sub-pipelines

include: 'snake/local.snake'

# {{{2 Metadata config

configfile: 'config.yaml'

if os.path.exists(config['_meta_db']):  # USUALLY data/meta.db
    meta_db_conn = sqlite3.connect(config['_meta_db'])

    config['library'] = pd.read_sql('SELECT * FROM library;', con=meta_db_conn, index_col='library_id')

    config['library_group'] = {}
    _library_group = pd.read_sql('SELECT * FROM library_group;', con=meta_db_conn)
    for library_group, d in _library_group.groupby('library_group'):
        config['library_group'][library_group] = {}
        config['library_group'][library_group]['all_libraries'] = d.library_id.to_list()

    d0 = pd.read_sql(('SELECT * FROM mgen_library '
                      'JOIN preparation USING (preparation_id) '
                      'JOIN stool USING (stool_id) '
                      'JOIN visit USING (visit_id) '
                      'JOIN library_group USING (library_id);'),
                     con=meta_db_conn)
    for library_group, d1 in d0.groupby('library_group'):
        # config['library_group'][library_group] = {}  # Already done above.
        config['library_group'][library_group]['subject'] = {}
        for subject_id, d2 in d1.groupby('subject_id'):
            config['library_group'][library_group]['subject'][subject_id] = d2.library_id.to_list()

    config['all_subjects'] = d0.subject_id.unique()

def library_to_external_id(library_id):
    return config['library'].loc[library_id].external_id

# {{{2 Configure default actions

rule all:
    output: []

# {{{1 Utility rules

rule drop_table_header:
    output: '{stem}.noheader.{suffix}'
    input: '{stem}.{suffix}'
    wildcard_constraints:
        suffix=noperiod
    shell: "sed '1,1d' {input} > {output}"

rule rename_fasta_with_long_suffix:
    output: '{stem}.fasta'
    input: '{stem}.fn'
    shell: hardlink_recipe

rule build_meta_db:
    output: 'data/meta.db'
    input:
        schema='schema.meta.sql',
        subject='meta/subject.noheader.tsv',
        visit='meta/visit.noheader.tsv',
        stool='meta/stool.noheader.tsv',
        preparation='meta/preparation.noheader.tsv',
        library='meta/library.noheader.tsv',
        library_group='meta/library_group.noheader.tsv',
    shell:
        r"""
        echo '
.bail on
PRAGMA foreign_keys = TRUE;
.read {input.schema}
.separator \t
.import {input.subject} subject
.import {input.visit} visit
.import {input.stool} stool
.import {input.preparation} preparation
.import {input.library} library
.import {input.library_group} library_group
             ' \
            | sqlite3 {output}
        """


rule initialize_project:
    shell:
        '''
        git config --local filter.dropoutput_ipynb.clean scripts/ipynb_output_filter.py
        git config --local filter.dropoutput_ipynb.smudge cat
        git config --local diff.daff-csv.command "daff.py diff --git"
        git config --local merge.daff-csv.name "daff.py tabular merge"
        git config --local merge.daff-csv.driver "daff.py merge --output %A %O %A %B"
        echo 'Please activate your environment and then run `pip install -r requirements.txt` or analogous.'
        '''

rule start_jupyter:
    shell: 'jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/'

rule start_ipython:
    threads: MAX_THREADS
    shell: limit_numpy_procs + 'ipython'

rule query_db:
    output: 'data/{db}.select_{query}.tsv'
    input:
        db='data/{db}.db',
        query='scripts/query/{query}.sql',
    shell:
        """
        sqlite3 -header -separator '\t' {input.db} < {input.query} > {output}
        """

# {{{1 Generalizable rules

rule bowtie_index_build:
    output:
        '{stem}.1.bt2',
        '{stem}.2.bt2',
        '{stem}.3.bt2',
        '{stem}.4.bt2',
        '{stem}.rev.1.bt2',
        '{stem}.rev.2.bt2'
    input: '{stem}.fn'
    shell:
        """
        bowtie2-build {input} {wildcards.stem}
        """

rule bowtie_index_build_from_gzipped:
    output:
        '{stem}.1.bt2',
        '{stem}.2.bt2',
        '{stem}.3.bt2',
        '{stem}.4.bt2',
        '{stem}.rev.1.bt2',
        '{stem}.rev.2.bt2'
    input: '{stem}.fn.gz'
    params: db=lambda w: f'{w.stem}'
    shell:
        """
        tmp=$(mktemp)
        echo $tmp
        zcat {input} > $tmp
        bowtie2-build $tmp {params.db}
        """

rule index_bam:
    output: 'data/{stem}.sort.bam.bai'
    input: 'data/{stem}.sort.bam'
    shell: 'samtools index {input} {output}'

rule make_diamond_db:
    output: 'ref/{db}.dmnd'
    input: 'ref/{db}.fa'
    shell:
        """
        diamond makedb -d ref/{db} < {input}
        """

rule diamond_search_fasta:
    output: 'data/{query}.{db}-blastp.tsv'
    input: fasta='data/{query}.fn', db='ref/{db}.dmnd'
    params:
        db='ref/{db}'
    threads: MAX_THREADS
    shell:
        """
        diamond blastx --threads {threads} --db {params.db} --query {input.fasta} > {output}
        """

rule diamond_search_fastq:
    output: 'data/{query}.{db}-blastp.tsv'
    input: fasta='data/{query}.fq.gz', db='ref/{db}.dmnd'
    params:
        db='ref/{db}'
    threads: 4
    shell:
        """
        seqtk seq -A {input.fasta} | diamond blastx --threads {threads} --db {params.db} --query - > {output}
        """

ruleorder: bowtie_index_build > bowtie_index_build_from_gzipped

rule count_seq_lengths_nucl:
    output: '{stem}.nlength.tsv'
    input: script='scripts/count_seq_lengths.py', seqs='{stem}.fn'
    shell:
        r"""
        {input.script} {input.seqs} > {output}
        """

# {{{1 Downloads

# {{{2 Download and organize reference data

rule download_illumina_adapters:
    output: 'raw/ref/illumina_adapters.fn'
    params:
        url='https://raw.githubusercontent.com/vsbuffalo/scythe/master/illumina_adapters.fa'
    resources:
        network_connections=1,
    shell: curl_recipe

rule link_illumina_adapters:
    output: 'ref/illumina_adapters.fn'
    input: 'raw/ref/illumina_adapters.fn'
    shell: alias_recipe

rule download_GRCh38_index:
    output: 'raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz'
    params:
        url='ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz'
    resources:
        network_connections=1,
    shell: curl_recipe

rule unpack_GRCh38_index:
    output:
        'raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz'
        'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2',
        'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.2.bt2',
        'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.3.bt2',
        'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.4.bt2',
        'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.1.bt2',
        'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.2.bt2'
    input: 'raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz'

rule alias_GRCh38_index_file:
    output: 'ref/GRCh38.{suffix}'
    input: 'raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.{suffix}'
    shell: alias_recipe

ruleorder: alias_GRCh38_index_file > bowtie_index_build

rule download_checkm_data:
    output: 'raw/ref/checkm_data_2015_01_16.tar.gz'
    params:
        url='https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz'
    resources:
        network_connections=1,
    shell:
        curl_recipe

rule unpack_checkm_data:
    output: 'raw/ref/checkm_db'
    input: 'raw/ref/checkm_data_2015_01_16.tar.gz'
    shell:
        """
        rm -rf {output}
        mkdir -p {output}
        tar -C {output} -xzf {input}
        """

rule alias_checkm_data:
    output: directory('ref/checkm_db')
    input: 'raw/ref/checkm_db'
    shell: alias_recipe

# {{{2 Download IBDMDB Raw Data

rule download_ibdmdb_metadata:
    output: 'raw/metadata/hmp2_metadata_2018-08-20.csv'
    params:
        url='ftp://ftp.broadinstitute.org/metadata/hmp2_metadata_2018-08-20.csv',
        user='public',
        password='hmp2_ftp',
    resources:
        network_connections=1
    shell:
        'curl --user {params.user}:{params.password} {params.url} > {output}'

rule extract_metadata_tables:
    output:
        subject='meta/subject.tsv',
        visit='meta/visit.tsv',
        stool='meta/stool.tsv',
        preparation='meta/preparation.tsv',
        library='meta/library.tsv',
    input:
        script='scripts/process_metadata_tables.py',
        raw='raw/metadata/hmp2_metadata_2018-08-20.csv',
    shell:
        """
        cat {input.raw} | {input.script} {output.subject} {output.visit} {output.stool} {output.preparation} {output.library}
        """

rule download_mgen_from_ibdmdb:
    output: temp('raw/2018-05-04-mgen/{stem}.tar')
    params:
        url='ftp.broadinstitute.org/raw/HMP2/MGX/2018-05-04',
        username='public',
        password='hmp2_ftp',
    resources:
        network_connections=1,
    log: 'log/{stem}.wget.log'
    shell:
        """
        wget --no-verbose -o {log} -P raw/2018-05-04-mgen ftp://{params.username}:{params.password}@{params.url}/{wildcards.stem}.tar
        """

rule unpack_mgen_from_ibdmdb:
    output:
        r1='raw/2018-05-04-mgen/{stem}_R1.fastq.gz',
        r2='raw/2018-05-04-mgen/{stem}_R2.fastq.gz',
    input:
        tar='raw/2018-05-04-mgen/{stem}.tar'
    shell:
        """
        tar -C raw/2018-05-04-mgen/ -xf {input.tar}
        """

# {{{1 Process Metagenomic Libraries

rule alias_mgen_library_r1:
    output: 'data/{library}.m.r1.fq.gz'
    input: lambda w: 'raw/2018-05-04-mgen/{}_R1.fastq.gz'.format(library_to_external_id(w.library))
    shell: alias_recipe

rule alias_mgen_library_r2:
    output: 'data/{library}.m.r2.fq.gz'
    input: lambda w: 'raw/2018-05-04-mgen/{}_R2.fastq.gz'.format(library_to_external_id(w.library))
    shell: alias_recipe

rule deduplicate_reads:
    output:
        r1=temp('data/{stem}.dedup.r1.fq.gz'),
        r2=temp('data/{stem}.dedup.r2.fq.gz')
    input:
        script='scripts/fastuniq_wrapper.sh',
        r1='data/{stem}.r1.fq.gz',
        r2='data/{stem}.r2.fq.gz'
    resources:
        mem_mb=24000
    shell: "{input.script} {input.r1} {input.r2} {output.r1} {output.r2}"

rule trim_adapters:
    output:
        fq=temp('data/{stem}.deadapt.{r}.fq.gz')
    input:
        adapters='ref/illumina_adapters.fn',
        fq='data/{stem}.{r}.fq.gz',
    log: 'log/{stem}.scythe.{r}.log'
    threads: 2
    shell:
        """
        scythe -a {input.adapters} {input.fq} >{output.fq} 2>{log}
        ! grep -Fxq 'Blank FASTA header or sequence in adapters file.' {log}
        """

rule quality_trim_reads:
    output:
        r1=temp('data/{stem}.qtrim.r1.fq.gz'),
        r2=temp('data/{stem}.qtrim.r2.fq.gz'),
        r3=temp('data/{stem}.qtrim.r3.fq.gz'),
    input:
        r1='data/{stem}.r1.fq.gz',
        r2='data/{stem}.r2.fq.gz',
    params:
        qual_type='sanger',
        qual_thresh=20
    shell:
        r"""
        sickle pe -t {params.qual_type} -q {params.qual_thresh} --gzip-output \
            --pe-file1 {input.r1} --pe-file2 {input.r2} \
            --output-pe1 {output.r1} --output-pe2 {output.r2} \
            --output-single {output.r3}
        """

rule filter_out_host:
    output:
        r1=temp('data/{stem}.hfilt.r1.fq.gz'),
        r2=temp('data/{stem}.hfilt.r2.fq.gz'),
    input:
        script='scripts/filter_out_mapping_reads.sh',
        r1='data/{stem}.r1.fq.gz',
        r2='data/{stem}.r2.fq.gz',
        index=['ref/GRCh38.1.bt2',
               'ref/GRCh38.2.bt2',
               'ref/GRCh38.3.bt2',
               'ref/GRCh38.4.bt2',
               'ref/GRCh38.rev.1.bt2',
               'ref/GRCh38.rev.2.bt2']
    params:
        index='ref/GRCh38',
    shell:
        """
        {input.script} {params.index} {input.r1} {input.r2} {output.r1} {output.r2}
        """

rule rename_cleaned_reads:
    output: 'data/{stem}.proc.{r}.fq.gz'
    input: 'data/{stem}.dedup.deadapt.qtrim.hfilt.{r}.fq.gz'
    shell: hardlink_recipe

# {{{1 Assembly

# 'sa' for Subject Assembly (co-assemble all libraries from that subject).
rule co_assemble_mgen_by_subject:
    output:
        fasta='data/{group}.{subject}.sa.fn',
        fastg='data/{group}.{subject}.sa.fg',
        dir=directory('data/{group}.{subject}.sa.megahit.d')
    input:
        r1=lambda w: [f'data/{library}.m.proc.r1.fq.gz'
                      for library in config['library_group'][w.group]['subject'][w.subject]],
        r2=lambda w: [f'data/{library}.m.proc.r2.fq.gz'
                      for library in config['library_group'][w.group]['subject'][w.subject]],
    params:
        r1_list=lambda w, input: ','.join(input.r1),
        r2_list=lambda w, input: ','.join(input.r2),
        k_list='21,41,61,81,101',
        k_max=101,
    threads: MAX_THREADS
    log: 'data/{group}.{subject}.sa.megahit.log'
    shell:
        r'''
        megahit \
            -1 {params.r1_list} \
            -2 {params.r2_list} \
            --k-list {params.k_list} \
            --out-dir {output.dir} \
            --num-cpu-threads {threads} \
            --verbose \
            2>&1 | tee {log}
        sed 's:k{params.k_max}_\(\S\+\).*$:\1:' {output.dir}/final.contigs.fa > {output.fasta}
        megahit_toolkit contig2fastg {params.k_max} {output.dir}/intermediate_contigs/k{params.k_max}.contigs.fa > {output.fastg}
        '''

# 'aa' for All Assembly (co-assemble all libraries in the group).
rule co_assemble_mgen_by_group:
    output:
        fasta='data/{group}.aa.fn',
        fastg='data/{group}.aa.fg',
        dir=directory('data/{group}.aa.megahit.d')
    input:
        r1=lambda w: [f'data/{library}.m.proc.r1.fq.gz'
                      for library in config['library_group'][w.group]['all_libraries']],
        r2=lambda w: [f'data/{library}.m.proc.r2.fq.gz'
                      for library in config['library_group'][w.group]['all_libraries']],
    params:
        r1_list=lambda w, input: ','.join(input.r1),
        r2_list=lambda w, input: ','.join(input.r2),
        k_list='21,41,61,81,101',
        k_max=101,
    threads: MAX_THREADS
    log: 'data/{group}.aa.megahit.log'
    shell:
        r'''
        megahit \
            -1 {params.r1_list} \
            -2 {params.r2_list} \
            --k-list {params.k_list} \
            --out-dir {output.dir} \
            --num-cpu-threads {threads} \
            --verbose \
            2>&1 | tee {log}
        sed 's:k{params.k_max}_\(\S\+\).*$:\1:' {output.dir}/final.contigs.fa > {output.fasta}
        megahit_toolkit contig2fastg {params.k_max} {output.dir}/intermediate_contigs/k{params.k_max}.contigs.fa > {output.fastg}
        '''


# 'ma' for 'merged assembly'
rule merge_subject_assemblies:
    output: merged_asmbl=directory('data/{group}.ma.d')
    input:
        subject_asmbl=lambda w: [f'data/{{group}}.{subject}.sa.fasta' for subject in config['library_group'][w.group]['subject']]
    params:
        input_args=lambda w, input: ' '.join(f'--trusted-contigs {contigs}' for contigs in input.subject_asmbl)
    threads: MAX_THREADS
    shell:
        """
        dummy_reads=$(mktemp --suffix=.fasta)
        echo ">dummy_1" > $dummy
        echo "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" >> $dummy
        echo ">dummy_2" >> $dummy
        echo "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" >> $dummy

        spades.py {params.input_args} -o {output} \
                --only-assembler --disable-rr --tmp-dir $TMPDIR \
                -t {threads} --12 $dummy_reads
        """

# {{{1 Assembly-Free Analysis

rule run_iggsearch_species:
    output:
        dir=directory('data/{stem}.iggmidas/iggsearch'),
        # bam='data/{stem}.iggmidas/iggsearch/mapped_reads.bam',
        # species='data/{stem}.iggmidas/iggsearch/species_profile.tsv',
    input:
        r1='data/{stem}.r1.fq.gz',
        r2='data/{stem}.r2.fq.gz'
    params:
        dir=lambda w: f'data/{w.stem}.iggmidas/iggsearch',
        runscript='$HOME/Projects/IGGsearch/run_iggsearch.py',
        db='/pollard/data/IGGdb/v1.0.0/iggsearch/iggdb_v1.0.0',
    log:
        read_alignment='data/{stem}.iggmidas/iggsearch/read_alignment.log'
    threads: 4
    shell:
        '''
        {params.runscript} search --threads {threads} --db_dir {params.db} --m1 {input.r1} --m2 {input.r2} --outdir {params.dir}
        '''

rule iggsearch_all_from_subject:
    output: touch('data/{group}.{subject}.iggsearch.touch')
    input:
        dir=lambda w: [f'data/{library}.m.proc.iggmidas/iggsearch'
                      for library in config['library_group'][w.group]['subject'][w.subject]],
