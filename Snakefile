# {{{1 Preamble

# import pandas as pd
# import math
# from lib.snake import curl_recipe, alias_recipe

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
# wildcard_constraints:  # TODO

# {{{1 Configuration

MAX_THREADS = 99
if 'MAX_THREADS' in config:
    MAX_THREADS = config['MAX_THREADS']

# {{{2 Include sub-pipelines

include: 'snake/local.snake'

# {{{2 Configure default actions

rule all:
    output: []

# {{{1 Utility rules

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

rule download_ibdmdb_metadata:
    output: 'raw/metadata/hmp2_metadata_2018-08-20.csv'
    params:
        url='ftp://ftp.broadinstitute.org/metadata/hmp2_metadata_2018-08-20.csv',
        user='public',
        password='hmp2_ftp',
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
