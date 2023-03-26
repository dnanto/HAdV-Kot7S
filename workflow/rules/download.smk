from pathlib import Path
from csv import DictReader
from datetime import datetime

from snakemake.utils import validate


configfile: "config/download.config.yaml"


validate(config, schema=Path(workflow.basedir) / "schemas" / "download.schema.yaml")

species = config["species"]
crossref = config["crossref"][species]

results = Path("results")
root_data = results / "data" / species
root_download = results / "download" / species


rule download:
    message:
        "download"
    output:
        ref=root_download / f"{species}.ref.zip",
        lib=root_download / f"{species}.lib.zip",
    threads: 1
    params:
        assembly=crossref["assembly"],
        taxon=crossref["taxon"],
    threads: 1
    shell:
        """
        datasets download genome accession {params.assembly:q} \
            --include genome,gff3,gbff \
            --filename {output.ref:q}
        datasets download virus genome taxon {params.taxon:q} \
            --complete-only \
            --include genome \
            --filename {output.lib:q}
        """


rule ref:
    message:
        "ref"
    input:
        zip=rules.download.output.ref,
    output:
        fna=root_data / "ref.fna",
        gbf=root_data / "ref.gbff",
        gff=root_data / "ref.gff",
    log:
        fna=root_data / "ref.fna.log",
        gbf=root_data / "ref.gbff.log",
        gff=root_data / "ref.gff.log",
    params:
        list="ncbi_dataset/data/genomic.fna",
    threads: 1
    shell:
        """
        unzip -p {input.zip:q} "ncbi_dataset/data/*/*_genomic.fna" 1> {output.fna:q} 2> {log.fna:q}
        unzip -p {input.zip:q} "ncbi_dataset/data/*/genomic.gbff"  1> {output.gbf:q} 2> {log.gbf:q}
        unzip -p {input.zip:q} "ncbi_dataset/data/*/genomic.gff"   1> {output.gff:q} 2> {log.gff:q}
        """


rule lib:
    message:
        "lib"
    input:
        zip=rules.download.output.lib,
    output:
        fna=root_data / "lib.fna",
    log:
        fna=root_data / "lib.fna.log",
    threads: 1
    shell:
        """
        unzip -p {input.zip:q} "ncbi_dataset/data/genomic.fna" 1> {output.fna:q} 2> {log.fna:q}
        """


rule meta:
    message:
        "meta"
    input:
        zip=rules.download.output.lib,
    output:
        tsv=root_data / "meta.tsv",
    log:
        tsv=root_data / "meta.tsv.log",
    params:
        list="ncbi_dataset/data/data_report.jsonl",
        script=Path(workflow.basedir) / "scripts" / "meta.jq",
        species=species,
    threads: 1
    shell:
        """
        unzip -p {input.zip:q} {params.list:q} | \
        jq --arg species {params.species} -f {params.script:q} | \
        jq -r -s '(.[0] | keys_unsorted), (.[] | [.[]]) | @tsv' \
            1> {output.tsv:q} \
            2> {log.tsv:q}
        """


rule filter:
    message:
        "filter"
    input:
        fna=rules.lib.output.fna,
        tsv=rules.meta.output.tsv,
    output:
        fna=root_data / "sbj.fna",
        tsv=root_data / "sbj.tsv",
    log:
        log=root_data / "sbj.log",
    params:
        id="accession",
        resolution="month",
    shell:
        """
        augur filter \
            --metadata {input.tsv:q} \
            --metadata-id-columns {params.id:q} \
            --sequences {input.fna:q} \
            --exclude-ambiguous-dates-by {params.resolution:q} \
            --output {output.fna:q} \
            --output-metadata {output.tsv} \
            2> {log.log:q}
        """


rule align:
    message:
        "align"
    input:
        fna=rules.filter.output.fna,
    output:
        fna=root_data / "msa.fna",
    log:
        fna=root_data / "msa.fna.log",
    params:
        reference=crossref["refseq"],
    threads: 8
    shell:
        """
        mafft --adjustdirection --thread {threads} {input.fna:q} 1> {output.fna:q} 2> {log.fna:q}
        """


rule recom:
    input:
        fna=rules.align.output.fna,
    output:
        tree=root_data / "gub.final_tree.tre",
        node=root_data / "gub.node_labelled.final_tree.tre",
        stat=root_data / "gub.per_branch_statistics.csv",
    params:
        prefix=root_data / "gub",
        bootstrap=10,
        iterations=1,
    threads: 8
    shadow:
        "shallow"
    shell:
        """
        run_gubbins.py {input.fna:q} \
            --prefix {params.prefix:q} \
            --threads {threads} \
            --model-fitter raxml \
            --bootstrap {params.bootstrap:q} \
            --transfer-bootstrap \
            --sh-test \
            --best-model \
            --iterations {params.iterations:q} \
            --extensive-search && \
        rm -rf "$(pwd)"/tmp*
        """


rule all:
    input:
        rules.download.output.lib,
        rules.download.output.ref,
        rules.ref.output.fna,
        rules.ref.output.gbf,
        rules.ref.output.gff,
        rules.lib.output.fna,
        rules.meta.output.tsv,
        rules.filter.output.fna,
        rules.filter.output.tsv,
        rules.align.output.fna,
        rules.recom.output.tree,
        rules.recom.output.node,
        rules.recom.output.stat,
