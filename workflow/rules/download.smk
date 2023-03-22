from pathlib import Path
from csv import DictReader
from datetime import datetime

from snakemake.utils import validate


configfile: "config/download.config.yaml"


validate(config, schema="../schemas/download.schema.yaml")

species = config["species"]
crossref = config["crossref"][species]

results = Path("results")
log = Path("logs") / datetime.today().isoformat().replace(":", "")


rule download:
    message:
        "download"
    output:
        ref=results / "download" / f"{species}.ref.zip",
        lib=results / "download" / f"{species}.lib.zip",
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
        fna=results / "data" / species / "ref.fna",
        gbf=results / "data" / species / "ref.gbff",
        gff=results / "data" / species / "ref.gff",
    log:
        err=log / f"{species}.ref.err.txt",
    params:
        list="ncbi_dataset/data/genomic.fna",
    threads: 1
    shell:
        """
        unzip -p {input.zip:q} "ncbi_dataset/data/*/*_genomic.fna" 1> {output.fna:q} 2>  {log.err:q}
        unzip -p {input.zip:q} "ncbi_dataset/data/*/genomic.gbff"  1> {output.gbf:q} 2>> {log.err:q}
        unzip -p {input.zip:q} "ncbi_dataset/data/*/genomic.gff"   1> {output.gff:q} 2>> {log.err:q}
        """


rule lib:
    message:
        "lib"
    input:
        zip=rules.download.output.lib,
    output:
        fna=results / "data" / species / "lib.fna",
    log:
        err=log / f"{species}.lib.err.txt",
    threads: 1
    shell:
        """
        unzip -p {input.zip:q} "ncbi_dataset/data/genomic.fna" 1> {output.fna:q} 2> {log.err:q}
        """


rule meta:
    message:
        "meta"
    input:
        zip=rules.download.output.lib,
    output:
        tsv=results / "data" / species / "meta.tsv",
    log:
        err=log / f"{species}.meta.err.txt",
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
            2> {log.err:q}
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
