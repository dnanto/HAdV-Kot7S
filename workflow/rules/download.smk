from pathlib import Path
from csv import DictReader
from datetime import datetime
from shlex import quote

from snakemake.utils import validate


configfile: "config/download.config.yaml"


validate(config, schema="../schemas/download.schema.yaml")

species = list(config["species"])

results = Path("results")
log = Path("logs") / datetime.today().isoformat().replace(":", "")


rule download:
    message:
        "download raw data"
    output:
        ref=results / "download" / "{species}" / "ref.zip",
        lib=results / "download" / "{species}" / "lib.zip",
    threads: 1
    log:
        out=log / "{species}.download.out.txt",
        err=log / "{species}.download.err.txt",
    params:
        root=results / "download",
        json=quote(json.dumps(config["species"])),
    threads: 1
    shell:
        """
        echo {params.json:1} | \
        jq -r 'to_entries | map([.key, .value.assembly, .value.taxon] | @tsv)[]' | \
        while read -r species assembly taxon; do
            mkdir -p {params.root:q}/"$species"
            echo "$species -> $assembly"
            datasets download genome accession "$assembly" \
                --include genome,gff3,gbff \
                --filename {params.root:q}/"$species/ref.zip"
            echo "$species -> $taxon"
            datasets download virus genome taxon "$taxon" \
                --complete-only \
                --include genome \
                --filename {params.root:q}/"$species/lib.zip"
        done > >(tee -a {log.out:q}) 2> >(tee -a {log.err:q} >&2)
        """


rule ref:
    message:
        "library"
    input:
        zip=rules.download.output.ref,
    output:
        fna=results / "data" / "{species}" / "ref.fna",
        gbf=results / "data" / "{species}" / "ref.gbff",
        gff=results / "data" / "{species}" / "ref.gff",
    log:
        err=log / "{species}.ref.err.txt",
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
        "library"
    input:
        zip=rules.download.output.lib,
    output:
        fna=results / "data" / "{species}" / "lib.fna",
    log:
        err=log / "{species}.lib.err.txt",
    threads: 1
    shell:
        """
        unzip -p {input.zip:q} "ncbi_dataset/data/genomic.fna" 1> {output.fna:q} 2> {log.err:q}
        """


rule meta:
    message:
        "metadata"
    input:
        zip=rules.download.output.lib,
    output:
        tsv=results / "data" / "{species}" / "meta.tsv",
    log:
        err=log / "{species}.meta.err.txt",
    params:
        list="ncbi_dataset/data/data_report.jsonl",
        script=Path(workflow.basedir) / "scripts" / "meta.jq",
        species=lambda wildcards: wildcards.species,
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
        expand(rules.download.output.lib, species=species),
        expand(rules.download.output.ref, species=species),
        expand(rules.ref.output.fna, species=species),
        expand(rules.ref.output.gbf, species=species),
        expand(rules.ref.output.gff, species=species),
        expand(rules.lib.output.fna, species=species),
        expand(rules.meta.output.tsv, species=species),
