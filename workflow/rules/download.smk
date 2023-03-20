from datetime import datetime
from pathlib import Path
from csv import DictReader

from snakemake.utils import validate


configfile: "config/download.config.yaml"


validate(config, schema="../schemas/download.schema.yaml")


results = Path("results")
log = Path("logs") / datetime.today().isoformat().replace(":", "")

with open(config["path"]) as file:
    species = [row["species"] for row in DictReader(file, delimiter="\t")]


def get_species(wildcards):
    output = checkpoints.download.get(**wildcards).output[0]
    return [ele.stem for ele in sorted(Path(output).glob("*.zip"))]


checkpoint download:
    message:
        "download raw data"
    input:
        config["path"],
    output:
        directory(results / "download"),
    threads: 1
    log:
        out=log / "download.out.txt",
        err=log / "download.err.txt",
    shell:
        """
        mkdir -p {output:q}
        while IFS=$'\t' read -r species taxon _; do
            echo "$species -> $taxon"
            datasets download virus genome taxon "$taxon" \
                --annotated \
                --complete-only \
                --include genome,annotation \
                --filename {output:q}/"$species.zip"
        done < <(tail -n +2 {input:q}) 1> {log.out:q} 2> {log.err:q}
        """


rule unzip:
    message:
        "unzip files: {wildcards.species}"
    input:
        results / "download" / "{species}.zip",
    output:
        results / "data" / "{species}" / "ncbi_dataset" / "data" / "data_report.jsonl",
    params:
        root=lambda wildcards: results / "data" / wildcards.species,
    threads: 1
    log:
        out=log / "unzip.{species}.out.txt",
        err=log / "unzip.{species}.err.txt",
    shell:
        """
        mkdir -p {params.root:q} && unzip -d {params.root:q} {input:q} 1> {log.out:q} 2> {log.err:q}
        """


rule meta:
    message:
        "metadata"
    input:
        lambda wildcards: expand(rules.unzip.output[0], species=get_species(wildcards)),
    output:
        results / "data" / "meta.tsv",
    log:
        err=log / "meta.err.txt",
    params:
        script=Path("workflow") / "scripts" / "meta.jq",
    threads: 1
    shell:
        """
        jq -f {params.script:q} {input:q} | \
            jq -r -s '(.[0] | keys_unsorted), (.[] | [.[]]) | @tsv' \
            1> {output:q} \
            2> {log.err:q}
        """
