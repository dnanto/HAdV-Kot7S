from datetime import datetime
from pathlib import Path

from snakemake.utils import validate


configfile: "config/download.config.yaml"
validate(config, schema="../schemas/download.schema.yaml")

results = Path("results")
log = Path("logs") / datetime.today().isoformat().replace(":", "")

def get_download_output(wildcards):
    output = checkpoints.download.get(**wildcards).output[0]
    return list(Path(output).glob("*.zip"))


checkpoint download:
    message: "download raw data"
    input: config["path"]
    output: directory(results / "download")
    threads: 1
    log:
        out=log / "download.out.txt",
        err=log / "download.err.txt"
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
        done < <(tail -n +2 {input:q}) > {log.out:q} 2> {log.err:q}
        """

rule unzip:
    message: "unzip files"
    input: get_download_output
    output: directory(results / "data")
    log:
        out=log / "unzip.out.txt",
        err=log / "unzip.err.txt"
    shell:
        """
        mkdir -p {output:q}
        for ele in {input:q}; do
            d={output:q}/"$(basename "${{ele/.zip/}}")"
            unzip -d "$d" "$ele" > {log.out:q} 2> {log.err:q}
        done
        """
