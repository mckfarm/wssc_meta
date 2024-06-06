rule get_cds:
    input:
         "results/megahit/{sample}/final.contigs.fa"
    output:
        protein = "results/pyrodigal/{sample}/cds_proteins.faa",
        nucl = "results/pyrodigal/{sample}/cds_nucl.fa",
        gff = "results/pyrodigal/{sample}/cds.gff"
    params:
        index_name = "results/coverage/{sample}/{sample}"
    threads: 10
    conda: "../envs/pyrodigal.yaml"
    resources:
        mem="30G",
        time="05:00:00"
    shell:
        """
        pyrodigal -i {input} -a {output.protein} -d {output.nucl} -f gff -o {output.gff} -p meta --jobs {threads}
        """

rule diamond_denit:
    input:
         "results/pyrodigal/{sample}/cds_proteins.faa"
    output:
        "results/diamond_contigs/{sample}/{sample}_uniref100_denit.tsv"
    params:
        db_loc = "/projects/p31629/resources/denit_genes/uniref100_denit"
    threads: 10
    conda: "../envs/diamond.yaml"
    resources:
        mem="50G",
        time="05:00:00"
    shell:
        """
        diamond blastp -d {params.db_loc} -q {input} -o {output} \
        --very-sensitive --iterate --verbose --threads {threads} 
        """