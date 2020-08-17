import os

seq_files = os.listdir("seqs")

filenames_stripped = [os.path.splitext(file)[0] for file in seq_files]

rule all:
    input:
        expand("filtered/{seq}_filtered.fasta",seq=filenames_stripped)

rule zorro:
    input:
        "seqs/{seq}.fasta"
    output:
        "masks/{seq}.mask"
    threads:
        1
    shell:
        "zorro -sample {input} > {output}"
        
rule filter:
    input:
        seq="seqs/{seq}.fasta",
        mask="masks/{seq}.mask"
    output:
        "filtered/{seq}_filtered.fasta"
    threads:
        1
    shell:
        """
        ~/bin/zorro.py -s {input.seq} -m {input.mask} -t 5 -o {output} -c
        """
