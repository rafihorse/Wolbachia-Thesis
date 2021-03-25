IDS = range(1, 1673)

rule all:
	input:
		expand("prot_{file}/paml_results.out", file=IDS)

rule macse_alignseqs:
        input:
                "all_prot.part-{file}.fasta"
        output:
                NT = "prot_{file}/all_prot.part-{file}_aligned_NT.fasta",
		AA = "prot_{file}/all_prot.part-{file}_aligned_AA.fasta"
        shell:
                "java -jar ~/bin/macse_v2.03.jar -prog alignSequences -seq {input} "
		"-out_NT {output.NT} -out_AA {output.AA}"

rule trim_alignments:
	input:
		"prot_{file}/all_prot.part-{file}_aligned_NT.fasta"
	output:
		"prot_{file}/all_prot.part-{file}.trimal"
	shell:
		"trimal -gt .85 -cons 60 -in {input} -out {output}"

rule macse_realignseqs:
        input:
                "prot_{file}/all_prot.part-{file}.trimal"
	output:
		NT = "prot_{file}/all_prot.part-{file}_refined_NT.fasta",
                AA = "prot_{file}/all_prot.part-{file}_refigned_AA.fasta"
	shell:
                "java -jar ~/bin/macse_v2.03.jar -prog refineAlignment -align {input} "
		"-out_NT {output.NT} -out_AA {output.AA}"

rule macse_exportaln:
        input:
                "prot_{file}/all_prot.part-{file}_refined_NT.fasta"
        output:
                NT = "prot_{file}/all_prot.part-{file}_noFSstop_NT.fasta",
        	AA = "prot_{file}/all_prot.part-{file}_noFSstop_AA.fasta"
	shell:
                "java -jar ~/bin/macse_v2.03.jar -prog exportAlignment -align {input} -codonForInternalStop NNN "
                "-codonForFinalStop --- -codonForInternalFS --- -charForRemainingFS - "
		"-out_NT {output.NT} -out_AA {output.AA}"

rule generate_phylip:
	input:
		"prot_{file}/all_prot.part-{file}_noFSstop_NT.fasta"
	output:
		"prot_{file}/all_prot.part-{file}.phy"
	shell:
		"/home/smnieves/bin/WritePhylip.py {input} {output}"

rule generate_tree:
	input:
		"prot_{file}/all_prot.part-{file}.phy"
	output:
		"prot_{file}/all_prot.part-{file}.phy.treefile"
	shell:
		"iqtree -s {input} -m GTR+G"

rule run_paml:
	input:
		seqfile="prot_{file}/all_prot.part-{file}.phy",
                treefile="prot_{file}/all_prot.part-{file}.phy.treefile"
	params:
		wdir="prot_{file}",
		model=0
	output:
		"prot_{file}/paml_results.out"
	shell:
		"PAML.py {params.wdir} {input.seqfile} {input.treefile} {params.model} {output}"
