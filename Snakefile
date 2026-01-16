SIMULATIONS = range(2)
BUCKET = "carbonplan-scratch/ew-workflows"

rule all:
    input:
        expand(f"outputs/sim_{{sim_id}}/result.txt", sim_id=SIMULATIONS)

rule run_simulation:
    output:
        f"outputs/sim_{{sim_id}}/result.txt"
    resources:
        vcpus=2,
        mem_mb=7700
    shell:
        """
	/srv/conda/envs/notebook/bin/pip install snakemake snakemake-storage-plugin-s3
        git clone --branch your-branch-name https://github.com/carbonplan/ew-workflows.git /app
        cd /app
        /srv/conda/envs/notebook/bin/pip install .
        echo "result from sim {wildcards.sim_id}" > {output}
        """
