# import pandas as pd
# params_df = pd.read_csv("simulations.csv")
# TK .to_dict?

SIMULATIONS = {'my-second-test':{'row_number': 1, 'dustrate': 0, 'dustrad': 5, 'spinrun': 'site_311a_pr9_spintuneup4',  }}

rule all:
    input:
        expand("outputs/sim_{sim_id}/result.txt", sim_id=SIMULATIONS.keys())

rule run_simulation:
    output:
        'outputs/sim_{sim_id}/result.txt'
    params:
        row_number=lambda x: SIMULATIONS[x.sim_id]['row_number'],
        dustrate=lambda x: SIMULATIONS[x.sim_id]['dustrate'],
        dustrad=lambda x: SIMULATIONS[x.sim_id]['dustrad'],
        spinrun=lambda x: SIMULATIONS[x.sim_id]['spinrun'],
    resources:
        vcpus=2,
    shell:
        r""" 
	 aws batch submit-job \
	  --job-name "scepter-test" \
	  --job-queue arn:aws:batch:us-west-2:631969445205:job-queue/ew-workflow-compute-env \
	  --job-definition test-ec2-job \
	  --container-overrides '{{"command": ["bash", "-c", "git clone --branch aws-batch-demo https://github.com/carbonplan/ew-workflows.git /home/jovyan/ew-workflows && pip install -e /home/jovyan/ew-workflows && python3 /home/jovyan/SCEPTER/rock_app_singleRun-szn.py --modeldir /home/jovyan/SCEPTER/ --outdir /home/jovyan/SCEPTER/scepter_output/ --default_dict singlerun_default --row_number {params.row_number} --dustrate {params.dustrate} --dustrad {params.dustrad} --spinrun {params.spinrun} --spindir s3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/ --duration 1 --dustsp gbas --singlerun_seasonality False --aws_save move --aws_bucket s3://carbonplan-scratch/SCEPTER/scepter_output_scratch --scepter_exec_name scepter --newrun_id {wildcards.sim_id}"]}}' > /dev/null
        """ 
