import subprocess


def submit_sbatch(slurm_script):

    job_id, errs = subprocess.Popen("sbatch {}".format(slurm_script), # submit array
            shell="True",  stdout=subprocess.PIPE, universal_newlines=True).communicate()
    print(job_id, end='')
    job_id = job_id.split(" ")[-1].strip()

    return job_id

def submit_salloc(cmd, dependency_id, time='20:00'):
    subprocess.Popen('salloc -t {} --dependency={} {}'.format(time, dependency_id, cmd), shell=True)
