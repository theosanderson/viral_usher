# 'build' subcommand
import docker
import os
import shutil
import subprocess
import sys
import tomllib

docker_image = 'angiehinrichs/viral_usher'
docker_platform = 'linux/amd64'
docker_workdir = '/data'

def check_docker_command():
    success = False
    try:
        result = subprocess.run(['docker', '--help'], check=True, capture_output=True)
        success = True
    except:
        pass
    return success

def read_config(config_path):
    with open(config_path, 'rb') as f:
        config = tomllib.load(f)
    return config

def check_refseq_acc(refseq_acc):
    if not refseq_acc or not isinstance(refseq_acc, str):
        raise ValueError("RefSeq accession must be a non-empty string.")
    if not refseq_acc.startswith("NC_") and not refseq_acc.startswith("NZ_"):
        raise ValueError("RefSeq accession must start with 'NC_' or 'NZ_'.")

def check_taxonomy_id(taxonomy_id):
    if not taxonomy_id or not isinstance(taxonomy_id, str):
        raise ValueError("Taxonomy ID must be a non-empty string.")
    if not taxonomy_id.isdigit():
        raise ValueError("Taxonomy ID must be a numeric string.")
    
def check_refseq_assembly(refseq_assembly):
    if not refseq_assembly or not isinstance(refseq_assembly, str):
        raise ValueError("RefSeq assembly must be a non-empty string.")
    if not refseq_assembly.startswith("GCA_") and not refseq_assembly.startswith("GCF_"):
        raise ValueError("RefSeq assembly must start with 'GCA_' or 'GCF_'.")
    
def parse_config(config_path):
    """Read the config file and validate its contents."""
    config = read_config(config_path)
    print(f"refseq_acc: {config['refseq_acc']}")
    print(f"taxonomy_id: {config['taxonomy_id']}")
    print(f"refseq_assembly: {config['refseq_assembly']}")
    try:
        check_refseq_acc(config['refseq_acc'])
        check_taxonomy_id(config['taxonomy_id'])
        check_refseq_assembly(config['refseq_assembly'])
    except (ValueError, FileNotFoundError) as e:
        print(f"Error in config file {config_path}: {e}", file=sys.stderr)
        raise
    return config

def handle_build(args):
    if not check_docker_command():
        print(f"docker command not working; please install Docker.")
        sys.exit(1)
    print(f"Running pipeline with config file: {args.config}")
    config = parse_config(args.config)
    workdir = config['workdir']
    user_id = os.getuid()
    group_id = os.getgid()
    docker_client = docker.from_env()
    config_dirname = os.path.dirname(args.config)
    config_basename = os.path.basename(args.config)
    if os.path.abspath(config_dirname) != os.path.abspath(workdir):
        #TODO make init write config file in workdir
        # -- or copy to a unique file name in workdir and use that in the docker run.
        print(f"Config file needs to be present in build directory ({workdir}); copying {args.config} into build directory.")
        shutil.copy(args.config, workdir)

    #TODO client.images.pull(docker_image)
    print(f"Running docker image {docker_image}")
    try:
        docker_client.containers.run(docker_image,
                                    platform=docker_platform,
                                    remove=True,
                                    user=':'.join([str(user_id), str(group_id)]),
                                    volumes=[ ':'.join([workdir, docker_workdir]) ],
                                    command=["viral_usher_build", "--config", config_basename])
    except (docker.errors.ContainerError) as e:
        print(f"docker container {docker_image} failed:\n{e}")
        sys.exit(1)
    print(f"Success -- you can view {workdir}/tree.json.gz in taxonium now.")
    #TODO check results, maybe save a log file?