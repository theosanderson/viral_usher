# 'build' subcommand
import docker
import os
import shutil
import subprocess
import sys
from .config import parse_config

docker_image = 'angiehinrichs/viral_usher'
docker_platform = 'linux/amd64'
docker_workdir = '/data'

def check_docker_command():
    try:
        result = subprocess.run(['docker', '--help'], check=True, capture_output=True)
    except:
        print(f"\nUnable to run 'docker --help' -- please install Docker from https://www.docker.com/ and try again.\n", file=sys.stderr)
        return False
    try:
        result = subprocess.run(['docker', 'images'], check=True, capture_output=True)
    except:
        print(f"\nUnable to run 'docker --images' -- please make sure the docker daemon is running (e.g. on a Mac, start the Docker app)\n", file=sys.stderr)
        return False
    return True

def handle_build(args):
    if not check_docker_command():
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

    print(f"Running docker image {docker_image}")
    try:
        docker_client.containers.run(docker_image,
                                    platform=docker_platform,
                                    remove=True,
                                    network_mode='host',
                                    user=':'.join([str(user_id), str(group_id)]),
                                    volumes=[ ':'.join([workdir, docker_workdir]) ],
                                    command=["viral_usher_build", "--config", config_basename])
    except (docker.errors.ContainerError) as e:
        print(f"docker container {docker_image} failed:\n{e}")
        sys.exit(1)
    print(f"Success -- you can view {workdir}/tree.json.gz in taxonium now.")
    #TODO check results, maybe save a log file?