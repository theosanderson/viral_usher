# 'build' subcommand: build a tree according to parameters in config file

import datetime
import docker
import filecmp
import importlib.metadata
import os
import shutil
import subprocess
import sys
from . import config

docker_platform = 'linux/amd64'
docker_workdir = '/data'


def check_docker_command():
    """Make sure we can run the docker command and connect to the docker server.  If there's a problem, return False with
    an error message; otherwise return True with no message."""
    try:
        subprocess.run(['docker', '--help'], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False, "Unable to run 'docker --help' -- please install Docker from https://www.docker.com/ and try again."
    try:
        subprocess.run(['docker', 'images'], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        return False, "Unable to run 'docker images' -- please make sure the docker daemon is running (e.g. on a Mac, start the Docker app)"
    return True, None


def pull_docker_image(client: docker.DockerClient, image_name: str):
    """Pull the docker image corresponding to this version of the viral_usher python package unless it is already present."""
    try:
        client.images.pull(image_name, platform=docker_platform)
    except (docker.errors.ContainerError) as e:
        print(f"Failed to pull docker container {image_name}:\n{e}", file=sys.stderr)
        sys.exit(1)
    except (TimeoutError):
        print("Timeout error while trying to pull docker image; maybe try again later?")
        sys.exit(1)


def get_docker_image(client: docker.DockerClient, args_docker_image: str):
    if args_docker_image:
        try:
            client.images.get(args_docker_image)
            print(f"Using docker image {args_docker_image}")
        except docker.errors.ImageNotFound:
            print(f"Docker image '{args_docker_image}' not found.  Please try again with a different --docker_image value.", sys.stderr)
        return args_docker_image
    else:
        viral_usher_version = importlib.metadata.version('viral_usher')
        image_name = config.DEFAULT_DOCKER_IMAGE + ':v' + viral_usher_version
        try:
            client.images.get(image_name)
            print(f"Using docker image {image_name}")
        except docker.errors.ImageNotFound:
            print(f"Docker image {image_name} is not present, pulling it from DockerHub...")
            pull_docker_image(client, image_name)
        return image_name


def maybe_copy_to_workdir(filepath, workdir):
    """Return filepath's relative path within workdir, copying the file into workdir if necessary."""
    basename = os.path.basename(filepath)
    dirname = os.path.dirname(filepath)
    dirname_abs = os.path.abspath(dirname)
    workdir_abs = os.path.abspath(workdir)
    if dirname_abs == workdir_abs:
        return basename
    else:
        if dirname_abs.startswith(workdir_abs + '/'):
            relpath = dirname_abs.removeprefix(workdir_abs + '/')
            relname = relpath + '/' + basename
            return relname
        else:
            workdir_file = workdir + '/' + basename
            if os.path.exists(workdir_file):
                if not filecmp.cmp(filepath, workdir_file, shallow=False):
                    print(f"File {filepath} is specified, but workdir file {workdir_file} already exists and is not the same", file=sys.stderr)
                    sys.exit(1)
            else:
                shutil.copy(filepath, workdir)
            return basename


def rewrite_config(config_contents, workdir):
    """Make a new timestamped config file in workdir that contains the given config settings"""
    now = datetime.datetime.now()
    new_name = 'local_config.' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.toml'
    new_path = workdir + '/' + new_name
    config.write_config(config_contents, new_path)
    return new_name


def localize_config(args_config, config_contents):
    """When running in docker, only the workdir will be available (as the current directory).
    So absolute paths in the config need to be converted to relative paths, and if input files
    are not already in the workdir, then they need to be copied there."""
    workdir = config_contents['workdir']
    # If extra_fasta is given, make sure it will be available in workdir and that
    # the config used by the docker image script uses its relative path in workdir.
    config_changed = False
    extra_fasta = config_contents.get('extra_fasta', '')
    if extra_fasta:
        extra_fasta_rel = maybe_copy_to_workdir(extra_fasta, workdir)
        if extra_fasta_rel != extra_fasta:
            config_contents['extra_fasta'] = extra_fasta_rel
            config_changed = True
    if config_changed:
        config_rel = rewrite_config(config_contents, workdir)
    else:
        config_rel = maybe_copy_to_workdir(args_config, workdir)
    return config_rel


def run_in_docker(workdir, config_rel, args_docker_image, args_update):
    """Start up a docker container that runs the pipeline."""
    docker_client = docker.from_env()
    user_id = os.getuid()
    group_id = os.getgid()
    docker_image = get_docker_image(docker_client, args_docker_image)
    command = ["viral_usher_build", "--config", config_rel]
    if args_update:
        command.append("--update")
    try:
        container = docker_client.containers.run(
            docker_image,
            platform=docker_platform,
            remove=True,
            network_mode='host',
            user=':'.join([str(user_id), str(group_id)]),
            volumes=[':'.join([workdir, docker_workdir])],
            command=command,
            detach=True,  # Run container in background so we can stream logs
            stdout=True,
            stderr=True,
            tty=True
        )
        # Stream logs in real time
        for line in container.logs(stream=True):
            print(line.decode(), end='')
        exit_code = container.wait()['StatusCode']
        if exit_code != 0:
            print(f"docker container {docker_image} failed with exit code {exit_code}", file=sys.stderr)
            sys.exit(exit_code)
    except (docker.errors.ContainerError) as e:
        print(f"docker container {docker_image} failed:\n{e}", file=sys.stderr)
        sys.exit(1)
    except (TimeoutError):
        print("Timeout error while trying to run docker; maybe try again later?")
        sys.exit(1)


def handle_build(args):
    """Set up input files in workdir, run pipeline in docker and check final results."""
    ok, error = check_docker_command()
    if not ok:
        print(f"\n{error}\n", file=sys.stderr)
        sys.exit(1)
    print(f"Running pipeline with config file: {args.config}")
    config_contents = config.parse_config(args.config)
    config_rel = localize_config(args.config, config_contents)
    workdir = config_contents['workdir']
    run_in_docker(workdir, config_rel, args.docker_image, args.update)
    if os.path.exists(f"{workdir}/tree.jsonl.gz"):
        print(f"Success -- you can view {workdir}/tree.jsonl.gz using https://taxonium.org/ now.")
    else:
        print(f"The pipeline exited without an error status but the expected result file {workdir}/tree.jsonl.gz is not there.  Please file an issue in GitHub: https://github.com/AngieHinrichs/viral_usher/issues/new", file=sys.stderr)
        sys.exit(1)
