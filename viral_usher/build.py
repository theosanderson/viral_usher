# 'build' subcommand: build a tree according to parameters in config file

import datetime
import docker
import filecmp
import os
import shutil
import subprocess
import sys
from .config import parse_config, write_config

docker_image = 'angiehinrichs/viral_usher'
docker_platform = 'linux/amd64'
docker_workdir = '/data'


def check_docker_command():
    try:
        subprocess.run(['docker', '--help'], check=True, capture_output=True)
    except:
        print("\nUnable to run 'docker --help' -- please install Docker from https://www.docker.com/ and try again.\n", file=sys.stderr)
        return False
    try:
        subprocess.run(['docker', 'images'], check=True, capture_output=True)
    except:
        print("\nUnable to run 'docker --images' -- please make sure the docker daemon is running (e.g. on a Mac, start the Docker app)\n", file=sys.stderr)
        return False
    return True


def image_created_within_last_day(image: docker.models.images.Image):
    """Return True if the Docker image object was created within the past 24 hours"""
    created_at_str = image.attrs["Created"]
    created_at = datetime.datetime.strptime(created_at_str[:26], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=datetime.timezone.utc)
    now = datetime.datetime.now(datetime.timezone.utc)
    return (now - created_at) <= datetime.timedelta(days=1)


def maybe_pull_docker_image(client: docker.DockerClient, image_name: str):
    """If image does not already exist locally, or is more than a day old, pull it."""
    need_pull = False
    try:
        image = client.images.get(image_name)
    except docker.errors.ImageNotFound:
        print(f"Docker image {image_name} is not present, pulling it from DockerHub...")
        need_pull = True
    if image:
        if not image_created_within_last_day(image):
            print(f"Pulling latest version of docker image {image_name} from DockerHub...")
            need_pull = True
    if need_pull:
        client.images.pull(image_name)


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


def rewrite_config(config, workdir):
    """Make a new timestamped config file in workdir that contains the given config settings"""
    now = datetime.datetime.now()
    new_name = 'local_config.' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.toml'
    new_path = workdir + '/' + new_name
    write_config(config, new_path)
    return new_name


def handle_build(args):
    if not check_docker_command():
        sys.exit(1)
    print(f"Running pipeline with config file: {args.config}")
    config = parse_config(args.config)
    workdir = config['workdir']
    user_id = os.getuid()
    group_id = os.getgid()
    docker_client = docker.from_env()
    # The config path passed to the script in the docker image must be relative to workdir
    # because that is the only directory visible to the script.
    config_rel = maybe_copy_to_workdir(args.config, workdir)
    # If extra_fasta is given, make sure it will be available in workdir and that
    # the config used by the docker image script uses its relative path in workdir.
    config_changed = False
    extra_fasta = config.get('extra_fasta', '')
    if extra_fasta:
        extra_fasta_rel = maybe_copy_to_workdir(extra_fasta, workdir)
        if extra_fasta_rel != extra_fasta:
            config['extra_fasta'] = extra_fasta_rel
            config_changed = True
    if config_changed:
        config_rel = rewrite_config(config, workdir)

    print(f"Running docker image {docker_image}")
    try:
        maybe_pull_docker_image(docker_client, docker_image)
        container = docker_client.containers.run(
            docker_image,
            platform=docker_platform,
            remove=True,
            network_mode='host',
            user=':'.join([str(user_id), str(group_id)]),
            volumes=[':'.join([workdir, docker_workdir])],
            command=["viral_usher_build", "--config", config_rel],
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
    if os.path.exists(f"{workdir}/tree.jsonl.gz"):
        print(f"Success -- you can view {workdir}/tree.jsonl.gz using https://taxonium.org/ now.")
    else:
        print(f"The docker command ran successfully but the expected result file {workdir}/tree.jsonl.gz is not there.  Please file an issue in GitHub: https://github.com/AngieHinrichs/viral_usher/issues/new", file=sys.stderr)
        sys.exit(1)
