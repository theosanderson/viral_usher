# 'build' subcommand: build a tree according to parameters in config file

import datetime
import docker
import filecmp
import os
import shutil
import subprocess
import sys
from . import config

docker_image = 'angiehinrichs/viral_usher'
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
        return False, "Unable to run 'docker --images' -- please make sure the docker daemon is running (e.g. on a Mac, start the Docker app)"
    return True, None


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
        image = None
        need_pull = True
    if image:
        if not image_created_within_last_day(image):
            print(f"Pulling latest version of docker image {image_name} from DockerHub...")
            need_pull = True
    if need_pull:
        client.images.pull(image_name, platform=docker_platform)


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


def run_in_docker(workdir, config_rel):
    """Start up a docker container that runs the pipeline."""
    docker_client = docker.from_env()
    user_id = os.getuid()
    group_id = os.getgid()
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
    run_in_docker(workdir, config_rel)
    if os.path.exists(f"{workdir}/tree.jsonl.gz"):
        print(f"Success -- you can view {workdir}/tree.jsonl.gz using https://taxonium.org/ now.")
    else:
        print(f"The pipeline exited without an error status but the expected result file {workdir}/tree.jsonl.gz is not there.  Please file an issue in GitHub: https://github.com/AngieHinrichs/viral_usher/issues/new", file=sys.stderr)
        sys.exit(1)
