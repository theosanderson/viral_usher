import subprocess

already_checked = False

def check_conda():
    """Check if conda is installed and available in PATH."""
    global already_checked
    if already_checked:
        return True
    already_checked = True
    try:
        subprocess.run(["conda", "--version"], check=True)
    except FileNotFoundError:
        raise RuntimeError("Conda not found in PATH. Please install Anaconda or Miniconda.")
    except subprocess.CalledProcessError:
        raise RuntimeError("Conda command failed. Please check your conda installation.")
    return True

def check_env_name(env_name):
    if not env_name or not isinstance(env_name, str):
        raise ValueError("Environment name arg must be a non-empty string.")

def check_env(env_name):
    """Check if the conda environment exists."""
    check_env_name(env_name)
    check_conda()
    try:
        result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True, check=True)
        envs = result.stdout.splitlines()
        for line in envs:
            if env_name in line:
                return True
    except subprocess.CalledProcessError:
        raise RuntimeError("Failed to list conda environments.")
    return False

def create_env(env_name):
    """Create a new conda environment."""
    check_env_name(env_name)
    check_conda()
    if not check_env(env_name):
        subprocess.run(["conda", "create", "-n", env_name, "-y"], check=True)
    else:
        print(f"Environment '{env_name}' already exists.")

def activate_env(env_name):
    """Activate the conda environment."""
    check_env_name(env_name)
    check_conda()
    try:
        subprocess.run(["conda", "activate", env_name], shell=True, check=True)
    except subprocess.CalledProcessError:
        raise RuntimeError(f"Failed to activate conda environment '{env_name}'. Please activate it manually.")

def install_tools(tools, env_name):
    check_env_name(env_name)
    check_conda()
    # Create environment if not exists
    if not check_env(env_name):
        create_env(env_name)
    activate_env(env_name)

    install_cmd = ["conda", "install", "-n", env_name, "-c", "bioconda", "-y"] + tools
    try:
        subprocess.run(install_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error installing tools: {e}")
        raise
