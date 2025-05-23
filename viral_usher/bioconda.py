import subprocess

def install_tools(tools, env_name="viral_usher-env"):
    try:
        subprocess.run(["conda", "--version"], check=True)
    except FileNotFoundError:
        raise RuntimeError("Conda not found in PATH.")

    # Create environment if not exists
    subprocess.run(["conda", "create", "-n", env_name, "-y"], check=True)

    # Install tools
    install_cmd = ["conda", "install", "-n", env_name, "-c", "bioconda", "-y"] + tools
    subprocess.run(install_cmd, check=True)
    