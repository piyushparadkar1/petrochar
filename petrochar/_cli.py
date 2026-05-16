"""
Command-line entry point for petrochar.

Invoked by the `petrochar` console script (defined in pyproject.toml
[project.scripts]).  Launches the Streamlit app via subprocess so that
`petrochar` at the shell prompt is equivalent to `streamlit run app.py`.

Usage
-----
    petrochar              # launch with default port 8501
    petrochar --port 8080  # launch on a custom port

Notes
-----
The app.py file is located at the repository root.  The entry point
locates it relative to this module's package directory.
"""

import pathlib
import subprocess
import sys


def main() -> None:
    """Launch the petrochar Streamlit application."""
    app_path = pathlib.Path(__file__).parent.parent / "app.py"
    if not app_path.exists():
        print(
            f"petrochar: could not find app.py at {app_path}.\n"
            "Run `petrochar` from the repository root directory, or reinstall "
            "the package with `pip install -e .`.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Pass any extra arguments through to streamlit (e.g., --port, --server.*)
    extra_args = sys.argv[1:]

    cmd = [sys.executable, "-m", "streamlit", "run", str(app_path)] + extra_args
    try:
        subprocess.run(cmd, check=True)
    except KeyboardInterrupt:
        pass
    except subprocess.CalledProcessError as exc:
        sys.exit(exc.returncode)


if __name__ == "__main__":
    main()
