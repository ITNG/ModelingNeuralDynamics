"""Git clean/diff filter: retain notebook content without runtime metadata."""

import argparse
import json
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", help="Read a file instead of stdin (Git diff).")
    parser.add_argument("--install", action="store_true", help="Install the local Git filter.")
    args = parser.parse_args()
    if args.install:
        command = ".venv/bin/python scripts/strip_notebook.py"
        for key, value in {
            "filter.nbstripout.clean": command,
            "filter.nbstripout.smudge": "cat",
            "filter.nbstripout.required": "true",
            "diff.ipynb.textconv": command,
        }.items():
            subprocess.run(["git", "config", "--local", key, value], check=True)
        return

    if args.path:
        with open(args.path, encoding="utf-8") as source:
            notebook = json.load(source)
    else:
        notebook = json.load(sys.stdin)
    if notebook.get("nbformat") != 4:
        raise ValueError("The notebook Git filter requires notebook format 4.")
    notebook["metadata"] = {}
    for index, cell in enumerate(notebook["cells"]):
        cell["metadata"] = {}
        if cell["cell_type"] == "code":
            cell["execution_count"] = None
            cell["outputs"] = []
        if "id" in cell:
            cell["id"] = str(index)
    json.dump(notebook, sys.stdout, ensure_ascii=False, indent=1, sort_keys=True)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
