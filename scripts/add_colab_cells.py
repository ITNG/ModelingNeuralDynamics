"""One-off (and reusable-for-new-chapters) script: prepend a Colab badge
and a Colab-only install cell to a brian/*.ipynb notebook.

Usage: python scripts/add_colab_cells.py brian/chapter01.ipynb [more.ipynb ...]
"""
import json
import sys
from pathlib import Path

REPO = "ITNG/ModelingNeuralDynamics"


def badge_cell(notebook_path):
    rel = notebook_path.as_posix()
    url = f"https://colab.research.google.com/github/{REPO}/blob/main/{rel}"
    badge = f"[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({url})"
    return {"cell_type": "markdown", "metadata": {}, "source": [badge]}


def install_cell():
    # NOTE: deliberately avoids the `%pip install ...` IPython magic form.
    # tests/matlab_ref.py::load_notebook_as_module executes notebook code
    # cells via plain exec(compile(source, ...)) (no IPython kernel), and
    # compile() parses the whole cell regardless of the `if` guard below,
    # so an IPython magic there is a SyntaxError even though it never runs
    # outside Colab. subprocess + sys.executable is valid plain Python and
    # behaves identically in Colab.
    source = (
        "import subprocess\n"
        "import sys\n"
        "if \"google.colab\" in sys.modules:\n"
        "    subprocess.run([sys.executable, \"-m\", \"pip\", \"install\", \"-q\", \"modelingneuraldynamics\"], check=True)\n"
    )
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


def add_colab_cells(path):
    path = Path(path)
    nb = json.loads(path.read_text())
    already = any(
        "colab.research.google.com" in "".join(cell.get("source", []))
        for cell in nb["cells"]
    )
    if already:
        print(f"{path}: already has a Colab badge, skipping")
        return
    # Keep the notebook's original title cell first; insert the badge and
    # install cells right after it rather than before it.
    nb["cells"] = nb["cells"][:1] + [badge_cell(path), install_cell()] + nb["cells"][1:]
    if nb.get("nbformat_minor", 0) < 5:
        for cell in nb["cells"]:
            cell.pop("id", None)
    path.write_text(json.dumps(nb, indent=1))
    print(f"{path}: added Colab badge + install cell")


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        add_colab_cells(arg)
