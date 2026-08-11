import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = REPO_ROOT / "python"
CHAPTERS = (
    "01_Modeling_a_Single_Neuron",
    "03_The_Classical_HH_ODEs",
    "04_Numerical_Solution_of_HH_ODEs",
    "05_The_Simple_Model_of_Neurons_in_Rodent_Brains",
    "07_Linear_Integrate_and_Fire_(LIF)_Neurons",
    "08_Quadratic_Integrate_and_Fire_(QIF)_and_Theta_Neurons",
    "10_The_Slow_Fast_Phase_Plane",
    "11_The_Saddle_Node_Bifurcation",
    "12_Two_Dimensional_Bifurcation_Analysis",
    "13_Hopf_Bifurcations",
    "14_Model_Neurons_of_Bifurcation_Type_2",
    "15_Canard_Explosions",
    "16_Model_Neurons_of_Bifurcation_Type_3",
    "17_Frequency_Current_Curves",
    "18_Bistability_Resulting_from_Rebound_Firing",
    "19_Bursting",
    "20_Chemical_Synapses",
    "21_Gap_Junctions",
    "22_A_Wilson_Cowan_Model_of_an_Oscillatory_E-I_Network",
    "23_Entrainment_by_Excitatory_Input_Pulses",
    "24_Synchronization_by_Fast_Recurrent_Excitation",
    "25_Phase_Response_Curves_(PRCs)",
    "26_Phase_Locking_of_Two_Oscillators",
    "27_Phase_Locking_with_Delays",
    "28_Weakly_Coupled_Oscillators",
    "29_Stability_of_the_Synchronous_State",
    "30_The_PING_Model_of_Gamma_Rhythms",
    "31_ING_Rhythms",
    "32_M_Current_PING_and_Poisson_PING",
    "33_M_Current_PING_and_PINB",
    "34_Nested_Gamma_Theta_Rhythms",
    "35_Periodic_Inhibition",
    "36_F_I_Curves_Pulsed_Excitation",
    "37_Thresholding_in_PING",
    "38_Gamma_Coherence",
    "39_Short-Term_Depression_and_Facilitation",
    "40_Spike_Timing-Dependent_Plasticity(STDP)",
)
REQUIRED_HEADINGS = (
    "## Overview",
    "## Core ideas",
    "## Essential model",
    "## Code examples",
    "## What to look for",
    "## Suggested order",
    "## Prerequisites and related chapters",
    "## Running the examples",
)
LINK_RE = re.compile(r"\[[^\]]+\]\((?:<([^>]+)>|([^)]+))\)")


def chapter_number(chapter: str) -> str:
    return chapter.split("_", 1)[0]


def notebook_path(chapter: str) -> Path:
    return PYTHON_ROOT / f"chapter{chapter_number(chapter)}.ipynb"


def is_converted(chapter: str) -> bool:
    """True once a chapter's main.py scripts have been replaced by a
    chapterNN.ipynb notebook -- its guide then moves from
    python/<chapter>/README.md to a flat python/chapterNN.md next to the
    notebook, and the python/<chapter>/ directory is deleted entirely.
    """
    return notebook_path(chapter).is_file()


def guide_path(chapter: str) -> Path:
    if is_converted(chapter):
        return PYTHON_ROOT / f"chapter{chapter_number(chapter)}.md"
    return PYTHON_ROOT / chapter / "README.md"


def example_directories(chapter: str) -> list[Path]:
    if is_converted(chapter):
        return []
    return sorted(
        path
        for path in (PYTHON_ROOT / chapter).iterdir()
        if path.is_dir()
        and path.name != "__pycache__"
        and not path.name.startswith(".")
    )


@pytest.mark.parametrize("chapter", CHAPTERS)
def test_chapter_guide_has_required_structure(chapter: str) -> None:
    path = guide_path(chapter)
    assert path.is_file(), f"Missing chapter guide: {path.relative_to(REPO_ROOT)}"
    text = path.read_text(encoding="utf-8")
    missing = [heading for heading in REQUIRED_HEADINGS if heading not in text]
    assert not missing, f"{path.relative_to(REPO_ROOT)} missing {missing}"


@pytest.mark.parametrize("chapter", CHAPTERS)
def test_chapter_guide_covers_every_example(chapter: str) -> None:
    path = guide_path(chapter)
    text = path.read_text(encoding="utf-8")
    for example in example_directories(chapter):
        assert f"`{example.name}`" in text, f"Missing description for {example.name}"
        assert f"({example.name}/)" in text, f"Missing link for {example.name}"


@pytest.mark.parametrize("chapter", CHAPTERS)
def test_chapter_guide_relative_links_resolve(chapter: str) -> None:
    path = guide_path(chapter)
    text = path.read_text(encoding="utf-8")
    for match in LINK_RE.finditer(text):
        target = unquote(match.group(1) or match.group(2)).split("#", 1)[0]
        if not target or urlsplit(target).scheme or target.startswith("mailto:"):
            continue
        assert (path.parent / target).resolve().exists(), (
            f"Broken link in {path.relative_to(REPO_ROOT)}: {target}"
        )


def test_inventory_matches_all_python_chapter_directories() -> None:
    actual = {
        path.name
        for path in PYTHON_ROOT.iterdir()
        if path.is_dir() and re.match(r"^\d{2}_", path.name)
    }
    non_converted = {chapter for chapter in CHAPTERS if not is_converted(chapter)}
    assert non_converted == actual - {"09_Spike_Frequency_Adaptation"}
    for chapter in CHAPTERS:
        if is_converted(chapter):
            assert not (PYTHON_ROOT / chapter).exists(), (
                f"{chapter}/ should be deleted once converted to {notebook_path(chapter).name}"
            )
            assert guide_path(chapter).is_file()
    assert len(CHAPTERS) == 37


def test_python_index_links_every_chapter() -> None:
    index = PYTHON_ROOT / "README.md"
    assert index.is_file()
    text = index.read_text(encoding="utf-8")
    for chapter in CHAPTERS:
        if is_converted(chapter):
            assert f"({guide_path(chapter).name})" in text
        else:
            encoded = chapter.replace("(", "%28").replace(")", "%29")
            assert f"({encoded}/README.md)" in text
    assert "Chapters 2 and 6" in text
    assert "Chapter 9" in text
    assert "(09_Spike_Frequency_Adaptation/README.md)" not in text


def test_root_readme_links_python_guide() -> None:
    text = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    assert "(python/README.md)" in text


def test_all_documentation_relative_links_resolve() -> None:
    paths = [PYTHON_ROOT / "README.md", REPO_ROOT / "README.md"]
    paths.extend(guide_path(chapter) for chapter in CHAPTERS)
    for path in paths:
        text = path.read_text(encoding="utf-8")
        for match in LINK_RE.finditer(text):
            target = unquote(match.group(1) or match.group(2)).split("#", 1)[0]
            if not target or urlsplit(target).scheme or target.startswith("mailto:"):
                continue
            assert (path.parent / target).resolve().exists(), (
                f"Broken link in {path.relative_to(REPO_ROOT)}: {target}"
            )
