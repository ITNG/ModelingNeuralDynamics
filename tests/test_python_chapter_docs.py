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


def guide_path(chapter: str) -> Path:
    return PYTHON_ROOT / chapter / "README.md"


def example_directories(chapter: str) -> list[Path]:
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
