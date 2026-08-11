from pathlib import Path

import pytest


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Group Brian2 notebook tests without editing every test module."""
    for item in items:
        if "_brian_" in Path(str(item.path)).name:
            item.add_marker(pytest.mark.brian2)
