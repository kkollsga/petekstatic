from __future__ import annotations


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "perf: release-gated wall-clock regression budgets; run explicitly for "
        "static property workflow performance validation.",
    )
