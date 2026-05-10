"""Smoke tests: package imports and exposes expected submodules."""

import mag_rescue


def test_version():
    assert isinstance(mag_rescue.__version__, str)
    assert mag_rescue.__version__


def test_submodules():
    assert hasattr(mag_rescue, "pp")
    assert hasattr(mag_rescue, "tl")
    assert hasattr(mag_rescue, "pl")
