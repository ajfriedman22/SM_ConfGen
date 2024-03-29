"""
Unit and regression test for the SM_ConfGen package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import SM_ConfGen


def test_SM_ConfGen_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "SM_ConfGen" in sys.modules
