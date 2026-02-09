"""Project modules.

This repo uses a src-layout. For local runs (without installation), add the repo root to PYTHONPATH,
or add the `src/` directory to `sys.path` in notebooks/scripts.

Main entry points:
- DecisionMakingNet (network.py)
- Tool utilities (tools.py)
- Motor readout utilities (motor.py)
"""

from .network import DecisionMakingNet  # noqa: F401
from .tools import Tool  # noqa: F401
