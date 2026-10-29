"""Option-pricing research package with stable historical module names."""
from pathlib import Path

_package = Path(__file__).parent
__path__.extend([str(_package / "code" / "tests"), str(_package / "code")])
