# regfans/__init__.py
from importlib.metadata import PackageNotFoundError, version

from .vectorconfig import VectorConfiguration

try:
    __version__ = version("regfans")
except PackageNotFoundError:   # running from a source tree without an install
    __version__ = "0.0.0+unknown"

__all__ = ["VectorConfiguration", "__version__"]
