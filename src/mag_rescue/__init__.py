"""mag-rescue: ARIBA-driven virulence/AMR profiling on Klebsiella short-read genomes."""

from importlib.metadata import version

from . import pl, pp, tl

__all__ = ["pl", "pp", "tl"]
__version__ = version("mag-rescue")
