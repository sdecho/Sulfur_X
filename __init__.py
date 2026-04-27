from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("Sulfur_X")
except PackageNotFoundError:
    __version__ = "unknown"