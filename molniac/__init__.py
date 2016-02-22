try:
    # python 3 ?
    from molniac.molniac import load
except ImportError:
    # python 2 ?
    from molniac import load

__all__ = []
