import pandas as pd


def pandas_index_intersect(*argc) -> pd.Index:
    """
    Computes pandas index intersection for arbitrary set of pandas tables
    :param argc: pandas tables or series
    :return: pandas index
    """
    en = enumerate(argc)
    ret = next(en)[1].index
    for _, series in en:
        ret = ret.intersection(series.index)
    return ret
