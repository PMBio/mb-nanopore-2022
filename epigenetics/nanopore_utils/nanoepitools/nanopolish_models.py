import pandas as pd


def read_nanopolish_model(filename):
    return pd.read_csv(
        filename,
        sep="\t",
        comment="#",
        names=["kmer", "level_mean", "level_stdv", "std", "std_stdv"],
    ).set_index("kmer")


def read_nanopolish_training_summary(filename):
    ret = pd.read_csv(
        filename,
        sep="\t",
        header=0,
        names=[
            "model",
            "kmer",
            "matches",
            "skips",
            "stays",
            "num_events",
            "was_trained",
            "level_mean",
            "level_stdv",
        ],
        dtype={"num_events": int},
    )
    return ret[
        ["kmer", "num_events", "was_trained", "level_mean", "level_stdv"]
    ].set_index("kmer")
