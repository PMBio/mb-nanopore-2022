import matplotlib.pyplot as plt
from nanoepitools.alignment_quality import AlignmentTypeTotals


def plot_alignment_type_totals(totals: AlignmentTypeTotals):
    plt.pie(
        [
            totals.total_chimeric,
            totals.total_nonchimeric_multi,
            totals.total_inverse_repeats,
            totals.total_only_supp,
            totals.total_single_alignment,
        ],
        labels=[
            "Chimeric same strand",
            "Multiple primary alignments",
            "Chimeric different strand",
            "Only Supplementary",
            "Single primary alignment",
        ],
    )
