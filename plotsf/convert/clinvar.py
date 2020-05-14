from typing import Dict
import re
import pandas as pd

# TODO: docstring and refactor name format
review_stars: Dict[str, int] = {
    "no assertion criteria provided": 0,
    "criteria provided, conflicting interpretations": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}

# TODO: docstring and refactor name format
clinsig_abbreviations: Dict[str, str] = {
    "Conflicting interpretations of pathogenicity": "C",
    "Benign": "B",
    "Likely benign": "LB",
    "Benign/Likely benign": "LB",
    "Uncertain significance": "VUS",
    "Likely pathogenic": "LP",
    "Pathogenic/Likely pathogenic": "LP",
    "Pathogenic": "P",
    "Risk factor": "R",
}


def read_clinvar_tabular_file(fname: str) -> pd.DataFrame:
    """Read a tabular file downloaded from ClinVar into a pandas DataFrame.

    Values in the original file will be converted to use a more concise internal format.
    For details of the conversion see :py:func:`convert_clinvar_tabular_row`.

    Parameters
    ----------
    fname: str
        File path for the tabular data file downloaded from the ClinVar website.

    Returns
    -------
    pd.DataFrame
        Data frame formatted for plotting.

    """
    try:
        raw_df = pd.read_csv(fname, sep="\t")
    except IOError as e:  # TODO: make this better
        raise e

    try:
        result_df = raw_df.apply(convert_clinvar_tabular_row, axis=1)
    except KeyError as e:
        raise e  # TODO: make this informative

    result_df.dropna(
        axis="index", how="any", inplace=True
    )  # drops any row that failed to convert

    return result_df


def convert_clinvar_tabular_row(row: pd.Series) -> pd.Series:
    """Convert a row from ClinVar tabular text output to internal format.

    # TODO: add details of conversion

    Parameters
    ----------
    row

    Returns
    -------

    """
    hgvs = row.loc["Name"]
    nt_match = re.search(r"c\.(-?\d+)[ACGT]>[ACGT]", hgvs)
    if nt_match is not None:
        nt_position = int(nt_match.group(1))
    else:
        nt_position = None
    aa_match = re.search(r"p\.[A-Z][a-z]{2}(-?\d+)[A-Z][a-z]{2}", hgvs)
    if aa_match is not None:
        aa_position = int(aa_match.group(1))
    else:
        aa_position = None
    gene = row.loc["Gene(s)"]
    clinsig = clinsig_abbreviations[
        row.loc["Clinical significance (Last reviewed)"].split("(")[0]
    ]
    stars = review_stars[row.loc["Review status"]]
    return pd.Series(
        {
            "gene": gene,
            "hgvs": hgvs,
            "nt_position": nt_position,
            "aa_position": aa_position,
            "clinsig": clinsig,
            "stars": stars,
        }
    )
