import json
from plotsf.convert.clinvar import read_clinvar_tabular_file
from plotsf.genemodel.plot import plot_figure

if __name__ == "__main__":
    # this is a really basic test example
    with open("example/PTEN.json", "r") as handle:
        gene = json.load(handle)
    data = read_clinvar_tabular_file("example/PTEN_clinvar.txt")
    fig = plot_figure(data, gene)
    fig.savefig("example/PTEN_test.pdf")
