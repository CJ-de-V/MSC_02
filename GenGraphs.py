import numpy as np
from matplotlib import pyplot as plt
from matplotlib import scale as scl
import pandas as pd
import sys
from matplotlib.backends.backend_pdf import PdfPages
from itertools import cycle

# makes a cycle that lets us iterate over which colors to use
from distinctipy import distinctipy

colors = distinctipy.get_colors(10)
print(colors)
number = 1


# Pandas dataframe, x dependent & y independent variable, variable on legend and number of figure for plotting
# figure numbers have to be distinct
def loglogplot(datfram, Xdata, Ydata, legend, xerror=0, yerror=0):
    global number
    plt.figure(number)
    number += 1
    legset = datfram[legend].unique()

    numcolz = len(legset)  # number of different colors
    for i in range(len(legset)):
        ndata = df.loc[df[legend] == legset[i]]
        x = np.array(ndata[Xdata].tolist())
        y = np.array(ndata[Ydata].tolist())
        plt.loglog(x, y, 'o-', base=2, c=colors[i])

        if xerror != 0:
            plt.errorbar(x, y, xerr=np.array(ndata[xerror].tolist()), c=colors[i])
        if yerror != 0:
            plt.errorbar(x, y, yerr=np.array(ndata[yerror].tolist()), c=colors[i])

    plt.title("Log-Log plot of " + Xdata + " vs " + Ydata + " for various " + legend)
    plt.xlabel("Log(" + Xdata + ")")
    plt.ylabel("Log(" + Ydata + ")")
    plt.legend(legset, title=legend)


def save_multi_image(filename):
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()


df = pd.read_csv('finaldata.csv')

# plot ln Lp vs ln k
# print(df)
loglogplot(df, 'Kp', 'lp', 'Np', yerror='sd_lp')
# plot ln N vs ln Rg^2

loglogplot(df, 'Np', 'Rg', 'Kp', yerror='sd_Rg')

# plot ln R_g vs ln lp

#loglogplot(df, 'R_g', 'lp', 'k', xerror='R_g_sd', yerror='lp_sd')

# plot ln R_g vs ln lp

#loglogplot(df, 'R_g', 'lp', 'N', xerror='lp_sd', yerror='R_g_sd')

# save_multi_image("FinalGraphs/Results.pdf")

plt.show()
