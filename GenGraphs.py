import numpy as np
from matplotlib import pyplot as plt
from matplotlib import scale as scl
import pandas as pd
import sys
from matplotlib.backends.backend_pdf import PdfPages
from itertools import cycle

# makes a cycle that lets us iterate over which colors to use
from distinctipy import distinctipy

colors = distinctipy.get_colors(25)
# print(colors)
number = 1


# Pandas dataframe, x dependent & y independent variable, variable on legend and number of figure for plotting
# figure numbers have to be distinct
def loglogplot(datfram, Xdata, Ydata, legend, legend2=0, xerror=0, yerror=0):
    global number
    plt.figure(number)
    number += 1
    legset = datfram[legend].unique()

    numcolz = len(legset)  # number of different colors
    for i in range(len(legset)):
        ndata = df.loc[df[legend] == legset[i]]
        x = np.array(ndata[Xdata].tolist())
        y = np.array(ndata[Ydata].tolist())
        plt.loglog(x, y, '.', base=2, c=colors[i])  # . formerly o-

        if xerror != 0:
            plt.errorbar(x, y, xerr=np.array(ndata[xerror].tolist()), c=colors[i])
        if yerror != 0:
            plt.errorbar(x, y, yerr=np.array(ndata[yerror].tolist()), c=colors[i])

    plt.title("Log-Log plot of " + Xdata + " vs " + Ydata + " for various " + legend)
    plt.xlabel("Log(" + Xdata + ")")
    plt.ylabel("Log(" + Ydata + ")")
    plt.legend(legset, title=legend)


# Identical to the above, except it is not a log log plot.

def notloglogplot(datfram, Xdata, Ydata, legend, legend2=0, xerror=0, yerror=0):
    global number
    plt.figure(number)
    number += 1
    legset = datfram[legend].unique()

    numcolz = len(legset)  # number of different colors
    for i in range(len(legset)):
        ndata = df.loc[df[legend] == legset[i]]
        x = np.array(ndata[Xdata].tolist())
        y = np.array(ndata[Ydata].tolist())
        plt.plot(x, y, '.', c=colors[i])  # . formerly o-

        if xerror != 0:
            plt.errorbar(x, y, xerr=np.array(ndata[xerror].tolist()), c=colors[i])
        if yerror != 0:
            plt.errorbar(x, y, yerr=np.array(ndata[yerror].tolist()), c=colors[i])

    plt.title("Plot of " + Xdata + " vs " + Ydata + " for various " + legend)
    plt.xlabel("" + Xdata + "")
    plt.ylabel("" + Ydata + "")
    plt.legend(legset, title=legend)


def save_multi_image(filename):
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()


def boundede2ecalc(k, n, l):
    return (np.sqrt(2 / k ** 2 + l ** 2*n
                    - (4 * l * np.sqrt(3 * n + k ** 2 * l ** 2 * n ** 2 / 4)) / (3 * k * np.pi)
                    - (4 * np.arctan(3 / (l * k * np.sqrt(3 * n + k ** 2 * l ** 2 * n ** 2 / 4)))) / (k ** 2 * np.pi)))


df = pd.read_csv('finaldata.csv')
# Np	Nt	Kp	Kt	lp	Rg	zd	lb	sd_lp	sd_Rg	sd_zd
# data columns for faster copypasting

loglogplot(df, 'Rg', 'lp', 'Kp')

loglogplot(df, 'Rg', 'lp', 'Np')

loglogplot(df, 'zd', 'Nt', 'Kt')

loglogplot(df, 'zd', 'Np', 'Nt')

loglogplot(df, 'lp', 'zd', 'Kp', )

loglogplot(df, 'Kt', 'Rg', 'Np', )

notloglogplot(df, 'Np', 'E2E', 'Kp', )

lb = df['lb'].mean()
Xn = np.linspace(1, 128, 128)
Kn = np.linspace(0,64,256)
lin1 =plt.plot(Xn, lb * Xn, label='rigid rod')  # Plots rigid rod's End to End distance
lin2 =plt.plot(Xn, boundede2ecalc(10, Xn, lb), label = 'strong spring bound gaussian')  # Plots bounded polymer's End to End distance
#plt.plot(Xn, boundede2ecalc(0.0000001, Xn, lb))  # Plots bounded polymer's End to End distance
lin3 =plt.plot(Xn, np.sqrt(Xn)*lb, label= 'free gaussian')  # Plots Free polymer's End to End distance
plt.legend()

notloglogplot(df, 'Nt', 'zd', 'Kt')
plt.plot(Xn, 3/np.sqrt(32)/Xn, label = 'Gaussian Chain bounded by spring')
# Plots bounded gaussian's average extension, does not really work that well since we don't have a good grasp on what
# the spring constant of these polymers are, first we need to map them to it, it is however definitely more related
# to the arguments of the FENE bond and the number of polymers, additionally does not take into account the rest length
# all things considered a pretty crude and mostly bad comparison... also this is COM not exact
plt.plot(Xn, Xn, label = 'Rod bounded by WEAK spring')  # Plots bounded rod's average extension
# critical comment: see the mathematica file, but increasing the spring strength sees this line drop down.
# currently an exact formulation is not in place. Fitting a spring to the FENE and treating the polymer as a
# string of series springs should allow us to fit a spring constant to our polymer parameters.
# this can be used for the previous ones too, AND we can probably adjust this to have a nonzero base extension.
# in short do some basic analytical work for the extension distance

plt.legend()

save_multi_image("Results.pdf")

plt.show()
