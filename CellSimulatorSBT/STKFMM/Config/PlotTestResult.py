import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

import argparse
import json


def parseError(error):
    '''
    convert error from a list of dict to numpy 2D array
    '''
    if not error:
        return None

    result = []
    for item in error:
        result.append([item["drift"], item["driftL2"],
                       item["errorRMS"], item["errorL2"], item["errorMaxRel"]])
    return np.transpose(np.array(result))


def plotRecord(ax, multOrder, treeTime, runTime, name, error):
    dim = error.shape[2]
    for i in range(dim):
        ax.semilogy(multOrder, error[:, 1, i], 'x', label=str(i))
    ax.set_prop_cycle(None)  # reset color cycle
    for i in range(dim):
        ax.semilogy(multOrder, error[:, 3, i], '--o')

    ax.set_title(name)
    ax.legend(loc='upper left', ncol=3,
              title=r'component error: x $\delta_{L2}$, o $\epsilon_{L2}$')
    ax.set_xlabel(r"$p$")
    ax.set_ylabel("Error")
    ax.set_ylim(1e-15, 1)
    axt = ax.twinx()
    axt.plot(multOrder, runTime, '-', label=r'$t_{run}$')
    axt.plot(multOrder, treeTime, '-', label=r'$t_{tree}$')
    axt.set_ylabel("time")
    axt.set_ylim(0, max(runTime)*1.5)
    axt.legend(title='time')

    ax.tick_params(axis="x", direction="in")
    ax.tick_params(axis="y", direction="in")
    axt.tick_params(axis="y", direction="in")

    return


def plotData(data):
    '''
    data should be records for the same kernel
    '''
    data.sort(key=lambda k: k['multOrder'])
    kernel = data[0]["kernel"]
    multOrder = []
    treeTime = []
    runTime = []
    error = dict()
    if data[0]["errorVerify"]:
        error["Verification Error"] = []
    if data[0]["errorConvergence"]:
        error["Convergence Error"] = []
    if data[0]["errorTranslate"]:
        error["Translation Error"] = []

    for record in data:
        multOrder.append(record["multOrder"])
        treeTime.append(record["treeTime"])
        runTime.append(record["runTime"])
        if record["errorVerify"]:
            error["Verification Error"].append(
                parseError(record["errorVerify"]))
        if record["errorConvergence"]:
            error["Convergence Error"].append(
                parseError(record["errorConvergence"]))
        if record["errorTranslate"]:
            error["Translation Error"].append(
                parseError(record["errorTranslate"]))

    nax = len(error.keys())
    fig = plt.figure(
        figsize=(6.0*nax, 5.0), dpi=150, constrained_layout=True)
    axs = fig.subplots(nrows=1, ncols=nax, squeeze=False)
    index = 0
    for k in error.keys():
        error[k] = np.array(error[k])
        ax = axs.flat[index]
        plotRecord(ax, multOrder, treeTime, runTime, kernel+' '+k, error[k])
        index += 1
    plt.savefig('Test_'+kernel+'.png')


parser = argparse.ArgumentParser()
parser.add_argument("logfile")
args = parser.parse_args()
print("Parsing "+args.logfile)

f = open(args.logfile,)
logs = json.load(f)

kernelset = set()
for record in logs:
    kernelset.add(record["kernel"])
print("Found logs for kernels: ", kernelset)

for kernel in kernelset:
    data = []
    for record in logs:
        if record["kernel"] == kernel:
            data.append(record)
    plotData(data)
