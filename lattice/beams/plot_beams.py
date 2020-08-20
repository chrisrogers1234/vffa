import math
import numpy
import matplotlib
import matplotlib.pyplot
import matplotlib.mlab

def gaussian_dist(x, mean, rms, norm):
    norm = norm/(2.*math.pi)**0.5 
    exponent = -0.5*(x-mean)**2/rms**2
    result = norm*math.exp(exponent)
    return result

def read_file(file_name, mean_ke):
    mass = 938.27208816
    fin = open(file_name)
    keys = ["n", "x", "x’", "y", "y’", "dPhi", "dW"]
    my_data = {}
    for key in keys:
        my_data[key] = []
    for line in fin.readlines():
        var = [float(word) for word in line.split()]
        for i, key in enumerate(keys):
            my_data[key].append(var[i])
    my_data["p"] = [((dW+mean_ke+mass)**2-mass**2)**0.5 for dW in my_data["dW"]]
    mean_p = numpy.mean(my_data["p"])
    my_data["dp_over_p"] = [(p-mean_p)/mean_p for p in my_data["p"]]
    return my_data

def plot(data, title):
    n_bins = 50
    bin_width = 0.02
    figure = matplotlib.pyplot.figure(figsize=(10, 10))
    axes = figure.add_subplot(2, 1, 1)
    bins = [-0.01 +bin_width*i/n_bins for i in range(n_bins+1)]
    mean = numpy.mean(data["dp_over_p"])
    rms = numpy.std(data["dp_over_p"])
    n_tot = len(data["dp_over_p"])
    print("Plotting", n_tot, "events with mean", mean, "rms", rms)
    gauss = [gaussian_dist(x, mean, rms, n_tot/n_bins*bin_width/rms) for x in bins]
    axes.hist(data["dp_over_p"], bins)
    axes.plot(bins, gauss)
    ylim = axes.get_ylim()
    axes.plot([mean-3*rms, mean-3*rms], [-1, n_tot], "g--")
    axes.plot([mean+3*rms, mean+3*rms], [-1, n_tot], "g--")
    axes.set_ylim(ylim)
    axes.set_xlabel("dp/p")
    axes = figure.add_subplot(2, 1, 2)
    bins = [-0.06 +0.12*i/n_bins for i in range(n_bins+1)]
    axes.hist(data["dW"], bins)
    axes.set_xlabel("dW [MeV]")
    figure.suptitle(title)
    figure.savefig(title+"_plot.png")

def main():
    for file_name in ["FETS-v7d.coord", "FETS-vFFA-BT.coord"]:
        data = read_file(file_name, 3)
        plot(data, file_name)
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")


if __name__ == "__main__":
    main()
