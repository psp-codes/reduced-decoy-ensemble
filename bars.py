# Developed by: Shehu Lab
# Generates bar charts for the results of the lRMSD distances

import matplotlib.pyplot as plt
import configparser as cp
import matplotlib as mpl


# helper function
def FloatsFromString(string_with_t_sep_floats):
    x = []
    for y in string_with_t_sep_floats.split('\t'):
            x.append(float(y))
    return x


############################################################

# read configuration file
conf = cp.ConfigParser()
conf.read("configs/bar_charts.ini")

# paths
inputDir = conf['input']['inputDir']
outputDir = conf['output']['outputDir']

# setup for figure
mpl.rcParams['mathtext.default'] = 'regular'

pdf_figName="min_RMSD_comparison.pdf"
png_figName="min_RMSD_comparison.png"

posMetricOfInterest = 0

nameMetricOfInterest = 'Minimum ' + 'CA' + ' lRMSD to Native Structure (' + r'$\AA$' + ')'

plt.figure(figsize=(4,2.5))
plt.grid(True, alpha=0.5)
plt.tick_params(axis='both', which='major', labelsize=4)
ax=plt.gca()
plt.ylabel(nameMetricOfInterest,fontsize=6)
plt.xlabel("Target", fontsize=6)

colors=["#d53e4f", "#1a9850", "#542788", "#a6611a", "#3288bd", "#ff7f00"]

bar_width = 0.11

systems_array = [ "T0859-D1", "T0898-D2", "T0886-D1", "T0897-D1", "T0892-D2",
                   "T0957s1-D1", "T0953s1-D1", "T0953s2-D3", "T0960-D2", "T1008-D1"]

plt.xticks(range(len(systems_array)), systems_array, rotation=45)

file_dict = {}  # create dictionary

x_loc = 0
for test_case in systems_array:
        inputFilename = inputDir + str(test_case) + ".dat"

        # read data as list of dictionaries
        with open(inputFilename) as f:
            file_dict[test_case] = {}
            for line in f:
                key, value = line.split(": ")
                file_dict[test_case][key] = value.rstrip()

        # extract metrics for each algorithm, for this current test_case
        orig_metrics = FloatsFromString(file_dict[test_case]["Original"])
        hier_metrics = FloatsFromString(file_dict[test_case]["Hierarchical"])
        km_metrics = FloatsFromString(file_dict[test_case]["K-means"])
        gmm_metrics = FloatsFromString(file_dict[test_case]["GMM"])
        gromos_metrics = FloatsFromString(file_dict[test_case]["gromos"])
        trunc_metrics = FloatsFromString(file_dict[test_case]["Truncation"])
        
        # draw the bars for this current test_case
        ax = plt.subplot(111)
        orig_rects = ax.bar(x_loc,   orig_metrics[posMetricOfInterest],  width=bar_width, color=colors[0], align='center')
        hier_rects = ax.bar(x_loc + 1 * bar_width, hier_metrics[posMetricOfInterest], width=bar_width, color=colors[1], align='center')
        km_rects = ax.bar(x_loc + 2 * bar_width, km_metrics[posMetricOfInterest], width=bar_width, color=colors[2], align='center')
        gmm_rects = ax.bar(x_loc + 3 * bar_width, gmm_metrics[posMetricOfInterest],  width=bar_width, color=colors[3], align='center')
        gromos_rects = ax.bar(x_loc + 4 * bar_width, gromos_metrics[posMetricOfInterest], width=bar_width, color=colors[4], align='center')
        trunc_rects = ax.bar(x_loc + 5 * bar_width, trunc_metrics[posMetricOfInterest], width=bar_width, color=colors[5], align='center')

        x_loc +=1
        
# finish up figure business
ax.autoscale(tight=True)

clustering_methods = ['Original', 'Hierarchical', 'K-means', 'GMM', 'Gmx-cluster', 'Truncation']

ax.legend(clustering_methods, bbox_to_anchor=(1.01, 1.015), fontsize='xx-small')
ax.set_ybound(lower=0, upper=15)

plt.savefig(outputDir + pdf_figName, dpi=300, transparent='True', bbox_inches='tight', pad_inches=0.05)
plt.savefig(outputDir + png_figName, dpi=300, transparent='True', bbox_inches='tight', pad_inches=0.05)
 
