from matplotlib import pyplot as plt



# Histogram Analysis
def histogram(pd_table, xlabel, ylabel, title, legend_title, cell_generations, bins):
    # Extract data
    x_values = []
    y_values = []
    for gen in cell_generations:
        cell_gen_data = pd_table[pd_table['Cell Generation'] == gen]

        cell_speeds = []
        for s in list(cell_gen_data['Cell Speed (um/s)']):
            cell_speeds.append(s)
        [array, bins, patches] = plt.hist(cell_speeds, weights=np.ones(len(cell_speeds)) / len(cell_speeds), bins=bins, label=str(gen), histtype='step')
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        x_values.append(bins)
        y_values.append(array)

    # Figure details
    plt.figure(figsize=(8, 6), dpi=80)
    plt.clf()
    for i in range(len(x_values)):
        plt.plot(x_values[i][1:], y_values[i], label=str(cell_generations[i]))
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.xticks(fontsize=12, ticks=[0, 2.5E-4, 5E-4, 7.5E-4, 10E-4, 12.5E-4, 15E-4, 17.5E-4], labels=["0", "2.5E-4", "5E-4", "7.5E-4", "10E-4", "12.5E-4", "15E-4", "17.5E-4"])
    plt.yticks(fontsize=12)
    plt.gcf().subplots_adjust(left=0.2)
    plt.gcf().subplots_adjust(bottom=0.2)

    plt.title(title)
    plt.grid(True)
    plt.legend(title=legend_title)
    plt.show()

#  Persistence Time
def persistence_time(pd_table, ylabel, title, tf, ts):
    cells_accounted = {}
    for t in np.arange(0, tf, ts).tolist():
        time_data = pd_table[pd_table['Time (hr)'] == t*ts]

        for c in list(time_data['Cell ID']):
            if c in cells_accounted:
                row = time_data[time_data['Cell ID'] == c]
                row_list = row.values.tolist()[0]
                if cells_accounted[c][2] == row_list[3]:
                    cells_accounted[c][1] += 1
                    if cells_accounted[c][1] > cells_accounted[c][0]:
                        cells_accounted[c][0] = cells_accounted[c][1]
                else:
                    cells_accounted[c][1] = 0
            else:
                row = time_data[time_data['Cell ID'] == c]
                row_list = row.values.tolist()[0]
                cells_accounted[c] = [0, 0, row_list[3]]
    
    persistence_times = []
    for c in cells_accounted.keys():
        persistence_times.append(cells_accounted[c][0])

    plt.boxplot(msd)
    plt.ylabel(ylabel, fontsize=16)
    plt.title(title, fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.gcf().subplots_adjust(left=0.18)
    plt.gcf().subplots_adjust(bottom=0.18)

    plt.show()