import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def plot_info(full_data, x_axis, y_axis, title):
    kys = full_data.keys()
    colors = cm.rainbow(np.linspace(0,1, len(kys)))
    clr = 0

    for i in kys:
        plt.plot(range(len(full_data[i][y_axis])), full_data[i][y_axis], label=i,
                    color=colors[clr], marker=None, linestyle='-')
        clr = clr + 1
   
    plt.title(title, fontsize= 22)
    plt.xlabel(x_axis, fontsize= 22)
    plt.ylabel(y_axis, fontsize= 22)
    plt.legend(fontsize=15)
    #plt.show()111
    
    plt.savefig(title + '.png')
    plt.close()
    
def main():
    df = pd.read_csv("data.csv");
    print(df)
    column_names = list(df)
    methods = pd.unique(df[column_names[0]])
    noise_levels = pd.unique(df[column_names[1]])
    print("column names")
    print(column_names)
    print("methods")
    print(methods)
    print("noise levels")
    print(noise_levels)
    statistical_info = {}
    statistical_error_info = {}
    for i in methods:
        new_info = pd.DataFrame(columns=column_names[1:(len(column_names))])
        new_info_std = pd.DataFrame(columns=column_names[1:(len(column_names))])
        for j in noise_levels:
            stat = df[df[column_names[0]] == i][(df[column_names[1]] > j - 0.1) & (df[column_names[1]] < j + 0.1)].mean()
            std = df[df[column_names[0]] == i][(df[column_names[1]] > j - 0.1) & (df[column_names[1]] < j + 0.1)].std()
            count_info = df[df[column_names[0]] == i][(df[column_names[1]] > j - 0.1) & (df[column_names[1]] < j + 0.1)].shape
            cols = list(std)
            print('COLS')
            print(cols)
            new_info = new_info.append(stat, ignore_index=True)
            new_info_std = new_info_std.append(std, ignore_index=True)
        print(new_info)
        statistical_info[i] = new_info
        statistical_error_info[i] = new_info_std
    print('ERROR INFO')
    kys = statistical_info.keys()
    for i in kys:
        print("**************")
        print(i)
        print(statistical_info[i])
    print('std INFO')
    kys = statistical_error_info.keys()
    for i in kys:
        print("**************")
        print(i)
        print(statistical_error_info[i])
        
    plot_info(statistical_info, column_names[1], column_names[2], 'error_rot')
    plot_info(statistical_info, column_names[1], column_names[3], 'trans_rot')
    plot_info(statistical_info, column_names[1], column_names[4], 'time')
    plot_info(statistical_error_info, column_names[1], column_names[2], 'error_rot_std')
    plot_info(statistical_error_info, column_names[1], column_names[3], 'trans_rot_std')
    plot_info(statistical_error_info, column_names[1], column_names[4], 'time_std')
if __name__ == "__main__":
    main()

