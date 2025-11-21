import matplotlib.pyplot as plt
from formatting import format_cursor
import design0 as d0

def init_plots(loads, span, geometry, level, b):
    fos_types = ["flexural tension", "flexural compression", "shear", "cement shear",
                "case 1", "case 2", "case 3", "case 4"]

    env = {}
    for fos in fos_types:
        all_pts = []
        for i in range(span+1):
            all_pts.append([i,0])
        env[fos] = all_pts #[[i,0],[i,0]]
    
    # env = []
    # for i in range(span+1):
    #     env.append([i, 0])

    load_mag = []
    for load_pair in loads:
        load_mag.append([load_pair[0], (load_pair[1]-loads[0][1])])

    y_bar = d0.calculate_centroidal_axis(geometry)
    I = d0.second_moment_of_area(geometry, y_bar)
    I_list = [[I, geometry, (0,span), level]]

    for i in range(span-856):
        all_fos = d0.simulation_safety_factors(load_mag, span, I_list, b)
        # print(all_fos)
        # for fos in fos_types:
            # fos_value = all_fos[fos]??
            # env[fos].append([i, fos_value])

        # for j in range(len(env)):
        #     env[j][1] = all_fos
        for pair in load_mag:
            pair[1] += 1
    return env

def plot_env(data):
    plt.axhline(0,color='black')

    for fos_type in data:
        x = []
        y = []
        for pt in data[fos_type]:
            x.append(pt[0])
            y.append(pt[1])
            print(x,y)
        plt.plot(x, y, label=fos_type)

    for pt in data:
        for y_val in pt[1]:
            plt.plot(pt[0], y_val)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(10, 73.73), 6.27, 1.27], "A3": [(83.73, 73.73), 6.27, 1.27], "A4": [(10, 1.27), 1.27, 72.46], "A5": [(88.73, 1.27), 1.27, 72.46], "A6": [(0, 75), 100, 1.27]}
    loads = [(400/6, 172), (400/6, 348), (400/6, 512), (400/6, 688), (400/6, 852), (400/6, 1028)]
    # loads = [[67.5, 172], [67.5, 348], [67.5, 512], [67.5, 688], [91.0, 852], [91.0, 1028]]
    span = 1200
    level = 1
    b = 2.54
    env = init_plots(loads, span, geometry, level, b)
    plot_env(env)