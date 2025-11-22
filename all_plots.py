import matplotlib.pyplot as plt
from formatting import format_cursor
import design0 as d0

def init_plots(loads, span, geometry, level):
    env_max = []
    env_min = []
    for i in range(5):
        env_max_pts = []
        env_min_pts = []
        for j in range(span+1):
            env_max_pts.append([j, 0])
            env_min_pts.append([j, float('inf')])
        env_max.append(env_max_pts)
        env_min.append(env_min_pts)

    load_mag = []
    for load_pair in loads:
        load_mag.append([load_pair[0], (load_pair[1]-loads[0][1])])
    
    for i in range(span-856):
        [sf, bm, comp, tens, glue] = d0.calculate_envs(load_mag, span, geometry, level)
        env_list = [sf, bm, comp, tens, glue]
        # print(comp)
        max_min_list = []
        for k in range(len(env_list)):
            for pt in env_list[k]:
                if abs(pt[1]) > abs(env_max[k][pt[0]][1]):
                    env_max[k][pt[0]][1] = abs(pt[1])
                if abs(pt[1]) < abs(env_min[k][pt[0]][1]):
                    env_min[k][pt[0]][1] = abs(pt[1])
            max_min_list.append((env_max[k], env_min[k]))

        for pair in load_mag:
            pair[1] += 1
    sfe, bme, comp_e, tens_e, glue_e = max_min_list

    return sfe, bme, comp_e, tens_e, glue_e

def plot_env(data, type):
    # print(data)
    labels = {
        "sfe": ["Shear Force Envelope", "Span (mm)", "Shear Force (N)"],
        "bme": ["Bending Moment Envelope", "Span (mm)", "Bending Moment (Nmm)"],
        "comp_e": ["Flexural Compression Stress Envelope", "Span (mm)", "Flexural Compression Stress (MPa)"],
        "tens_e": ["Flexural Tension Stress Envelope", "Span (mm)", "Flexural Tension Stress (MPa)"],
        "glue_e":["Shear Glue Stress Envelope", "Span (mm)", "Shear Glue Stress (MPa)"]
    }
    title, xlabel, ylabel = labels.get(type, ("", "", ""))

    legends = {
        "sfe": ["Max Shear Force", "Min Shear Force"],
        "bme": ["Max Bending Moment", "Min Bending Moment"],
        "comp_e": ["Max Flexural Compression Stress", "Min Flexural Compression Stress"],
        "tens_e": ["Max Flexural Tension Stress", "Min Flexural Tension Stress"]
    }
    max_label, min_label = legends.get(type, ("max", "min"))

    plt.axhline(0,color='black')
    points = []
    max_x = None
    max_y = None
    all_data = []
    for i, set in enumerate(data):
        if i == 0:
            label = max_label
        else:
            label = min_label
   
        x = []
        y = []
        data_set = []
        for pt in set:
            x.append(pt[0])
            if type == "comp":
                y.append(-pt[1])
            else:
                y.append(pt[1])
        data_set.append([x,y])
        all_data.append(data_set)
        plt.plot([0, x[0]], [0, y[0]], color='steelblue')
        plt.plot([x[-1], span], [y[-1], 0], color='steelblue')
        plt.plot(x, y, label=label)
        
        point = plt.plot(x, y, 'o', alpha=0)[0]
        points.append(point)

        set_max_y = max(y, key=abs)
        set_max_x = x[y.index(set_max_y)]
        max_idxs = []
        for i in range(len(y)):
            if abs(y[i]) == abs(set_max_y):
                max_idxs.append(i)

        last_idx = max_idxs[-1]
        set_max_x = x[last_idx]
        set_max_y = y[last_idx]

        if max_y is None or abs(set_max_y) >= abs(max_y):
            max_y = set_max_y
            max_x = set_max_x

    if max_x is not None and max_y is not None:
        plt.text(max_x, max_y,
                 f"({max_x:.0f}, {max_y:.4f})",
                 fontsize=10, ha='left', va='bottom', color='red')
    
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if type == "bme":
        plt.gca().invert_yaxis()
    plt.grid(True, color='lightgrey', linewidth=0.7)

    format_cursor(points)

    plt.show()
    return all_data

if __name__ == "__main__":
    geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(10, 73.73), 6.27, 1.27], "A3": [(83.73, 73.73), 6.27, 1.27], "A4": [(10, 1.27), 1.27, 72.46], "A5": [(88.73, 1.27), 1.27, 72.46], "A6": [(0, 75), 100, 1.27]}
    loads = [(400/6, 172), (400/6, 348), (400/6, 512), (400/6, 688), (400/6, 852), (400/6, 1028)]
    # loads = [[67.5, 172], [67.5, 348], [67.5, 512], [67.5, 688], [91.0, 852], [91.0, 1028]]
    freight = 500/(1.38+1.1+1)
    kN = [[freight*1.38,0], [freight*1.38,176], [freight*1, 340], [freight*1, 516], [freight*1.1, 680], [freight*1.1,856]]    
    # print(kN)
    span = 1260
    I_list = [[418480.7, geometry, (0, 1260), 1]]
    level = 1
    
    # sfe_max, sfe_min, bme, flex_comp, flex_tens = init_plots(loads,span, geometry)
    sfe, bme, comp_e, tens_e, glue_e = init_plots(loads, span, geometry, level)

    plot_env(sfe, "sfe")
    # plot_env(bme, "bme")
    # plot_env(comp_e, "comp_e")
    # plot_env(tens_e, "tens_e")
    # plot_env(glue_e, "glue_e")