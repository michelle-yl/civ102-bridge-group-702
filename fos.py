import matplotlib.pyplot as plt
from formatting import format_cursor
import design0 as d0

def plot_env(span, I_list, b):
    labels = [
        "FOS Shear", "FOS Flexural Tension", "FOS Flexural Compression", 
        "FOS Shear Glue Tabs", "FOS  Shear Sheets", 
        "FOS Plate Buckling Case 1", "FOS Plate Buckling Case 2",
        "FOS Plate Buckling Case 3", "FOS Plate Buckling Case 4"
    ]

    fos_list = d0.simulation_safety_factors_across_bridge(span, I_list, b)
    # print(fos_list[8])
    
    env = []
    for i in range(len(fos_list)):
        env.append([])
        for point in fos_list[i]:
            env[i].append([point[0], point[1]])

    for j in range(len(env)):
        fos_type = env[j]
        x = []
        y = []
        for pt in fos_type:
            x.append(pt[0])
            y.append(abs(pt[1]))
        if x and y:
            plt.plot(x, y, label=labels[j])
    plt.legend()
    plt.show()

if __name__ == "__main__":
    geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(88.73, 1.27), 1.27, 72.46], "A3": [(10, 1.27), 1.27, 72.46], "A4": [(83.73, 73.73), 6.27, 1.27], "A5": [(10, 73.73), 6.27, 1.27], "A6": [(0, 75), 100, 1.27]}
    loads = [[400/6, 172], [400/6, 348], [400/6, 512], [400/6, 688], [400/6, 852], [400/6, 1028]]    # loads = [[67.5, 172], [67.5, 348], [67.5, 512], [67.5, 688], [91.0, 852], [91.0, 1028]]
    span = 1200
    I_list = [[418480.7, geometry, (0, 1200), 1]]
    level = 1
    b = 2.54
    plot_env(span, I_list, b)
    