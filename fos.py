import matplotlib.pyplot as plt
from formatting import format_cursor
import design0 as d0

def plot_env(span, I_list, b):
    fos_list = d0.simulation_safety_factors_across_bridge(span, I_list, b)
    print(fos_list)
    
    env = []
    for k in range(len(fos_list)):
        env.append([])
        for point in fos_list[k]:
            env[k].append([point[0], point[1]])
    # print(env)
    
    # env = []
    # for i in range(9):
    #     env.append([])

    # for point in fos_list[k]:
    #     env[k].append([point[0], point[1]])

    for k in range(len(env)):
        fos_type = env[k]
        x = []
        y = []
        for pt in fos_type:
            x.append(pt[0])
            y.append(pt[1])
        if x and y:
            plt.plot(x, y)

    # for fos_type in env:
    #     x = []
    #     y = []
    #     for pt in fos_type:
    #         x.append(pt[0])
    #         y.append(pt[1])
    #     plt.plot(x,y)
    plt.show()

if __name__ == "__main__":
    geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(10, 73.73), 6.27, 1.27], "A3": [(83.73, 73.73), 6.27, 1.27], "A4": [(10, 1.27), 1.27, 72.46], "A5": [(88.73, 1.27), 1.27, 72.46], "A6": [(0, 75), 100, 1.27]}
    loads = [[400/6, 172], [400/6, 348], [400/6, 512], [400/6, 688], [400/6, 852], [400/6, 1028]]    # loads = [[67.5, 172], [67.5, 348], [67.5, 512], [67.5, 688], [91.0, 852], [91.0, 1028]]
    span = 1200
    I_list = [[418480.7, geometry, (0, 1200), 1]]
    level = 1
    b = 2.54
    plot_env(span, I_list, b)
    