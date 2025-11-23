# fos.py generates factor of safety envelopes for shear force, flexural tension, flexural compression, shear in glue tabs, shear in glue sheets, and all 4 cases of plate buckling

# Importing "matplotlib" to generate graphs, "format_cursor" for graph hover tool, "BDP_calcs" file for calculation functions
import matplotlib.pyplot as plt
from formatting import format_cursor
import BDP_calcs as bdp

# Generates and plots all envelopes 
def plot_env(span, geometry, b, layers):
    # Parameters:
        # span (int): length of bridge (mm)
        # geometry (dict): dictionary of components with their location and dimensions
        #   geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ..., An: [anchor, width, height]}
        #       An is the nth component of the cross-section.  Components are listed from bottom to top, left to right of cross-section.
        #       Anchor = (x, y) --> coordinate of bottom-left corner of component.  y = 0 occurs at the bottom of the cross-section.
        #       width = width of component
        #       height = height of component
        # b (int): width of the beam's cross section at the centroidal axis (mm)
        # layers (int): the number of matboard layers making the top flange (1, 2, 3, ...)
    # Returns:
        # None (uses matplotlib to generate plot)

    #  Creates list of lists of cross-sectional features (cross-sections change for different sections of the span)
    #   I = [[I1, geometry1, (start, end)], [I2, geometry2, (start, end)]...]
    I_list = [[bdp.second_moment_of_area(geometry, bdp.calculate_centroidal_axis(geometry)), geometry, (0, span), layers]]
    labels = [
        "FOS Shear", "FOS Flexural Tension", "FOS Flexural Compression", 
        "FOS Shear Glue Tabs", "FOS Shear Sheets", 
        "FOS Plate Buckling Case 1", "FOS Plate Buckling Case 2",
        "FOS Plate Buckling Case 3", "FOS Plate Buckling Case 4"
    ]

    # Generates 9 lists of x-coordinate and FOS pairs using BDP_calcs.py function
    #   FOS_'type'_diagram = [[x1, FOS1], [x2, FOS2], ...] ---> 9 of these lists in a tuple
    fos_list = bdp.simulation_safety_factors_across_bridge(span, I_list, b)
    
    # Formats fos_list into a list of 9 lists, instead of tuple
    env = []
    points = []
    for i in range(len(fos_list)):
        env.append([])
        for point in fos_list[i]:
            env[i].append([point[0], point[1]])

    # Graphs 9 plots for all fos types
    for j in range(len(env)):
        fos_type = env[j]
        x = []
        y = []
        for pt in fos_type:
            x.append(pt[0])
            y.append(abs(pt[1]))
        if x and y:
            line_list = plt.plot(x, y, label=labels[j])
            line = line_list[0]
            color = line.get_color()  

            # tracks cursor with hidden plotted points
            point = plt.plot(x, y, 'o', alpha=0, color=color)[0]
            points.append(point)

    plt.title("Log Factory of Safety vs. Train Position")
    plt.xlabel("Position of Train's Left End (mm)")
    plt.ylabel("Log Factor of Safety")
    plt.yscale("log")
    plt.legend(loc="center left", bbox_to_anchor=(1.02, 1))
    plt.tight_layout(rect=[0, 0, 1.17, 1])
    plt.grid(True, color='lightgrey', linewidth=0.7)

    format_cursor(points)
    plt.show()

if __name__ == "__main__":
    span = 1260
    b = 2.54

    # # DESIGN0
    # geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(10, 1.27), 1.27, 72.46], "A3": [(88.73, 1.27), 1.27, 72.46], "A4": [(10, 73.73), 6.27, 1.27], "A5": [(83.73, 73.73), 6.27, 1.27], "A6": [(0, 75), 100, 1.27]}
    
    # IT14
    geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(10, 1.27), 6.27, 1.27], "A3": [(83.73, 1.27), 6.27, 1.27], "A4": [(10, 1.27), 1.27, 134.92], "A5": [(88.73, 1.27), 1.27, 134.92], "A6": [(10, 136.19), 6.27, 1.27], "A7": [(83.73, 136.19), 6.27, 1.27], "A8": [(0, 137.46), 100, 1.27], "A9": [(0, 138.73), 100, 1.27]}
    layers = 2

    # # LOAD 400N
    # loads = [[400/6, 172], [400/6, 348], [400/6, 512], [400/6, 688], [400/6, 852], [400/6, 1028]]    

    # LOAD 452N
    loads = [[67.5, 0], [67.5, 176], [67.5, 340], [67.5, 516], [91.0, 680], [91.0, 856]]

    # # LOAD 1020N
    # freight = 510/(1.38+1.1+1)
    # loads = [[freight*1.38, 0], [freight*1.38, 176], [freight*1, 340], [freight*1, 516], [freight*1.1, 680], [freight*1.1, 856]]    
    
    plot_env(span, geometry, b, layers)
    