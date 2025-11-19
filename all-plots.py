import matplotlib.pyplot as plt
import design0 as d0

def init_plots(loads, span, geometry):
    reaction_forces = d0.reaction_forces(loads, span)
    y_bar = d0.calculate_centroidal_axis(geometry)
    I = d0.second_moment_of_area(geometry, y_bar)
    I_list = [[I, geometry, (0,1200)]]
    sfe = d0.calculate_shear_force(loads, reaction_forces, span)
    bme = d0.calculate_BMD(sfe)
    flex_comp = d0.flexural_stress_diagram(bme, I_list)[0]
    flex_tens = d0.flexural_stress_diagram(bme, I_list)[1]

    return sfe, bme, flex_comp, flex_tens

def plot_env(data, type=None):
    x = [0]
    y = [0]
    for pt in data:
        x.append(pt[0])
        y.append(pt[1])
    x.append(1200)
    y.append(0)
    if type == "bme":
        plt.gca().invert_yaxis()
    plt.axhline(0)
    plt.plot(x, y)
    plt.show()

# def plot_flex_stress(BMD,I):

if __name__ == "__main__":
    geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(10, 73.73), 6.27, 1.27], "A3": [(83.73, 73.73), 6.27, 1.27], "A4": [(10, 1.27), 1.27, 72.46], "A5": [(88.73, 1.27), 1.27, 72.46], "A6": [(0, 75), 100, 1.27]}
    loads = [(67.5, 172), (67.5, 348), (67.5, 512), (67.5, 688), (91.0, 852), (91.0, 1028)]
    span = 1200
    sfe, bme, flex_comp, flex_tens = init_plots(loads,span, geometry)
    plot_env(sfe)
    plot_env(bme, type="bme")
    plot_env(flex_comp)
    plot_env(flex_tens)