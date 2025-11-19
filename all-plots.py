import matplotlib.pyplot as plt
import design0 as d0

def init_plots(loads, span):
    reaction_forces = d0.reaction_forces(loads, span)
    sfe = d0.calculate_shear_force(loads, reaction_forces, span)
    bme = d0.calculate_BMD(sfe)

    return sfe,bme

def plot_env(data):
    x = [0]
    y = [0]
    for pt in data:
        x.append(pt[0])
        y.append(pt[1])
    x.append(1200)
    y.append(0)
    plt.axhline(0)
    plt.plot(x, y)
    plt.show()

# def plot_flex_stress(BMD,I):

if __name__ == "__main__":
    loads = [(67.5, 172), (67.5, 348), (67.5, 512), (67.5, 688), (91.0, 852), (91.0, 1028)]
    span = 1200
    sfe, bme = init_plots(loads,span)
    plot_env(sfe)