import matplotlib.pyplot as plt
import design0 as d0

def plot_sfe(loads,span):
    reaction_forces = d0.reaction_forces(loads, span)
    sfe = d0.calculate_shear_force(loads, reaction_forces, span)
    x = [0]
    y = [0]
    for pt in sfe:
        x.append(pt[0])
        y.append(pt[1])
    x.append(1200)
    y.append(0)
    plt.axhline(0)
    plt.plot(x, y)
    plt.show()

if __name__ == "__main__":
    loads = [(67.5, 172), (67.5, 348), (67.5, 512), (67.5, 688), (91.0, 852), (91.0, 1028)]
    span = 1200
    plot_sfe(loads,span)