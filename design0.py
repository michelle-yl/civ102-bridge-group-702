
# geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}


# calculating centroidal axis
def areas(geometry):
    print(geometry.items())
    areas_dict = {}
    for key, value in geometry.items():
        width = value[1]
        height = value[2]
        area = width * height
        areas_dict[key] = area
    return areas_dict # areas = {A1: area1, A2: area2, ...}

def find_centroids(geometry):
    centroids = []
    for value in geometry.values():
        anchor = value[0]
        height = value[2]
        c_y = anchor[1] + height / 2
        centroids.append((c_y))
    return centroids

def calculate_centroidal_axis(geometry):
    areas_dict = areas(geometry)
    total_area = sum(areas_dict.values())
    centroids = find_centroids(geometry)
    list_areas = []
    for area in areas_dict.values():
        list_areas.append(area)
    y_bar = 0
    for i in range(len(centroids)):
        y_bar += list_areas[i] * centroids[i]
    y_bar /= total_area
    return y_bar

# calculating I
def inertia(component):
    return component[1] * component[2]**3 / 12


def second_moment_of_area(geometry, y_bar):
    I_total = 0
    list_areas = areas(geometry)
    components = list(geometry.values())
    centroids = find_centroids(geometry)
    for i in range(len(components)):
        I = inertia(components[i])
        d_sq = (centroids[i] - y_bar)**2
        I_total += I + list(list_areas.values())[i] * d_sq
    return I_total

# reaction forces

# loads = [(load1, position1), (load2, position2), ...]

def reaction_forces(loads, span):
    A_y = 0
    B_y = 0
    loc_A_y = span / 2 - 600
    loc_B_y = span / 2 + 600
    total_load = 0
    for load in loads:
        total_load += load[0]
    # find B_y taking moment about A
    M_A = 0
    for load in loads:
            M_A -= load[0]*(loc_A_y - load[1])
    print(M_A)
    B_y = M_A / 1200
    A_y = total_load - B_y
    return [(A_y, loc_A_y), (B_y, loc_B_y)]

def update_loads(loads, span):
    for load in loads:
        load[1] += 1
    # assume loads move left to right, leftmost load begins at x = 0


# Params:
# (list) loads: list of tuples (magnitude of load, location)
# (list) reaction_forces: list of tuples (magnitude of force, location)
# (int)  span: length of bridge

def calculate_shear_force(loads, reaction_forces, span):
    shear_force_diagram = []
    V = 0
    for x in range(span):
        for force in reaction_forces:
            if force[1] == x:
                V += force[0]
        for load in loads:
            if load[1] == x:
                V -= load[0]
        shear_force_diagram.append([x, V])
    return shear_force_diagram



# max shear stress at a location x
# I = [[I, geometry, (start, end)], ...]
def max_shear_stress(shear_force_diagram, I, b):
    SFD = shear_force_diagram[:]
    for value in SFD:
        value[1] = abs(value[1])
    
    V_maxs = []
    for i in I:
        slice_of_SFD = SFD[i[2][0]:i[2][1]]
        max_V = max(SFD, key = lambda x: x[1])
        V_maxs.append(max_V)   

    shear_stresses = []
    for i in range(len(V_maxs)):
        Q = calculate_Qmax(I[i][1])
        tau = V_maxs[i] * Q / (I[i][0] * b)
        shear_stresses.append([tau, I[i][0]])
    return shear_stresses

def calculate_Qmax(geometry):
    geo_below_ybar = geometry
    y_bar = calculate_centroidal_axis(geometry)
    for component in geo_below_ybar.values():
        if component[2] > y_bar:
            component[2] = y_bar
    y_bar_below = calculate_centroidal_axis(geo_below_ybar)
    area_below = areas(geo_below_ybar)
    return area_below * (y_bar - y_bar_below)
        

# find BMD
def calculate_BMD(SFD):
    bending_moment_diagram = []
    M = 0
    for x in SFD:
        M += x[1]
        BMD.append([x[0], M])
    return bending_moment_diagram

# safety factor

def safety_factor(applied_stress, type):
    # all stresses in MPa
    allowable_stresses = {tensile: 30, compressive: 6, shear: 4, cement_shear: 1.5} #cement_shear is actually 2, but that's only if properly cured
    return allowable_stresses[type] / applied_stress

def initialize_loads():
    return [(67.5, 0), (67.5, 176), (67.5, 340), (67.5, 516), (91.0, 680), (91.0, 856)]

if __name__ == "__main__":
    #geometry = {"A1": [(0, 75), 100, 1.27], "A2": [(10, 73.73), 6.27, 1.27], "A3": [(83.73, 73.73), 6.27, 1.27], "A4": [(10, 1.27), 1.27, 72.46], "A5": [(88.73, 1.27), 1.27, 72.46], "A6": [(10, 0), 80, 1.27]}
    #centroidal_axis = calculate_centroidal_axis(geometry)
    #print(f"Centroidal Axis (ȳ): {centroidal_axis:.4f} mm")
    #I = second_moment_of_area(geometry, centroidal_axis)
    #print(f"Second Moment of Area (I): {I:.4f} mm^4")
    #loads = [(67.5, 172), (67.5, 348), (67.5, 512), (67.5, 688), (91.0, 852), (91.0, 1028)]
    #span = 1200
    loads = [(50, 25), (100, 1275)]
    span = 1300
    A_y, B_y = reaction_forces(loads, span)
    print(f"Reaction Forces: A_y = {A_y:.2f} N, B_y = {B_y:.2f} N")
