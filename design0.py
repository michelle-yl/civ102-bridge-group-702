
# geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}


# calculating centroidal axis
def areas(geometry):
    # print(geometry.items())
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

# loads = [[load1, position1], [load2, position2], ...]

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
    # print(M_A)
    B_y = M_A / 1200
    A_y = total_load - B_y
    return [(A_y, loc_A_y), (B_y, loc_B_y)]

def update_loads(loads, direction):
    if direction == "right":
        for load in loads:
            load[1] += 1
    elif direction == "left":
        for load in loads:
            load[1] -= 1
    return loads



# Params:
# (list) loads: list of lists [magnitude of load, location]
# (list) reaction_forces: list of tuples (magnitude of force, location)
# (int)  span: length of bridge

def calculate_shear_force(loads, reaction_forces, span):
    shear_force_diagram = []
    V = 0
    for x in range(span+1):
        for force in reaction_forces:
            if force[1] == x and x != 1200:
                V += force[0]
        for load in loads:
            if load[1] == x:
                V -= load[0]
        shear_force_diagram.append([x, V])
    return shear_force_diagram

def calculate_envs(load_mag, span, geometry):
    react = reaction_forces(load_mag, span)
    sf = calculate_shear_force(load_mag, react, span)

    bm = calculate_BMD(sf)

    y_bar = calculate_centroidal_axis(geometry)
    I = second_moment_of_area(geometry, y_bar)
    I_list = [[I, geometry, (0,span)]]
    flex_stress = flexural_stress_diagram(bm, I_list)
    flex_comp = flex_stress[0]
    flex_tens = flex_stress[1]
    
    return [sf, bm, flex_comp, flex_tens]

# max shear stress at a location x
# I = [[I1, geometry, (start, end)], ...]
def shear_stress_diagram(shear_force_diagram, I_list, b):
    Q_maxs = []
    for i in range(len(I_list)):
        Q = calculate_Qmax(I_list[i][1])
        Q_maxs.append(Q)
    shear_stresses_diagram = []
    SFD = shear_force_diagram[:]
    for i in range(len(I_list)):
        for x in range(I_list[2][0], I_list[2][1]):
            shear_stress = SFD[x][1] * Q_maxs[i] / (I_list[i][0] * b)
            shear_stresses_diagram.append([x, shear_stress])
    return shear_stresses_diagram

# geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}
def calculate_Qmax(geometry):
    geo_below_ybar = geometry.copy()
    y_bar = calculate_centroidal_axis(geometry)
    for component in geo_below_ybar.values():
        if component[2] + component[0][1] > y_bar:
            component[2] = y_bar
    y_bar_below = calculate_centroidal_axis(geo_below_ybar)
    area_below = sum(areas(geo_below_ybar).values())
    return area_below * (y_bar - y_bar_below)
        

# find BMD
def calculate_BMD(SFD):
    bending_moment_diagram = []
    M = 0
    for x in SFD:
        M += x[1]
        bending_moment_diagram.append([x[0], M])
    return bending_moment_diagram

# find flexural stress diagram
# I = [[I, geometry, (start, end), layers], ...]
# geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}
# assume geometry is organized from bottom to top of cross-section.  The last component, An, is the topmost component.
def flexural_stress_diagram(BMD, I):

    flexural_compression_diagram = []
    flexural_tension_diagram = []
    
    for i in range(len(I)):
        geometry = I[i][1]
        upper_component_dimensions = list(geometry.values())[-1]
        height = upper_component_dimensions[0][1] + upper_component_dimensions[2]
        y_bar = calculate_centroidal_axis(geometry)
        y_compression = y_bar - height
        y_tension = y_bar

        for x in range(I[i][2][0], I[i][2][1]):

            sigma_compression = BMD[x][1] * y_compression / I[i][0]
            sigma_tension = BMD[x][1] * y_tension / I[i][0]

            flexural_compression_diagram.append([x, sigma_compression])
            flexural_tension_diagram.append([x, sigma_tension])

    return flexural_compression_diagram, flexural_tension_diagram

    # for x in BMD:

# calculate plate buckling stress

# geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}
def plate_buckling_stress(geometry, case, layers, a = None):
    if layers == 2:
        t = list(geometry.values())[-1][2] + list(geometry.values())[-2][2] # thickness of flange or web
    else:
        t = list(geometry.values())[-1][2]

    # case 1: buckling of compressive flange between webs
    if case == 1:
        b = list(geometry.values())[-1][1] - list(geometry.values())[-layers][0][0] * 2 # width of flange between webs
        k = 4
        sigma = k * (3.14159**2) * 4000 * (t / b)**2 / (12 * (1 - 0.2**2))
        return sigma
    
    # case 2: buckling of the tips of the compressive flange
    if case == 2:
        b = list(geometry.values())[-(layers+1)][0][0] # length of flange tip beyond web
        k = 0.425
        sigma = k * (3.14159**2) * 4000 * (t / b)**2 / (12 * (1 - 0.2**2))
        return sigma
    
    # case 3: buckling of the webs due to the flexural stresses
    if case == 3:
        t = 1.27 # thickness of web
        h = list(geometry.values())[-layers][2] # height of web
        k = 6
        sigma = k * (3.14159**2) * 4000 * (t / h)**2 / (12 * (1 - 0.2**2))
        return sigma
    
    # case 4: buckling of the webs due to the shear stresses
    if case == 4:
        t = 1.27 # thickness of web
        h = list(geometry.values())[-layers][2] # height of web
        k = 5
        tau = k * (3.14159**2) * 4000 * ((t / h)**2 + (t / a)**2) / (12 * (1 - 0.2**2))
        return tau




# shear glue stress
# I = [[I, geometry, (start, end), layers], ...]
def shear_glue_stress_diagram(shear_force_diagram, I, level): # level is an integer indicating which glue line to analyze, counting from the top down
    shear_glue_stresses_diagram = []
    for i in range(len(I)):
        geometry = I[i][1]
        layers = I[i][3]
        if level == 1 and layers == 1:
            upper_component_dimensions = list(geometry.values())[-level]
            d = upper_component_dimensions[2]/2 # d from centroid of area above glue line to top shear-stress free surface
            A = upper_component_dimensions[1] * upper_component_dimensions[2] * d
            b = list(geometry.values())[-2][1]*2 # assume the second-to-top component is the web with glue tab, so its width*2 is the glue line length
        
        if level == 1 and layers == 2:
            upper_component_dimensions = list(geometry.values())[-level]
            second_upper_component_dimensions = list(geometry.values())[-(level-1)]
            A = upper_component_dimensions[1] * upper_component_dimensions[2]
            d = upper_component_dimensions[2]/2 # d from centroid of area above glue line to top shear-stress free surface
            b = second_upper_component_dimensions[1]

        if level == 2 and layers == 2:
            upper_component_dimensions = list(geometry.values())[-level]
            second_upper_component_dimensions = list(geometry.values())[-(level-1)]
            d = upper_component_dimensions[0][1] + upper_component_dimensions[2] - calculate_centroidal_axis({"upper": upper_component_dimensions, "second_upper": second_upper_component_dimensions}) # d from centroid of area above glue line to top shear-stress free surface
            A = (upper_component_dimensions[1] * upper_component_dimensions[2]) + (second_upper_component_dimensions[1] * second_upper_component_dimensions[2])
            b = list(geometry.values())[-3][1]*2 # assume the third-to-top component is the web with glue tab, so its width*2 is the glue line length
        
        Q = A * d

        for x in range(I[i][2][0], I[i][2][1]):
            shear_glue_stress = shear_force_diagram[x][1] * Q / (I[i][0] * b)
            shear_glue_stresses_diagram.append([x, shear_glue_stress])

    return shear_glue_stresses_diagram


# safety factor

def safety_factor(applied_stress, type):
    # all stresses in MPa
    allowable_stresses = {"tensile": 30, "compressive": 6, "shear": 4, "cement_shear": 1.5} #cement_shear is actually 2, but that's only if properly cured
    return allowable_stresses[type] / applied_stress

# I = [[I_value, geometry, (start, end), layers], ...]
# In = [I_value, geometry, (start, end), layers]
def safety_factor_thin_plate(In, flexural_compression_diagram, shear_stress_diagram, a = None):
    geometry = In[1]
    layers = In[3]

    case_1_failure = plate_buckling_stress(geometry, 1, layers)
    max_flex_comp = max(flexural_compression_diagram, key = lambda x: abs(x[1]))[1]
    FOS1 = case_1_failure / abs(max_flex_comp)

    case_2_failure = plate_buckling_stress(geometry, 2, layers)
    FOS2 = case_2_failure / abs(max_flex_comp)

    case_3_failure = plate_buckling_stress(geometry, 3, layers)
    FOS3 = case_3_failure / abs(max_flex_comp)

    case_4_failure = plate_buckling_stress(geometry, 4, layers, a)
    max_shear = max(shear_stress_diagram, key = lambda x: abs(x[1]))[1]
    FOS4 = case_4_failure / abs(max_shear) # ***double check the theory behind these calculations***

    return FOS1, FOS2, FOS3, FOS4

def initialize_loads():
    return [(67.5, 0), (67.5, 176), (67.5, 340), (67.5, 516), (91.0, 680), (91.0, 856)]

def simulation_safety_factors(loads, span, I):
    min_safety_factors = {"flexural tension": 1000, "flexural compression": 1000, "shear": 1000, "cement shear": 1000, "case 1": 1000, "case 2": 1000, "case 3": 1000, "case 4": 1000} # [flexural tension, flexural compression, shear, cement shear, plate buckling case 1, case 2, case 3, case 4]
    for x in range(span - loads[-1][1]):
        # shear stress
        shear_stress_profile = shear_stress_diagram(calculate_shear_force(loads, reaction_forces(loads, span), span), I)
        max_shear = max(shear_stress_profile, key = lambda x: abs(x[1]))[1]
        FOS_shear = safety_factor(max_shear, "shear")
        if abs(FOS_shear) < min_safety_factors["shear"]:
            min_safety_factors["shear"] = abs(FOS_shear)
        
        # flexural stress and flexural tension

        BMD = calculate_BMD(shear_stress_profile)
        flex_comp, flex_tens = flexural_stress_diagram(BMD, I)
        max_flex_tens = max(flex_tens, key = lambda x: abs(x[1]))[1]
        FOS_flex_tens = safety_factor(max_flex_tens, "tensile")

        if abs(FOS_flex_tens) < min_safety_factors["flexural tension"]:
            min_safety_factors["flexural tension"] = abs(FOS_flex_tens)
        max_flex_comp = max(flex_comp, key = lambda x: abs(x[1]))[1]
        FOS_flex_comp = safety_factor(abs(max_flex_comp), "compressive")

        if FOS_flex_comp < min_safety_factors["flexural compression"]:
            min_safety_factors["flexural compression"] = FOS_flex_comp

        # cement shear stress
        reaction_forces_list = reaction_forces(loads, span)
        shear_glue_stress_profile = shear_glue_stress_diagram(calculate_shear_force(loads, reaction_forces_list, span), I, 1)
        max_cement_shear = max(shear_glue_stress_profile, key = lambda x: abs(x[1]))[1]
        FOS_cement_shear = safety_factor(max_cement_shear, "cement_shear")
        if abs(FOS_cement_shear) < min_safety_factors["cement shear"]:
            min_safety_factors["cement shear"] = abs(FOS_cement_shear)
        
        # plate buckling
        for i in range(len(I)):
            geometry = I[i][1]
            layers = I[i][3]
            FOS1, FOS2, FOS3, FOS4 = safety_factor_thin_plate(I[i], flex_comp, shear_stress_profile, 400) # I put a random a for now
            if abs(FOS1) < min_safety_factors["case 1"]:
                min_safety_factors["case 1"] = abs(FOS1)
            if abs(FOS2) < min_safety_factors["case 2"]:
                min_safety_factors["case 2"] = abs(FOS2)
            if abs(FOS3) < min_safety_factors["case 3"]:
                min_safety_factors["case 3"] = abs(FOS3)
            if abs(FOS4) < min_safety_factors["case 4"]:
                min_safety_factors["case 4"] = abs(FOS4)
        
        loads = update_loads(loads, "right")

    for x in range(span - loads[-1][1]):
        # shear stress
        shear_stress_profile = shear_stress_diagram(calculate_shear_force(loads, reaction_forces(loads, span), span), I)
        max_shear = max(shear_stress_profile, key = lambda x: abs(x[1]))[1]
        FOS_shear = safety_factor(max_shear, "shear")
        if abs(FOS_shear) < min_safety_factors["shear"]:
            min_safety_factors["shear"] = abs(FOS_shear)
        
        # flexural stress and flexural tension

        BMD = calculate_BMD(shear_stress_profile)
        flex_comp, flex_tens = flexural_stress_diagram(BMD, I)
        max_flex_tens = max(flex_tens, key = lambda x: abs(x[1]))[1]
        FOS_flex_tens = safety_factor(max_flex_tens, "tensile")

        if abs(FOS_flex_tens) < min_safety_factors["flexural tension"]:
            min_safety_factors["flexural tension"] = abs(FOS_flex_tens)
        
        max_flex_comp = max(flex_comp, key = lambda x: abs(x[1]))[1]
        FOS_flex_comp = safety_factor(abs(max_flex_comp), "compressive")

        if abs(FOS_flex_comp) < min_safety_factors["flexural compression"]:
            min_safety_factors["flexural compression"] = abs(FOS_flex_comp)

        # cement shear stress
        reaction_forces_list = reaction_forces(loads, span)
        shear_glue_stress_profile = shear_glue_stress_diagram(calculate_shear_force(loads, reaction_forces_list, span), I, 1)
        max_cement_shear = max(shear_glue_stress_profile, key = lambda x: abs(x[1]))[1]
        FOS_cement_shear = safety_factor(max_cement_shear, "cement_shear")
        if abs(FOS_cement_shear) < min_safety_factors["cement shear"]:
            min_safety_factors["cement shear"] = abs(FOS_cement_shear)
        
        # plate buckling
        for i in range(len(I)):
            geometry = I[i][1]
            layers = I[i][3]
            FOS1, FOS2, FOS3, FOS4 = safety_factor_thin_plate(I[i], flex_comp, shear_stress_profile, 400) # I put a random a for now
            if abs(FOS1) < min_safety_factors["case 1"]:
                min_safety_factors["case 1"] = abs(FOS1)
            if abs(FOS2) < min_safety_factors["case 2"]:
                min_safety_factors["case 2"] = abs(FOS2)
            if abs(FOS3) < min_safety_factors["case 3"]:
                min_safety_factors["case 3"] = abs(FOS3)
            if abs(FOS4) < min_safety_factors["case 4"]:
                min_safety_factors["case 4"] = abs(FOS4)
        
        loads = update_loads(loads, "left")

    return min_safety_factors

if __name__ == "__main__":
    geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(88.73, 1.27), 1.27, 72.46], "A3": [(10, 1.27), 1.27, 72.46], "A4": [(83.73, 73.73), 6.27, 1.27], "A5": [(10, 73.73), 6.27, 1.27], "A6": [(0, 75), 100, 1.27]}
    #centroidal_axis = calculate_centroidal_axis(geometry)
    #print(f"Centroidal Axis (ȳ): {centroidal_axis:.4f} mm")
    #I = second_moment_of_area(geometry, centroidal_axis)
    #print(f"Second Moment of Area (I): {I:.4f} mm^4")
    loads = [[67.5, 172], [67.5, 348], [67.5, 512], [67.5, 688], [91.0, 852], [91.0, 1028]]
    #loads = [(50, 25), (100, 1275)]
    span = 1300
    #A_y, B_y = reaction_forces(loads, span)
    I = [[13690000, geometry, (0, 1300), 1]]
    min_safety_factors = simulation_safety_factors(loads, span, I)
    print(min_safety_factors)
