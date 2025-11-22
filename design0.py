import copy
import math

# geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}


''' Calculating centroidal axis '''

# Calculating area of each component
def areas(geometry):

    # Parameters:
    # Geometry: dictionary of components with their location and dimensions.  Must be inputted by user.
    #   geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ..., An: [anchor, width, height]}
    #       An is the nth component of the cross-section.  Components are listed from bottom to top of cross-section.
    #       Anchor = (x, y) --> coordinate of bottom-left corner of component (mm).  y = 0 occurs at the bottom of the cross-section.
    #       width = width of component (mm)
    #       height = height of component (mm)
    # Returns a list of each component's area
    #   areas_list = [area1, area2, ...]

    areas_list = []  # defining empty list to store areas
    for key, value in geometry.items():
        width = value[1]  # width of component
        height = value[2]  # height of component
        area = width * height
        areas_list.append(area)  # add area to list
    return areas_list  # areas_list = [area1, area2, ...]

# Finding centroids of each component
def find_centroids(geometry):

    # Parameters:
    # Geometry: dictionary of components with their location and dimensions
    #   geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ..., An: [anchor, width, height]}
    #       An is the nth component of the cross-section.  Components are listed from bottom to top of cross-section.
    #       Anchor = (x, y) --> coordinate of bottom-left corner of component (mm).  y = 0 occurs at the bottom of the cross-section.
    #       width = width of component (mm)
    #       height = height of component (mm)
    # Returns a list of each component's centroid y-coordinate. y = 0 occurs at the bottom of the cross-section.
    #   centroids = [c_y1, c_y2, ..., c_yn]

    centroids = [] # defining empty list to store centroids
    for value in geometry.values():
        anchor = value[0]  # anchor = (x, y)
        height = value[2]  # height of component
        c_y = anchor[1] + height / 2 # centroid of rectangular component = height of component / 2 + y-coordinate of bottom of component
        centroids.append((c_y))  # add centroid to list
    return centroids  # centroids = [c_y1, c_y2, ..., c_yn]

# Calculating centroidal axis of entire cross-section
def calculate_centroidal_axis(geometry):

    # Parameters:
    # geometry: dictionary of components with their location and dimensions
    #   geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ..., An: [anchor, width, height]}
    #       An is the nth component of the cross-section.  Components are listed from bottom to top of cross-section.
    #       Anchor = (x, y) --> coordinate of bottom-left corner of component (mm).  y = 0 occurs at the bottom of the cross-section.
    #       width = width of component (mm)
    #       height = height of component (mm)
    # Returns the y-coordinate of the centroidal axis of the entire cross-section. y = 0 occurs at the bottom of the cross-section.

    areas_list = areas(geometry)  # calculate areas of each component
    total_area = sum(areas_list)  # calculate total area of cross-section
    centroids = find_centroids(geometry)  # calculate centroids of each component
    
    y_bar = 0  # initialize y_bar

    for i in range(len(centroids)): # for each component's centroid, add area * centroid to y_bar
        y_bar += areas_list[i] * centroids[i]
    y_bar /= total_area  # divide by total area to get centroidal axis
    return y_bar  # return centroidal axis of cross-section

''' Calulating second moment of area (I) '''

# Calculating I for one rectangular component about its own centroidal axis
def second_moment_area_component(component):
    
    # Parameters:
    # component: list of component dimensions
    #   component = [anchor, width, height]
    #       anchor = (x, y) --> coordinate of bottom-left corner of component (mm).  y = 0 occurs at the bottom of the cross-section.
    #       width = width of component (mm)
    #       height = height of component (mm)
    # Returns the second moment of area of the component about its own centroidal axis
    
    I = component[1] * component[2]**3 / 12  # component[1] = width, component[2] = height
    return I  # return I of component

# Calculating second moment of area of the entire cross-section about the centroidal axis
def second_moment_of_area(geometry, y_bar):

    # Parameters:
    # geometry: dictionary of components with their location and dimensions
    #   geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ..., An: [anchor, width, height]}
    #       An is the nth component of the cross-section.  Components are listed from bottom to top of cross-section.
    #       Anchor = (x, y) --> coordinate of bottom-left corner of component.  y = 0 occurs at the bottom of the cross-section.
    #       width = width of component
    #       height = height of component
    # y_bar: y-coordinate of centroidal axis of entire cross-section. y = 0 occurs at the bottom of the cross-section.
    # Returns the second moment of area of the entire cross-section about the centroidal axis

    I_total = 0  # initialize I_total

    list_areas = areas(geometry)  # list of areas of each component
    components = list(geometry.values())  # list of each component's location and dimensions
    centroids = find_centroids(geometry)  # list of each component's centroid's y-coordinate
    
    for i in range(len(components)):  # for each component, calculate I using parallel axis theorem
        I = second_moment_area_component(components[i])  # I of component about its own centroidal axis
        d_sq = (centroids[i] - y_bar)**2  # square of distance from component's centroid to cross-section's centroidal axis
        I_total += I + list_areas[i] * d_sq  # add I of component + area * d^2 to I_total (parallel axis theorem)
    
    return I_total  # return I of entire cross-section

''' Reaction forces '''

# loads = [[load1, position1], [load2, position2], ...]
# Calculate reaction forces at supports A and B
def reaction_forces(loads, span):

    # Parameters:
    # loads: list of lists of loads and their positions
    #   loads = [[load1, position1], [load2, position2], ...,]
    #       load = magnitude of nth load (N)
    #       position = location of nth load along span (mm).  x = 0 is at left end of bridge span.
    # span: length of bridge (mm)
    # Returns a list of tuples of the reaction forces at support A and support B
    #   reaction_forces = [(A_y, loc_A_y), (B_y, loc_B_y)]

    A_y = 0  # initialize reaction force at A
    B_y = 0  # initialize reaction force at B

    loc_A_y = span / 2 - 600  # support A is located 600 mm left of center of span
    loc_B_y = span / 2 + 600  # support B is located 600 mm right of center of span

    total_load = 0  # initialize total load

    for load in loads:  # find the total applied load by summing all loads
        total_load += load[0]
    
    # find B_y first by taking moment about A
    B_y = 0  # initialize reaction force at B

    for load in loads:
    # sign convention: counterclockwise moment is positive.
    # sum(moments about A) = 0
    # 0 = sum(moments due to loads to the left of A) - sum(moments due to loads to the right of A) + moment due to B_y
    # B_y * 1200 = sum(moments due to loads to the right of A) - sum(moments due to loads to the left of A)
    # All applied loads point downwards.  Signs of moments are flipped when B_y is isolated.
            
            B_y -= load[0]*(loc_A_y - load[1]) 
            # if load is to the left of A, loc_A_y - load[1] is positive, so moments due to loads to the left of A are subtracted.
            # if load is to the right of A, loc_A_y - load[1] is negative, so moments due to loads to the right of A are added.

    B_y = B_y / 1200  # divide by distance from A to B (1200 mm) to isolate B_y (N)
    A_y = total_load - B_y  # A_y + B_y = total applied load, so A_y = total_load - B_y (N)
    return [(A_y, loc_A_y), (B_y, loc_B_y)]  # returns a list of tuples describing the magnitude (N) and location (mm) of the reaction forces along the span


''' Updating location of loads to simulate train movement '''

# Increment the position of the loads by 1 mm
def update_loads(loads, direction):
    
    # Parameters:
    # loads: list of loads (N) and their positions along the span of the bridge (mm)
    #   loads = [[load1, position1], [load2, position2], ...]
    # direction: string type, either "left" or "right", indicates which direction the train is moving
    # Returns the list of loads with each load's location incremented 1 mm left or right.
    #   loads = [[load1, new position1], [load2, new position2], ...]

    if direction == "right":
        for load in loads:
            load[1] += 1  # shift each train wheel right 1 mm
    elif direction == "left":
        for load in loads:
            load[1] -= 1  # shift each train wheel left 1 mm
    return loads  # returns the list of loads with incremented locations


''' Finding SFD '''

# Finds SFD along the bridge for a given load distribution
def calculate_shear_force(loads, reaction_forces, span):

    # Parameters:
    # loads: list of loads (N) and their positions along the span of the bridge (mm)
    #   loads = [[load1, position1], [load2, position2], ...]
    # reaction_forces: list of tuples (magnitude of force (N), location (mm)) describing reaction force at each support
    # span: length of bridge (mm)
    # Returns a list of lists of location-shear force pairs
    #   shear_force_diagram = [[x1, V1], [x2, V2], ...]

    shear_force_diagram = []  # initialize SFD
    V = 0  # initialize shear force

    # sign convention for V: upwards is positive, downwards is negative.

    for x in range(span + 1):  # run through every mm increment along the bridge and check if there is a force there.  If so, update V.
        for force in reaction_forces:  # check if there is a reaction force at that x-value
            if force[1] == x and x != 1200:
                V += force[0]  # add reaction force, since reaction forces point upwards
        for load in loads:  # check if there is a load at that x-value
            if load[1] == x:
                V -= load[0]  # subtract load, since loads point downwards
        shear_force_diagram.append([x, V])  # add V (N) and its location (mm) to the diagram
    
    return shear_force_diagram  # shear_force_diagram = [[x1, V1], [x2, V2], ...] --> return a list of lists of shear force at each mm along the bridge

''' Finding shear stress diagram '''

# Find maximum Q for a given geometry using the area between y_bar and the bottom shear stress-free surface
def calculate_Qmax(geometry):

    # Parameters
    # geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}
    # Returns maximum Q for the cross-section (mm^3)

    geo_below_ybar = copy.deepcopy(geometry)  # initialize geometry below y_bar
    y_bar = calculate_centroidal_axis(geometry)

    for component in copy.deepcopy(geo_below_ybar).values():
        if component[0][1] >= y_bar:
            geo_below_ybar.pop(list(geo_below_ybar.keys())[list(geo_below_ybar.values()).index(component)])  # remove components above y_bar
        if component[2] + component[0][1] > y_bar:
            component[2] = y_bar - component[0][1]  # cut the heights of the webs to y-bar
    
    Q = 0  # initialize Q

    for component in geo_below_ybar.values():
        area = component[1] * component[2]
        d = (component[0][1] + component[2]/2) - y_bar  # distance between y-bar and centroid of area below y-bar
        Q += area * d
    
    return Q # (mm^3)

# Finds maximum shear stress along span of bridge, accounting for changing I values for different sections.
def shear_stress_diagram(shear_force_diagram, I_list, b):

    # Parameters:
    # shear_force_diagram: list of shear forces at each location along the span
    #   shear_force_diagram = [[x1, V1], [x2, V2], ...]
    #       x = location along span (mm)
    #       V = shear force at that location
    # I_list: list of lists of cross-sectional features (cross-sections change for different sections of the span)
    #   I_list = [[I1, geometry1, (start, end)], [I2, geometry2, (start, end)]...]
    #       I = second moment of area of cross-section (mm^4)
    #       geometry = geometry of cross-section (dictionary)
    #           geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ..., An: [anchor, width, height]}
    #       (start, end) = position along span where cross-section begins and position where it ends (mm)
    # b: width of the beam's cross section at the centroidal axis (mm)
    # Returns a list of lists of location-tau pairs

    Q_maxs = []  # initialize list of Q_maxs for each cross-section

    for i in range(len(I_list)):  # for each cross-section, calculate Q_max and add it to the list Q_maxs
        Q = calculate_Qmax(I_list[i][1])
        Q_maxs.append(Q)
    
    shear_stresses_diagram = []  # initialize shear_stresses_diagram
    SFD = shear_force_diagram[:]  # make a copy of shear_force_diagram to avoid accidentally altering values

    for i in range(len(I_list)):  # run through each cross-section segment
        for x in range(I_list[i][2][0], I_list[i][2][1]):  # for every x in the range of that cross-section's segment
            shear_stress = SFD[x][1] * Q_maxs[i] / (I_list[i][0] * b)  # calculate shear stress
            shear_stresses_diagram.append([x, shear_stress])  # add location and corresponding shear stress to shear_stresses_diagram
    return shear_stresses_diagram  # shear_stresses_diagram = [[x1, tau1], [x2, tau2], ...]

        
''' Finding bending moment diagram '''

# Find bending moment diagram
def calculate_BMD(SFD):
    
    # Parameters:
    # SFD: shear force diagram = [[x1, V1], [x2, V2], ...]
    # Returns bending moment diagram = [[x1, M1], [x2, M2], ...]

    bending_moment_diagram = []  # initialize BMD
    M = 0
    for x in SFD:
        M += x[1]  # moment = area under SFD (Nmm).  When incrementing by 1 mm, area = V * 1 --> x[1] = V --> delta M = x[1] * 1
        bending_moment_diagram.append([x[0], M])

    return bending_moment_diagram  # bending_moment_diagram = [[x1, M1], [x2, M2], ...]


''' Find flexural stress diagram '''

# find flexural stress diagram --- flexural stress for every point along span of beam
def flexural_stress_diagram(BMD, I_list):

    # Parameters:
    # BMD: bending moment diagram
    #   BMD = [[x1, M1], [x2, M2], ...]
    # I_list: list of cross-sections
    #   I_list = [[I, geometry, (start, end), layers], ...]
    #       geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ..., An: [anchor, width, height]}
    #       Assume geometry is organized from bottom to top of cross-section.  The last component, An, is the topmost component.
    # Returns 2 lists: One of location-flexural compression pairs and one of location-flexural tension pairs

    flexural_compression_diagram = []
    flexural_tension_diagram = []
    
    for i in range(len(I_list)):  # run through each cross-section

        upper_component_dimensions = list(geometry.values())[-1]
        height = upper_component_dimensions[0][1] + upper_component_dimensions[2]
        y_bar = calculate_centroidal_axis(geometry)
        y_compression = y_bar - height  # distance from y_bar to top of cross-section
        y_tension = y_bar  # distance from y_bar to bottom of cross-section
        
        for x in range(I[i][2][0], I[i][2][1]+1):  # run through the range of that segment
            
            sigma_compression = BMD[x][1] * y_compression / I_list[i][0]  # sigma = MY/I (MPa)
            sigma_tension = BMD[x][1] * y_tension / I_list[i][0]

            flexural_compression_diagram.append([x, sigma_compression])
            flexural_tension_diagram.append([x, sigma_tension])

    return flexural_compression_diagram, flexural_tension_diagram  # flexural_compression_diagram = [[x1, sigma_c1], [x2, sigma_c2], ...] and flexural_tension_diagram [[x1, sigma_t1], [x2, sigma_t2], ...]


''' Find plate buckling stresses '''

# calculate plate buckling stress for a given case (1, 2, 3, or 4)

# geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}
def plate_buckling_stress(geometry, case, layers, a = None):

    # Parameters:
    # geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}
    # case: 1, 2, 3, or 4
    # layers: the number of matboard layers making the top flange (1, 2, 3, ...)
    # a = spacing between diaphragms (mm)
    # Returns the failure plate buckling stress for a given case (MPa)

    t = 1.27 * layers # thickness of flange

    # case 1: buckling of compressive flange between webs
    if case == 1:
        b = list(geometry.values())[-1][1] - (list(geometry.values())[-(layers+2)][0][0] + 1.27) * 2 # width of flange between webs
        k = 4
        sigma = k * (3.14159**2) * 4000 * (t / b)**2 / (12 * (1 - 0.2**2))
        return sigma
    
    # case 2: buckling of the tips of the compressive flange
    if case == 2:
        b = list(geometry.values())[-(layers+2)][0][0] # length of flange tip beyond web
        k = 0.425
        sigma = k * (3.14159**2) * 4000 * (t / b)**2 / (12 * (1 - 0.2**2))
        print(sigma)
        return sigma
    
    # case 3: buckling of the webs due to the flexural stresses
    if case == 3:
        t = 1.27 # thickness of web
        h = list(geometry.values())[-(layers+3)][0][1] + list(geometry.values())[-(layers+3)][2] - calculate_centroidal_axis(geometry) # height of web above centroidal axis
        k = 6
        sigma = k * (3.14159**2) * 4000 * (t / h)**2 / (12 * (1 - 0.2**2))
        return sigma
    
    # case 4: buckling of the webs due to the shear stresses
    if case == 4:
        t = 1.27 # thickness of web
        h = list(geometry.values())[-(layers+3)][2] # height of web
        k = 5
        tau = (k * (3.14159**2) * 4000 * ((t / h)**2 + (t / a)**2)) / (12 * (1 - 0.2**2))
        return tau



''' Find glue shear stress '''

# glue shear stress diagram
def shear_glue_stress_diagram(shear_force_diagram, I_list, type):
    
    # Parameters:
    # shear_force_diagram = [[x1, V1], [x2, V2], ...]
    # I_list = [[I, geometry, (start, end), layers], ...]
    #   geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}
    # type: A string indicating the type of glue shear stress, either "glue tabs" (shear stress at junction between webs and flange) or "sheets" (shear stress between layered sheets)
    # Returns shear_glue_stresses_diagram = [[x1, sigma_glue1], [x2, sigma_glue2], ...]

    shear_glue_stresses_diagram = []

    for i in range(len(I)):

        geometry = I_list[i][1]
        layers = I_list[i][3]

        if type == "glue tabs":

            if layers == 1:
                level = 1
                upper_component_dimensions = list(geometry.values())[-level]
                d = upper_component_dimensions[0][1] + upper_component_dimensions[2]/2 - calculate_centroidal_axis(geometry) # distance from centroid of area above glue line to centroidal axis
                A = upper_component_dimensions[1] * upper_component_dimensions[2]
                b = list(geometry.values())[-2][1]*2 # the second-to-top component is the glue tab component, so its width*2 is the total glue line length
            
            if layers == 2:
                level = 2
                upper_component_dimensions = list(geometry.values())[-level]
                second_upper_component_dimensions = list(geometry.values())[-(level-1)]
                d = calculate_centroidal_axis({"upper": upper_component_dimensions, "second_upper": second_upper_component_dimensions}) - calculate_centroidal_axis(geometry) # distance from centroid of area above glue line to centroidal axis
                A = (upper_component_dimensions[1] * upper_component_dimensions[2]) + (second_upper_component_dimensions[1] * second_upper_component_dimensions[2])
                b = list(geometry.values())[-3][1]*2 # the third-to-top component is the glue tab component, so its width*2 is the glue line length
            
            if layers == 3:
                level = 3
                upper_component_dimensions = list(geometry.values())[-level]
                second_upper_component_dimensions = list(geometry.values())[-(level-1)]
                third_upper_component_dimensions = list(geometry.values())[-(level-2)]
                d = calculate_centroidal_axis({"upper": upper_component_dimensions, "second_upper": second_upper_component_dimensions, "third_upper": third_upper_component_dimensions}) - calculate_centroidal_axis(geometry) # distance from centroid of area above glue line to centroidal axis
                A = (upper_component_dimensions[1] * upper_component_dimensions[2]) + (second_upper_component_dimensions[1] * second_upper_component_dimensions[2]) + (third_upper_component_dimensions[1] * third_upper_component_dimensions[2])
                b = list(geometry.values())[-4][1]*2 # the fourth-to-top component is the glue tab component, so its width*2 is the glue line length
        
        if type == "sheets":

            if layers == 1: # if there's only one layer in the cross-section, there's no "sheet" type glue stress in that segment of the beam
                continue
            
            level = 2
            upper_component_dimensions = list(geometry.values())[-level]
            second_upper_component_dimensions = list(geometry.values())[-(level-1)]
            A = upper_component_dimensions[1] * upper_component_dimensions[2]
            d = upper_component_dimensions[0][1] + upper_component_dimensions[2]/2 - calculate_centroidal_axis(geometry) # distance from centroid of area above glue line to top shear-stress free surface
            b = second_upper_component_dimensions[1]
        
        Q = A * d

        for x in range(I[i][2][0], I[i][2][1] + 1):
            shear_glue_stress = shear_force_diagram[x][1] * Q / (I[i][0] * b)  # (MPa)
            shear_glue_stresses_diagram.append([x, shear_glue_stress])

    return shear_glue_stresses_diagram  # shear_glue_stresses_diagram = [[x1, sigma_glue1], [x2, sigma_glue2], ...]


''' Safety factors '''

# safety factors for tension, compression, shear, and glue shear
def safety_factor(applied_stress, type):

    # Parameters:
    # applied_stress = applied stress of interest (MPa)
    # type = string indicating what type of stress ("tensile", "compressive", "shear", "cement_shear")
    # returns safety factor

    allowable_stresses = {"tensile": 30, "compressive": 6, "shear": 4, "cement_shear": 2} 
    return allowable_stresses[type] / applied_stress # FOS = capacity / demand


# safety factors for thin plate buckling
def safety_factor_thin_plate(I_n, flexural_compression_diagram, shear_stress_diagram, a = None):

    # Parameters:
    # I_n: One cross-section from I_list
    #   I_n = [I_value, geometry, (start, end), layers]
    # flexural_compression_diagram = [[x1, sigma_c1], [x2, sigma_c2], ...] --> used for cases 1-3
    # shear_stress_diagram = [[x1, tau1], [x2, tau2], ...] --> used for case 4
    # Returns the 4 safety factors for the 4 plate buckling cases for a given cross-section

    geometry = I_n[1]
    layers = I_n[3]
    start, end = I_n[2][0], I_n[2][1]

    case_1_failure = plate_buckling_stress(geometry, 1, layers)
    max_flex_comp = max(flexural_compression_diagram[start:end], key = lambda x: abs(x[1]))[1] # use max flexural compression
    FOS1 = case_1_failure / abs(max_flex_comp)

    case_2_failure = plate_buckling_stress(geometry, 2, layers)
    FOS2 = case_2_failure / abs(max_flex_comp)

    case_3_failure = plate_buckling_stress(geometry, 3, layers)
    list_component_geometry = list(geometry.values())
    y_bar = calculate_centroidal_axis(geometry)
    y_max = list_component_geometry[-1][0][1] + list_component_geometry[-1][2] - y_bar
    y_top_of_web = list_component_geometry[-(layers+3)][0][1] + list_component_geometry[-(layers+3)][2] - y_bar

    max_flex_comp_case_3 = max_flex_comp/y_max * y_top_of_web

    FOS3 = case_3_failure / abs(max_flex_comp_case_3)

    case_4_failure = plate_buckling_stress(geometry, 4, layers, a)
    max_shear = max(shear_stress_diagram[start:end], key = lambda x: abs(x[1]))[1]
    FOS4 = case_4_failure / abs(max_shear) 

    return FOS1, FOS2, FOS3, FOS4


def initialize_loads():
    freight = 500/(1.38+1.1+1)
    return [[freight*1.38,0], [freight*1.38,176], [freight*1, 340], [freight*1,516], [freight*1.1, 680], [freight*1.1,856]]

def simulation_safety_factors(loads, span, I, b=2.54):
    min_safety_factors = {"flexural tension": math.inf, "flexural compression": math.inf, "shear": math.inf, "cement shear": math.inf, "case 1": math.inf, "case 2": math.inf, "case 3": math.inf, "case 4": math.inf} # [flexural tension, flexural compression, shear, cement shear, plate buckling case 1, case 2, case 3, case 4]

    for x in range(span - loads[-1][1]):
        # shear stress
        shear_stress_profile = shear_stress_diagram(calculate_shear_force(loads, reaction_forces(loads, span), span), I, b)
        max_shear = max(shear_stress_profile, key = lambda x: abs(x[1]))[1]
        FOS_shear = safety_factor(max_shear, "shear")
        if abs(FOS_shear) < min_safety_factors["shear"]:
            min_safety_factors["shear"] = abs(FOS_shear)
        
        # flexural stress and flexural tension

        BMD = calculate_BMD(calculate_shear_force(loads, reaction_forces(loads, span), span))
        flex_comp, flex_tens = flexural_stress_diagram(BMD, I)
        max_flex_tens = max(flex_tens, key = lambda x: abs(x[1]))[1]
        FOS_flex_tens = safety_factor(max_flex_tens, "tensile")

        if abs(FOS_flex_tens) < min_safety_factors["flexural tension"]:
            min_safety_factors["flexural tension"] = abs(FOS_flex_tens)
        max_flex_comp = max(flex_comp, key = lambda x: abs(x[1]))[1]
        print(max_flex_comp)
        FOS_flex_comp = safety_factor(abs(max_flex_comp), "compressive")

        if FOS_flex_comp < min_safety_factors["flexural compression"]:
            min_safety_factors["flexural compression"] = FOS_flex_comp

        # cement shear stress
        reaction_forces_list = reaction_forces(loads, span)
        shear_glue_stress_profile = shear_glue_stress_diagram(calculate_shear_force(loads, reaction_forces_list, span), I, "glue tabs")
        max_cement_shear = max(shear_glue_stress_profile, key = lambda x: abs(x[1]))[1]
        FOS_cement_shear = safety_factor(max_cement_shear, "cement_shear")
        if abs(FOS_cement_shear) < min_safety_factors["cement shear"]:
            min_safety_factors["cement shear"] = abs(FOS_cement_shear)

        shear_glue_stress_profile = shear_glue_stress_diagram(calculate_shear_force(loads, reaction_forces_list, span), I, "sheets")
        if shear_glue_stress_profile != []:
            max_cement_shear = max(shear_glue_stress_profile, key = lambda x: abs(x[1]))[1]
            FOS_cement_shear = safety_factor(max_cement_shear, "cement_shear")
            if abs(FOS_cement_shear) < min_safety_factors["cement shear"]:
                min_safety_factors["cement shear"] = abs(FOS_cement_shear)
        
        # plate buckling
        for i in range(len(I)):
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

        BMD = calculate_BMD(calculate_shear_force(loads, reaction_forces(loads, span), span))
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
        shear_glue_stress_profile = shear_glue_stress_diagram(calculate_shear_force(loads, reaction_forces_list, span), I, "glue tabs")
        max_cement_shear = max(shear_glue_stress_profile, key = lambda x: abs(x[1]))[1]
        FOS_cement_shear = safety_factor(max_cement_shear, "cement_shear")
        if abs(FOS_cement_shear) < min_safety_factors["cement shear"]:
            min_safety_factors["cement shear"] = abs(FOS_cement_shear)

        shear_glue_stress_profile = shear_glue_stress_diagram(calculate_shear_force(loads, reaction_forces_list, span), I, "sheets")
        if shear_glue_stress_profile != []:
            max_cement_shear = max(shear_glue_stress_profile, key = lambda x: abs(x[1]))[1]
            FOS_cement_shear = safety_factor(max_cement_shear, "cement_shear")
            if abs(FOS_cement_shear) < min_safety_factors["cement shear"]:
                min_safety_factors["cement shear"] = abs(FOS_cement_shear)
        
        # plate buckling
        for i in range(len(I)):
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

def simulation_safety_factors_across_bridge(span, I, b=None):
    
    FOS_shear_diagram = []
    FOS_flex_tens_diagram = []
    FOS_flex_comp_diagram = []
    FOS_cement_shear_diagram_glue_tabs = []
    FOS_cement_shear_diagram_sheets = []
    FOS_case_1_diagram = []
    FOS_case_2_diagram = []
    FOS_case_3_diagram = []
    FOS_case_4_diagram = []
    
    loads = initialize_loads()

    for x in range(span - loads[-1][1]+1):
        # shear stress
        shear_stress_profile = shear_stress_diagram(calculate_shear_force(loads, reaction_forces(loads, span), span), I, b)
        max_shear = max(shear_stress_profile, key = lambda x: abs(x[1]))[1]
        FOS_shear = safety_factor(max_shear, "shear")
        FOS_shear_diagram.append([x, FOS_shear])
        
        # flexural stress and flexural tension

        BMD = calculate_BMD(calculate_shear_force(loads, reaction_forces(loads, span), span))
        flex_comp, flex_tens = flexural_stress_diagram(BMD, I)
        max_flex_tens = max(flex_tens, key = lambda x: abs(x[1]))[1]
        FOS_flex_tens = safety_factor(max_flex_tens, "tensile")
        FOS_flex_tens_diagram.append([x, FOS_flex_tens])
        
        max_flex_comp = max(flex_comp, key = lambda x: abs(x[1]))[1]
        FOS_flex_comp = safety_factor(abs(max_flex_comp), "compressive")
        FOS_flex_comp_diagram.append([x, FOS_flex_comp])

        # cement shear stress
        reaction_forces_list = reaction_forces(loads, span)
        shear_glue_stress_profile_glue_tabs = shear_glue_stress_diagram(calculate_shear_force(loads, reaction_forces_list, span), I, "glue tabs")
        max_cement_shear_glue_tabs = max(shear_glue_stress_profile_glue_tabs, key = lambda x: abs(x[1]))[1]
        FOS_cement_shear_glue_tabs = safety_factor(max_cement_shear_glue_tabs, "cement_shear")
        FOS_cement_shear_diagram_glue_tabs.append([x, FOS_cement_shear_glue_tabs])

        shear_glue_stress_profile_sheets = shear_glue_stress_diagram(calculate_shear_force(loads, reaction_forces_list, span), I, "sheets")
        if shear_glue_stress_profile_sheets != []:
            max_cement_shear_sheets = max(shear_glue_stress_profile_sheets, key = lambda x: abs(x[1]))[1]
            FOS_cement_shear_sheets = safety_factor(max_cement_shear_sheets, "cement_shear")
            FOS_cement_shear_diagram_sheets.append([x, FOS_cement_shear_sheets])
        
        
        # plate buckling
        for i in range(len(I)):
            FOS1, FOS2, FOS3, FOS4 = safety_factor_thin_plate(I[i], flex_comp, shear_stress_profile, 400) # I put a random a for now
            FOS_case_1_diagram.append([x, FOS1])
            FOS_case_2_diagram.append([x, FOS2])
            FOS_case_3_diagram.append([x, FOS3])
            FOS_case_4_diagram.append([x, FOS4])
        
        loads = update_loads(loads, "right")
        # print(FOS_shear_diagram)
    return FOS_shear_diagram, FOS_flex_tens_diagram, FOS_flex_comp_diagram, FOS_cement_shear_diagram_glue_tabs, FOS_cement_shear_diagram_sheets, FOS_case_1_diagram, FOS_case_2_diagram, FOS_case_3_diagram, FOS_case_4_diagram


def calculate_envs(load_mag, span, geometry, level):
    react = reaction_forces(load_mag, span)
    sf = calculate_shear_force(load_mag, react, span)

    bm = calculate_BMD(sf)

    y_bar = calculate_centroidal_axis(geometry)
    I = second_moment_of_area(geometry, y_bar)
    I_list = [[I, geometry, (0,span), level]]
    comp, tens = flexural_stress_diagram(bm, I_list)

    glue = shear_glue_stress_diagram(sf, I_list, level)
    
    return [sf, bm, comp, tens, glue]
    
if __name__ == "__main__":
    geometry = {"A1": [(10, 0), 80, 1.27], "A2": [(88.73, 1.27), 1.27, 72.46], "A3": [(10, 1.27), 1.27, 72.46], "A4": [(83.73, 73.73), 6.27, 1.27], "A5": [(10, 73.73), 6.27, 1.27], "A6": [(0, 75), 100, 1.27]}
    #centroidal_axis = calculate_centroidal_axis(geometry)
    #print(f"Centroidal Axis (ȳ): {centroidal_axis:.4f} mm")
    #I = second_moment_of_area(geometry, centroidal_axis)
    #print(f"Second Moment of Area (I): {I:.4f} mm^4")
    #loads = [[67.5, 172], [67.5, 348], [67.5, 512], [67.5, 688], [91.0, 852], [91.0, 1028]]
    loads = [[400/6, 172], [400/6, 348], [400/6, 512], [400/6, 688], [400/6, 852], [400/6, 1028]]
    freight = 500/(1.38+1.1+1)
    kN = [[freight*1.38,0], [freight*1.38,176], [freight*1, 340], [freight*1,516], [freight*1.1, 680], [freight*1.1,856]]
    #loads = [(50, 25), (100, 1275)]
    span = 1260
    b = 2.54
    #A_y, B_y = reaction_forces(loads, span)
    #I = [[418480.7, geometry, (0, 1200), 1]]
    #min_safety_factors = simulation_safety_factors(loads, span, I)
    #print(min_safety_factors)
    #geometry2 = {"A1": [(10, 0), 80, 1.27], "A2": [(10, 1.27), 1.27, 72.46], "A3": [(85, 1.27), 1.27, 72.46], "A4": [(10, 73.73), 6.27, 1.27], "A5": [(83.73, 73.73), 6.27, 1.27], "A6": [(0, 75), 100, 1.27], "A7": [(0, 77.54), 100, 1.27]}
    #I = [[second_moment_of_area(geometry2, calculate_centroidal_axis(geometry2)), geometry2, (0, 1260), 2]]
    #min_safety_factors = simulation_safety_factors(kN, span, I)
    #print(min_safety_factors)
    #print(calculate_centroidal_axis(geometry2))
    #print(second_moment_of_area(geometry2, calculate_centroidal_axis(geometry2)))

    #geometry3 = {"A1": [(10, 0), 80, 1.27], "A2": [(10, 1.27), 1.27, 72.46], "A3": [(85, 1.27), 1.27, 72.46], "A4": [(10, 73.73), 6.27, 1.27], "A5": [(83.73, 73.73), 6.27, 1.27], "A6": [(0, 75), 100, 1.27], "A7": [(0, 77.54), 100, 1.27], "A8": [(0, 78.81), 100, 1.27]}
    #I = [[second_moment_of_area(geometry3, calculate_centroidal_axis(geometry3)), geometry3, (0, 1260), 3]]
    #min_safety_factors = simulation_safety_factors(kN, span, I)
    #print(min_safety_factors)
    #print(calculate_centroidal_axis(geometry3))
    #print(second_moment_of_area(geometry3, calculate_centroidal_axis(geometry3)))

    geometry4 = {"A2": [(10, 1.27), 1.27, 72.46], "A3": [(85, 1.27), 1.27, 72.46], "A4": [(10, 73.73), 6.27, 1.27], "A5": [(83.73, 73.73), 6.27, 1.27], "A6": [(0, 75), 100, 1.27], "A7": [(0, 77.54), 100, 1.27], "A8": [(0, 78.81), 100, 1.27]}
    I = [[second_moment_of_area(geometry4, calculate_centroidal_axis(geometry4)), geometry4, (0, 1260), 3]]
    min_safety_factors = simulation_safety_factors(kN, span, I)
    print(min_safety_factors)
    print(calculate_centroidal_axis(geometry4))
    print(second_moment_of_area(geometry4, calculate_centroidal_axis(geometry4)))
    