
# geometry = {A1: [anchor, width, height], A2: [anchor, width, height], ...}


# calculating centroidal axis
def areas(geometry):
    # areas = {A1: area1, A2: area2, ...}
    print(geometry.items())
    areas = {}
    for key, value in geometry.items():
        width = value[1]
        height = value[2]
        area = width * height
        areas[key] = area
    return areas

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



if __name__ == "__main__":
    geometry = {"A1": [(0, 75), 100, 1.27], "A2": [(10, 73.73), 6.27, 1.27], "A3": [(83.73, 73.73), 6.27, 1.27], "A4": [(10, 1.27), 1.27, 72.46], "A5": [(88.73, 1.27), 1.27, 72.46], "A6": [(10, 0), 80, 1.27]}
    centroidal_axis = calculate_centroidal_axis(geometry)
    print(f"Centroidal Axis (ȳ): {centroidal_axis:.4f} mm")
    I = second_moment_of_area(geometry, centroidal_axis)
    print(f"Second Moment of Area (I): {I:.4f} mm^4")
