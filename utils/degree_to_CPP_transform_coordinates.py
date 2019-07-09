import numpy as np

center_lon = 18.183754
center_lat = -97.796768


def cpp_transform(new_lat, new_lon):
    """
    :param new_lat: in decimal degree
    :param new_lon: in decimal degree
    :return: meter
    """

    r_earth = 6378206.4
    deg_2_rad = np.pi / 180.0
    org_lat = center_lat  # SLAM0 and SFEA0
    org_lon = center_lon
    x_dis = r_earth * (new_lon - org_lon) * deg_2_rad * np.cos(org_lat * deg_2_rad)
    y_dis = new_lat * deg_2_rad * r_earth

    return x_dis, y_dis

#lat_deg = [30, 29, 29, 29, 29, 29, 27]
#lat_cen = [1.6, 27, 46.1, 28.8, 21.4, 18.6, 34.8]
#lon_deg = [90, 91, 93, 94, 94, 94, 97]
#lon_cen = [6.8, 20.3, 20.6, 55.1, 43.5, 47.6, 13]

lat_deg = [28, 27, 27]
lat_cen = [25.6, 50.2, 34.8]
lon_deg = [96, 97, 97]
lon_cen = [19.8, 2.3, 13]

lat = [lat_deg[i] + x / 60. for i, x in enumerate(lat_cen)]
lon = [- lon_deg[i] - x / 60. for i, x in enumerate(lon_cen)]

a = [cpp_transform(lat[i], lon[i]) for i in range(len(lat_deg))]


test_node_x = [x[0] for x in a]
test_node_y = [x[1] for x in a]

print "==========================================="
print test_node_x

print "==========================================="
print test_node_y
