import numpy as np

input_file = "fort.14_closed_small_gulf"
H0 = 10
center_lon = 18.183754
center_lat = -97.796768
output_file = "fort.14_closed_small_gulf_cartesian_" + str(H0) + "m"
node = 0


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


with open(input_file) as f:
    for i, this_line in enumerate(f):
        if i == 1:
            node = int(this_line.split()[1])

with open(input_file, 'r') as f, open(output_file, 'w+') as f1:
    for i, this_line in enumerate(f):
        tmp = this_line.split()
        if 2 <= i <= node + 1:
            x_distance, y_distance = cpp_transform(new_lat=float(tmp[2]), new_lon=float(tmp[1]))
            tmp[1] = str(x_distance)
            tmp[2] = str(y_distance)
            if float(tmp[3]) < H0:
                tmp[3] = str(H0)
        tmp_str = " ".join(tmp)
        f1.write(tmp_str)
        f1.write('\n')
