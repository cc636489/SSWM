

from fenics import RectangleMesh, Mesh, Point


def make_mesh(domain_dict):

    # check if the key and the value are entered right.
    key = list(domain_dict.keys())
    if len(key) != 1:
        raise TypeError("enter more than one domain! not allowed! Exit program...")

    # check if the value is as expected.
    value = domain_dict.get(key[0])
    if key[0] == "rectangle":
        if len(value) != 6:
            raise TypeError("expecting 6 domain parameters! Exit program...")
    elif key[0] == "importfile":
        if not isinstance(value, str):
            raise TypeError("expecting strings of the mesh full name! Exit program...")
    else:
        raise TypeError("not an implemented type of domain! Exit program...")

    # Make the mesh.
    if key[0] == "rectangle":
        mesh_object = RectangleMesh(Point(value[0], value[1]), Point(value[2], value[3]), value[4], value[5])
    elif key[0] == "importfile":
        mesh_object = Mesh(value)
    else:
        raise TypeError("not an implemented type of domain! Exit program...")

    return mesh_object
