from makemesh import MakeMesh

def test1_builtin_mesh():
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    MakeMesh(domain)
    return 1


def test2_import_mesh():
    domain = {"importfile": "inlet.xml"}
    MakeMesh(domain)
    return 1


def main():

    red = '\033[91m'
    green = '\033[92m'
    blue = '\033[94m'
    bold = '\033[1m'

    flag = test1_builtin_mesh()
    if flag:
        print blue+bold+"Test1: "+bold+green+"PASS"
    else:
        print blue+bold+"Test1: "+bold+red+"FAIL"

    flag = test2_import_mesh()
    if flag:
        print blue+bold+"Test2: "+bold+green+"PASS"
    else:
        print blue+bold+"Test2: "+bold+red+"FAIL"


if __name__ == '__main__':
    main()
