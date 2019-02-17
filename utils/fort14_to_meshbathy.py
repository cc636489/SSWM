from netCDF4 import *
import numpy as np


# Boundary Mark "3": open boundary node
# Boundary Mark "1": land boundary node
# Boundary Mark "2": island boundary node

inputfile = "fort.14_inlet_fine_mesh"
outputfile = "inlet_fine_mesh.msh"
bathymetryfile = "inlet_dgswem_compare_finer.nc"
DEBUG = True
Cartesian = True
center_lat = 18.183754
center_lon = -97.796768


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


######################################
# Reading necessary information.
######################################
# initialize reading
nelem = 0
nnode = 0
OpenbcNodes = []
OpenbcPointer = []
MainlandNodes = []
MainlandPointer = []
IslandNodes = []
IslandPointer = []

# read basic info
with open(inputfile) as finput:

    for i, thisline in enumerate(finput):
        # read node and element
        if i == 1:
            nelem = int(thisline.split()[0])
            nnode = int(thisline.split()[1])
            if DEBUG:
                print "nelem, nnode = ", nelem, ", ", nnode

        # read lat lon min max
        if i == 2:
            words = [float(x) for x in thisline.split()]
            lonmin = lonmax = words[1]
            latmin = latmax = words[2]
            bathymin = bathymax = words[3]

        # reset lat lon min max
        if (i > 2) and (i <= nnode+1):
            words = [float(x) for x in thisline.split()]
            if words[1] < lonmin:
                lonmin = words[1]
            if words[1] > lonmax:
                lonmax = words[1]
            if words[2] < latmin:
                latmin = words[2]
            if words[2] > latmax:
                latmax = words[2]
            if words[3] < bathymin:
                bathymin = words[3]
            if words[3] > bathymax:
                bathymax = words[3]

        # read bc info
        if "Number of open boundaries" in thisline or "NOPE" in thisline:
            nopenbc = int(thisline.split()[0])

            if DEBUG:
                print "number of openbc = ", nopenbc

        if "open boundary nodes" in thisline or "NETA" in thisline:
            TotalOpenbcNodes = int(thisline.split()[0])

            if DEBUG:
                print "total openbc nodes = ", TotalOpenbcNodes

        if "for open boundary" in thisline or "no. of nodes" in thisline:
            OpenbcPointer.append(i+1)
            OpenbcPointer.append(i+int(thisline.split()[0]))

            if DEBUG:
                print "openbc pointer is = ", OpenbcPointer

        if "Number of land boundaries" in thisline or "NBOU" in thisline:
            nlandbc = int(thisline.split()[0])

            if DEBUG:
                print "number of landbc = ", nlandbc

        if "land boundary nodes" in thisline or "NVEL" in thisline:
            TotalLandbcNodes = int(thisline.split()[0])

            if DEBUG:
                print "total landbc nodes = ", TotalLandbcNodes

        if "for land boundary" in thisline or "number of nodes" in thisline:
            words = thisline.split()

            if words[1] == '0':
                MainlandPointer.append(i+1)
                MainlandPointer.append(i+int(words[0]))

                if DEBUG:
                    print "mainland pointer is = ", MainlandPointer

            elif words[1] == '1':
                IslandPointer.append(i+1)
                IslandPointer.append(i+int(words[0]))

                if DEBUG:
                    print "island pointer is = ", IslandPointer

            else:
                print "Can't transfer this type of land boundaries!"

# read openbc nodes
OpenbcNNODEList = [0]

for k in range(nopenbc):
    OpenbcNNODEList.append(OpenbcPointer[2*k+1]-OpenbcPointer[2*k]+1)

with open(inputfile) as finput:
    for i, thisline in enumerate(finput):
        if (i > nnode+nelem+3) and (i <= nnode+nelem+3+nopenbc+TotalOpenbcNodes) and ("for open boundary" not in thisline) and ("no. of nodes" not in thisline):
            OpenbcNodes.append(int(thisline.split()[0]))

# read mainland nodes and island nodes
nMainland = len(MainlandPointer) / 2
nIsland = len(IslandPointer) / 2
MainlandNNODEList = [0]
IslandNNODEList = [0]

for k in range(nMainland):
    MainlandNNODEList.append(MainlandPointer[2*k+1]-MainlandPointer[2*k]+1)
    with open(inputfile) as finput:
        for i, thisline in enumerate(finput):
            if (i >= MainlandPointer[2 * k]) and (i <= MainlandPointer[2 * k + 1]):
                MainlandNodes.append(int(thisline.split()[0]))

for k in range(nIsland):
    IslandNNODEList.append(IslandPointer[2*k+1]-IslandPointer[2*k]+1)
    with open(inputfile) as finput:
        for i, thisline in enumerate(finput):
            if (i >= IslandPointer[2 * k]) and (i <= IslandPointer[2 * k + 1]):
                IslandNodes.append(int(thisline.split()[0]))

# calculate boundary element number
bcnelem = TotalOpenbcNodes + TotalLandbcNodes - nopenbc - nlandbc
if DEBUG:
    print 'total land bc is %d, total mainland bc is %d, total island bc is %d' % (nlandbc, nMainland, nIsland)
    print 'total open boundary node is %d, length of nodelist is %d' % (TotalOpenbcNodes, len(OpenbcNodes))
    print 'total land boundary node is %d, length of mainland is %d, length of island is %d, sum is %d' % \
          (TotalLandbcNodes, len(MainlandNodes), len(IslandNodes), len(MainlandNodes)+len(IslandNodes))


##################################
# Writing mesh output information
##################################
foutput = open(outputfile, "w")
foutput.write("$MeshFormat\n")
foutput.write("2.2 0 8\n")
foutput.write("$EndMeshFormat\n")
foutput.write("$Nodes\n")
foutput.write(str(nnode)+"\n")

# write node info
with open(inputfile) as finput:

    for i, thisline in enumerate(finput):

        if (i > 1) and (i <= nnode+1):

            # write to mesh file
            words = thisline.split()
            if Cartesian:
                x_distance, y_distance = float(words[1]), float(words[2])
            else:
                x_distance, y_distance = cpp_transform(new_lat=float(words[2]), new_lon=float(words[1]))

            string = words[0] + "  " + str(x_distance) + "  " + str(y_distance) + "  " + "0.0" + "\n"
            if DEBUG:
                print string
            foutput.write(string)

        if i > nnode+1: break

foutput.write("$EndNodes\n")
foutput.write("$Elements\n")
foutput.write(str(nelem+bcnelem)+"\n")

# write elem info
with open(inputfile) as finput:

    for i, thisline in enumerate(finput):

        if (i >= nnode+2) and (i <= nnode+nelem+1):
            words = thisline.split()
            string = words[0]+"  2  1  0  "+words[2]+"  "+words[3]+"  "+words[4]+"\n"
            if DEBUG:
                print string
            foutput.write(string)

        if i > nnode+nelem+1: break

# write open bc info
currIndex = 0
for i in range(len(OpenbcNNODEList)-1):
    currIndex += OpenbcNNODEList[i]

    for j in range(currIndex, currIndex + OpenbcNNODEList[i+1]-1):
        nelem += 1
        string = str(nelem)+"  1  1  3  "+str(OpenbcNodes[j])+"  "+str(OpenbcNodes[j+1])+"\n"

        if DEBUG:
            print string

        foutput.write(string)

# write mainland bc info
currIndex = 0
for i in range(len(MainlandNNODEList)-1):
    currIndex += MainlandNNODEList[i]

    for j in range(currIndex, currIndex + MainlandNNODEList[i+1]-1):
        nelem += 1
        string = str(nelem)+"  1  1  1  "+str(MainlandNodes[j])+"  "+str(MainlandNodes[j+1])+"\n"

        if DEBUG:
            print string

        foutput.write(string)

# write island bc info
currIndex = 0
for i in range(len(IslandNNODEList)-1):
    currIndex += IslandNNODEList[i]

    for j in range(currIndex, currIndex + IslandNNODEList[i+1]-1):
        nelem += 1
        string = str(nelem)+"  1  1  2  "+str(IslandNodes[j])+"  "+str(IslandNodes[j+1])+"\n"

        if DEBUG:
            print string

        foutput.write(string)

# close output file
foutput.write("$EndElements")
foutput.close()


########################################
# write bathymetry netcdf file
########################################
rootgrp = Dataset(bathymetryfile, "w", format="NETCDF4")

rootgrp.createDimension("nnode", nnode)

x_coord = rootgrp.createVariable("x_coord", np.double, ("nnode",))
x_coord.long_name = "x_coord"
x_coord.units = "meters"

y_coord = rootgrp.createVariable("y_coord", np.double, ("nnode",))
y_coord.long_name = "y_coord"
y_coord.units = "meters"

longitude = rootgrp.createVariable("lon", np.double, ("nnode",))
longitude.long_name = "longitude"
longitude.units = "degree"
longitude.minmax = [lonmin, lonmax]

latitude = rootgrp.createVariable("lat", np.double, ("nnode",))
latitude.long_name = "latitude"
latitude.units = "degree"
latitude.minmax = [latmin, latmax]

bathymetry = rootgrp.createVariable("bathy", np.double, ("nnode",))
bathymetry.long_name = "bathymetry"
bathymetry.units = "meters"
bathymetry.minmax = [bathymin, bathymax]

with open(inputfile) as finput:

    for i,thisline in enumerate(finput):

        if (i > 1) and (i <= nnode+1):
            words = thisline.split()
            if Cartesian:
                x_distance, y_distance = float(words[1]), float(words[2])
            else:
                x_distance, y_distance = cpp_transform(new_lat=float(words[2]), new_lon=float(words[1]))
            x_coord[i - 2] = x_distance
            y_coord[i - 2] = y_distance
            longitude[i - 2] = float(words[1])
            latitude[i - 2] = float(words[2])
            bathymetry[i - 2] = float(words[3])
            print "the bathymetry in ", i-2, "is ", bathymetry[i-2]

        if i > nnode+1: break

x_coord.minmax = [np.min(x_coord), np.max(x_coord)]
y_coord.minmax = [np.min(y_coord), np.max(y_coord)]

rootgrp.close()




