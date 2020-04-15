#!/usr/bin/env python

# Copyright (C) 2014 Garth N. Wells
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

""" Module for converting DOLFIN XML to XDMF/HDF5."""

import getopt
import os
import sys
import time

import xml.etree.cElementTree as etree
import h5py
import numpy

def write_xml(filename_xdmf=None, filename_h5=None, geometric_dim=None,
              num_vertices=None, num_cells=None, celltype=None):
    "Write XDMF meta-data XML files"

    xdmf_tree   = etree.Element('Xdmf')
    xdmf_tree.attrib['Version'] = '2.0'
    xdmf_tree.attrib['xmlns:xi'] = 'http://www.w3.org/2001/XInclude'

    # Bassix XDMF data
    domain = etree.SubElement(xdmf_tree, 'Domain')
    grid   = etree.SubElement(domain, 'Grid')
    grid.attrib['Name'] = 'mesh'
    grid.attrib['GridType'] = 'Uniform'

    # Set num cells and cell type
    topology = etree.SubElement(grid, 'Topology')
    topology.attrib['NumberOfElements'] = str(num_cells)
    topology.attrib['TopologyType'] = str(celltype)

    # Set HDF5 topology data
    topology_data = etree.SubElement(topology, 'DataItem')
    topology_data.attrib['Format'] = 'HDF'
    topology_data.attrib['Dimensions'] = str(num_cells) + ' ' + str(geometric_dim + 1)
    topology_data.text = filename_h5 + ':/mesh/topology'

    # Set XDMF geometry type
    gtype = {1: 'X', 2: 'XY', 3: 'XYZ'}
    geometry = etree.SubElement(grid, 'Geometry')
    geometry.attrib['GeometryType'] = gtype[geometric_dim]

    # Set HDF5 coordinate data
    geometry_data = etree.SubElement(geometry, 'DataItem')
    geometry_data.attrib['Format'] = 'HDF'
    geometry_data.attrib['Dimensions'] = str(num_vertices) + ' ' + str(geometric_dim)
    geometry_data.text = filename_h5 + ':/mesh/coordinates'

    # Write to file
    etree.ElementTree(xdmf_tree).write(filename_xdmf)


def convert(filename_xml=None, filename_xdmf=None, filename_h5=None):
    "Parse XML file and convert to XDMF/HDF5"

    # Open HDF5 file and creat mesh group
    hdf5_file = h5py.File(filename_h5, 'w')
    mesh_grp = hdf5_file.create_group("mesh")

    VERTEX_FLUSH=100000
    CELL_FLUSH=200000

    vertex_tags = ("v0", "v1", "v2", "v3")
    x_tags      = ("x", "y", "z")

    # Parse DOLFIN XML file line-by-line
    vertex_counter = 0
    cell_counter = 0
    for event, elem in etree.iterparse(filename_xml, events=('start', 'end')):
        if event == 'start':
            if elem.tag == 'mesh':
                # Get cell type and geometric dimensions
                celltype = elem.get('celltype')
                gdim = int(elem.get('dim'))
                root = elem

            elif elem.tag == 'vertices':
                # Get number of mesh vertices
                vertex_timer = time.clock()
                num_vertices = int(elem.get('size'))
                print "Number of vertices: {0}".format(num_vertices)

                # Create coordinate buffer and HDF5 dataset
                dim = min(VERTEX_FLUSH, num_vertices)
                x_buffer = numpy.zeros((dim, gdim), dtype=numpy.float64)
                x_dataset = mesh_grp.create_dataset("coordinates", (dim, gdim), \
                                                    maxshape=(num_vertices, gdim), \
                                                    dtype=numpy.float64)
                # Store vertex root (to enable clearing while reading)
                vroot = elem

            elif elem.tag == 'cells':
                # Get number of mesh cells
                cell_timer = time.clock()
                num_cells = int(elem.get('size'))
                print "Number of cells: {0}".format(num_cells)

                # Create topology buffer and HDF5 dataset
                dim = min(CELL_FLUSH, num_cells)
                t_buffer = numpy.zeros((dim , gdim + 1), dtype=numpy.uint64)
                topology_dset = mesh_grp.create_dataset("topology", (dim, gdim + 1), \
                                                        maxshape=(num_cells, gdim + 1), \
                                                        dtype=numpy.uint64)
                # Store cell root (to enable clearing while reading)
                croot = elem

        elif event == 'end':
            if elem.tag == 'vertex':
                # Process vertex
                index = int(elem.get('index'))
                local_index = index % VERTEX_FLUSH
                for i in range(gdim):
                    x_buffer[local_index, i] = float(elem.get(x_tags[i]))
                vroot.clear()
                vertex_counter += 1

                # Flush NumPy buffer to HDF5 if threhold is reached
                if (vertex_counter % VERTEX_FLUSH) == 0:
                    progress = (float(vertex_counter)/float(num_vertices))
                    print '  Processing vertices ({0:.0f}%)'.format(100*progress)
                    x_dataset.resize(vertex_counter, axis=0)
                    x_dataset[vertex_counter-VERTEX_FLUSH:vertex_counter] = x_buffer[:]

            elif elem.tag == celltype:
                # Process cell
                index = int(elem.get('index'))
                local_index = index % CELL_FLUSH
                for i in range(gdim + 1):
                    t_buffer[local_index, i] = int(elem.get(vertex_tags[i]))
                croot.clear()
                cell_counter += 1

                # Flush NumPy buffer to HDF5 if threhold is reached
                if (cell_counter % CELL_FLUSH) == 0:
                    progress = (float(cell_counter)/float(num_cells))
                    print '  Processing cells ({0:.0f}%)'.format(100*progress)
                    topology_dset.resize(cell_counter, axis=0)
                    topology_dset[cell_counter-CELL_FLUSH:cell_counter] = t_buffer[:]

            elif elem.tag == "vertices":
                # Flush remaining coordinates to HDF5
                x_buffer_size = vertex_counter % VERTEX_FLUSH
                x_dataset.resize(vertex_counter, axis=0)
                x_dataset[vertex_counter-x_buffer_size:vertex_counter] = x_buffer[:x_buffer_size]

                progress = (float(vertex_counter)/float(num_vertices))
                print '  Processing vertices ({0:.0f}%)'.format(100*progress)
                vtime = time.clock() - vertex_timer
                if vtime > 0.0:
                    print '  Vertex processing rate (k/s): {0}'.format(int((num_vertices/vtime)/1000))

            elif elem.tag == "cells":
                # Flush remaining topology data to HDF5
                t_buffer_size = cell_counter % CELL_FLUSH
                topology_dset.resize(cell_counter, axis=0)
                topology_dset[cell_counter-t_buffer_size:cell_counter] = t_buffer[:t_buffer_size]

                progress = (float(cell_counter)/float(num_cells))
                print '  Processing cells ({0:.0f}%)'.format(100*progress)
                ctime = time.clock() - cell_timer
                if ctime > 0.0:
                    print '  Cell processing rate (k/s): {0}'.format(int((num_cells/ctime)/1000))

    # Add cell type
    topology_dset.attrs['celltype'] = celltype

    # close HDF5 file
    hdf5_file.close()

    # Write meta-data to XML file
    write_xml(filename_xdmf=filename_xdmf, filename_h5=filename_h5, \
              geometric_dim=gdim, num_vertices=num_vertices, \
              num_cells=num_cells, celltype=celltype)

def main(argv):
    "Main function"

    print "Convert DOLFIN XML mesh file to XDMF/HDF5"

    # Get command-line arguments
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["help", "input=", "output="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    # Get options
    iformat = None
    oformat = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()

    # Check that we got two filenames
    if not len(args) == 2:
        usage()
        sys.exit(2)

    # Get filenames
    ifile = args[0]
    ofile = args[1]

    # Get output filename and extension
    ofilename, ofile_extension = os.path.splitext(ofile)

    # Can only convert to XDMF/HDF5
    if ofile_extension != ".xdmf":
        raise RuntimeError("Unable to convert to format \'{0}\'. Output file extension must be \'xdmf\'.".format(ofile_extension))

    # HDF5 storage filename
    filename_h5 = ofilename + '.h5'

    convert(filename_xml=ifile, filename_xdmf=ofile, filename_h5=filename_h5)


def usage():
    "Display usage"
    print """\
Usage: dolfin-convert-xdmf input.xml output.xdmf
"""

if __name__ == "__main__":
    main(sys.argv[1:])