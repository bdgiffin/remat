import sys
import os
import math
import numpy as np
import ezdxf
import pyexodus

def read_dxf(filename,layer_names):
    # attempt to open the DXF file
    try:
        doc = ezdxf.readfile(filename)
    except IOError:
        print(f"Not a DXF file or a generic I/O error.")
        sys.exit(1)
    except ezdxf.DXFStructureError:
        print(f"Invalid or corrupted DXF file.")
        sys.exit(2)

    # get the main model space object containing all entities
    msp = doc.modelspace()

    # get all lines
    lines = msp.query("LINE")
        
    # get all lines on specified layers
    if (len(layer_names) == 0):
        # if not specified on the command line, try to use all layers in the file
        layer_names = [None]*len(doc.layers)
        for i, layer in enumerate(doc.layers):
            layer_names[i] = layer.dxf.name

    # truncate layer list to include only those layers that contain lines
    new_layer_names = []
    for layer_name in layer_names:
        if (len(lines.query('*[layer=="'+layer_name+'"]')) > 0):
            new_layer_names.append(layer_name)
    layer_names = new_layer_names
                    
    # get all lines on specified layers
    line_count  = 0
    layer_lines = [None]*len(layer_names)
    layer_connectivity = [None]*len(layer_names)
    for i, layer_name in enumerate(layer_names):
        layer_lines[i] = lines.query('*[layer=="'+layer_name+'"]')
        layer_connectivity[i] = np.zeros(shape=(len(layer_lines[i]),2),dtype=int)
        line_count = line_count + len(layer_lines[i])
                        
    # populate data
    points = np.zeros(shape=(2*line_count, 3))
    line_count = 0
    for j, connectivity in enumerate(layer_connectivity):
        for i, line in enumerate(layer_lines[j]):
            points[2*(line_count+i)][:]   = line.dxf.start
            points[2*(line_count+i)+1][:] = line.dxf.end
            connectivity[i][:]            = [2*(line_count+i),2*(line_count+i)+1]
        line_count = line_count + len(layer_lines[j])
                                
    # get min/max coordinates, and set default proximity tolerance
    min_coordinates = np.array([+sys.float_info.max]*3)
    max_coordinates = np.array([-sys.float_info.max]*3)
    for point in points:
        for i in range(3):
            min_coordinates[i] = min(min_coordinates[i],point[i])
        for i in range(3):
            max_coordinates[i] = max(max_coordinates[i],point[i])
    box_dimensions = max_coordinates - min_coordinates
    tolerance = max(box_dimensions)*1.0e-6
                                            
    # consolidate duplicate nodes
    point_ids = -np.ones(2*line_count,dtype=int)
    num_unique = 0
    unique_points = np.zeros(shape=(2*line_count, 3))
    for i, ipoint in enumerate(points):
        if (point_ids[i]==-1):
            point_ids[i] = num_unique
            unique_points[num_unique,:] = ipoint
            num_unique += 1
        for j in range(i,len(point_ids)):
            if (np.linalg.norm(ipoint-points[j,:]) < tolerance):
                point_ids[j] = point_ids[i]
    points = np.resize(unique_points,(num_unique,3))
                                                            
    # map all original point ids in the connectivity array to be expressed in terms of unique node ids
    for connectivity in layer_connectivity:
        for i in range(connectivity.shape[0]):
            for j in range(connectivity.shape[1]):
                connectivity[i][j] = point_ids[connectivity[i][j]]
                                                                        
    # look for any points that intersect with existing lines, and subdivide into multiple elements

    # loop over all nodes
    for i in range(points.shape[0]):
        # loop over all elements
        for k, connectivity in enumerate(layer_connectivity):
            for j in range(layer_connectivity[k].shape[0]):
                if (not (layer_connectivity[k][j,0]==i) and not (layer_connectivity[k][j,1]==i)):
                    xa = points[layer_connectivity[k][j,0],:]
                    xb = points[layer_connectivity[k][j,1],:]
                    dxbar = xa-points[i,:]
                    dx = xb-xa
                    xi = -np.dot(dxbar,dx)/np.dot(dx,dx)
                    d = np.linalg.norm(dxbar+xi*dx)
                    if (d < tolerance):
                        if ((xi < 1.0-1.0e-6) and (xi > 1.0e-6)):
                            # create new element
                            layer_connectivity[k] = np.append(layer_connectivity[k],[layer_connectivity[k][j,:]],axis=0)
                            layer_connectivity[k][j,1]  = i
                            layer_connectivity[k][-1,0] = i

    # return values from read_dxf
    return points, layer_connectivity, layer_names


def refine_mesh(points,layer_connectivity,layer_element_sizes):
    # loop over all elements
    for i, connectivity in enumerate(layer_connectivity):
        element_size = layer_element_sizes[i]
        Nelements = connectivity.shape[0]
        for j in range(Nelements):
            xa = points[connectivity[j,0],:]
            xb = points[connectivity[j,1],:]
            dx = xb-xa
            d = np.linalg.norm(dx)
            if ((element_size > 0.0) and (d > element_size)):
                # sub-divide into smaller elements, as needed
                Nsub_elements = math.ceil(d/element_size)
                new_points = np.zeros(shape=(Nsub_elements-1,3))
                new_connectivity = np.zeros(shape=(Nsub_elements-1,2),dtype=int)
                for k in range(Nsub_elements-1):
                    new_points[k,:] = xa + (dx*(k+1))/float(Nsub_elements)
                    new_connectivity[k,0] = points.shape[0]+k
                    new_connectivity[k,1] = points.shape[0]+k+1
                new_connectivity[Nsub_elements-2,1] = connectivity[j,1]
                layer_connectivity[i][j,1] = points.shape[0]
                points = np.append(points,new_points,axis=0)
                layer_connectivity[i] = np.append(layer_connectivity[i],new_connectivity,axis=0)

    # return values from refine_mesh
    return points, layer_connectivity


def write_exodus(filename,points,layer_connectivity,layer_names):
    # count up all layers with non-zero numbers of lines
    num_blocks = 0
    num_elements = 0
    for connectivity in layer_connectivity:
        num_blocks = num_blocks + 1
        num_elements = num_elements + connectivity.shape[0]

    # create Exodus file
    exo = pyexodus.exodus(file=filename, mode='w', array_type='numpy', title=None, numDims=3, numNodes=points.shape[0], numElems=num_elements, numBlocks=num_blocks, numNodeSets=0, numSideSets=0, io_size=0, compression=None)

    # put node coordinates
    exo.put_coords(xCoords=points[:,0],yCoords=points[:,1],zCoords=points[:,2])

    # put element block data
    num_blocks = 0
    for connectivity in layer_connectivity:
        num_blocks = num_blocks + 1
        exo.put_elem_blk_info(id=num_blocks, elemType='BAR2', numElems=connectivity.shape[0], numNodesPerElem=connectivity.shape[1], numAttrsPerElem=0)
            
        # put element block 1 connectivity (shift indicies by 1 to ensure compatibility with 1-based indexing in Exodus files, whereas Python uses 0-based indexing)
        exo.put_elem_connectivity(id=num_blocks, connectivity=connectivity, shift_indices=1, chunk_size_in_mb=128)

    # close the Exodus object (avoids error message)
    exo.close()
