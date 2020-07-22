import numpy as np
from stl import mesh
import math
import stl
import scipy.constants as sc


def get_total_correction2(c,size_sfl_mll,factor_size_sfl,f_sfl,delta,pxsize):
    a_h=factor_size_sfl*size_sfl_mll/2
    c_mod=c*10**9/(f_sfl**3)
    phase_line_h_axis=np.linspace(-int(factor_size_sfl*size_sfl_mll)/2,int(factor_size_sfl*size_sfl_mll)/2,int(factor_size_sfl*size_sfl_mll/pxsize))
    phase_line_h=c_mod*(phase_line_h_axis**3 - a_h**2 * phase_line_h_axis)
    total_correction_pre=-phase_line_h-np.amin(phase_line_h)
    L=np.amax(total_correction_pre)/2
    n=np.amax(phase_line_h_axis)
    line=((-L)/(n*(1+1/np.sqrt(3))))*phase_line_h_axis-L/(np.sqrt(3)+1)
    if c_mod>0:
        total_correction=total_correction_pre+line
    else:
        total_correction=total_correction_pre-line-(L-L/(np.sqrt(3)+2))
    wavelength=sc.c*sc.h/(17480*sc.e)
    corr_fac=1E6*wavelength/(2*np.pi*delta)
    total_correction=total_correction*corr_fac
    return(phase_line_h_axis,total_correction)

def mk_mesh(vertices,faces):
    meshs=mesh.Mesh(np.zeros(len(faces),dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
                           for j in range(3):
                            #print(i,f[j])
                            #print(vertices[f[j]])
                            meshs.vectors[i][j]=vertices[f[j]]
    return(meshs)


def mk_vertices_faces(input_array,xdim,ydim, N_px_long, x_mid=0, y_mid=0, z_mid=0):
    array_new = np.zeros((N_px_long, N_px_long))
    array_new[array_new.shape[0] - input_array.shape[0]:array_new.shape[0], 0:input_array.shape[1]] = input_array
    array_new = array_new[:, :-1]
    vertices = []
    faces = []
    i = 0
    # make front vertices
    for i1 in range(0, array_new.shape[0], 1):
        for i2 in range(0, array_new.shape[1], 1):
            vertices.append((array_new[i1, i2] + z_mid, (i2) * xdim / array_new.shape[1] + x_mid,
                             (-i1) * ydim / array_new.shape[0] + y_mid))
    # make front faces
    for i1 in range(0, array_new.shape[0], 1):
        for i2 in range(0, array_new.shape[1], 1):
            if i1 >= 1:
                if i2 >= 0 and i2 < (array_new.shape[1] - 1):
                    faces.append((i, i - array_new.shape[1], i - array_new.shape[1] + 1))
                if i2 >= 1:
                    faces.append((i - 1, i - array_new.shape[1], i))
            i += 1
    # other side
    array_flip = array_new
    offset = len(vertices)
    i1 = 0
    i2 = 0
    i = 0
    # make back vertices
    for i1 in range(0, array_new.shape[0], 1):
        for i2 in range(0, array_new.shape[1], 1):
            vertices.append((-array_flip[i1, i2] - z_mid, (i2) * xdim / array_flip.shape[1] + x_mid,
                             (-i1) * ydim / array_flip.shape[0] + y_mid))

    # make back faces
    for i1 in range(0, array_new.shape[0], 1):
        for i2 in range(0, array_new.shape[1], 1):
            if i1 >= 1:
                if i2 >= 0 and i2 < (array_new.shape[1] - 1):
                    faces.append(
                        np.flip(np.add(offset, (i, i - array_new.shape[1], i - array_new.shape[1] + 1)), axis=0))
                if i2 >= 1:
                    faces.append(np.flip(np.add(offset, (i - 1, i - array_new.shape[1], i)), axis=0))
            i += 1

    # make_sides
    for i in range(0, array_flip.shape[1] - 1):
        if i < (array_flip.shape[1] - 1):
            #       #top
            faces.append(np.flip((i, i + 1, offset + i), axis=0))
            faces.append(np.flip((i + 1, offset + 1 + i, offset + i), axis=0))
            # bottom
            short = (array_new.shape[0] - 1) * array_new.shape[1] + i
            faces.append((short, short + 1, short + offset))
            faces.append(np.flip((short + offset, short + offset + 1, short + 1), axis=0))
        if i == (array_flip.shape[1] - 1):
            # top
            faces.append(np.flip((i, i + 1, offset + i), axis=0))
            # bottom
            short = (array_new.shape[0] - 1) * array_new.shape[1] + i
            faces.append((short, short + 1, short + offset))

    i = 0
    for i in range(0, array_flip.shape[0] - 1, 1):
        # left
        faces.append(
            np.flip(((array_new.shape[1]) * i, offset + i * (array_new.shape[1]), (i + 1) * (array_new.shape[1])),
                    axis=0))
        faces.append(np.flip(((array_new.shape[1]) * (i + 1), offset + i * (array_new.shape[1]),
                              (i + 1) * (array_new.shape[1]) + offset), axis=0))

    for i in range(1, array_flip.shape[0], 1):
        #    #right
        faces.append(
            ((array_new.shape[1]) * i - 1, offset + (array_new.shape[1]) * i - 1, (array_new.shape[1]) * (i + 1) - 1))
        faces.append(((array_new.shape[1]) * (i + 1) - 1, offset + (array_new.shape[1]) * i - 1,
                      offset + ((array_new.shape[1]) * (i + 1)) - 1))
        # faces.append(np.flip(((array_new.shape[1])*i-1,offset+(array_new.shape[1])*i-1,(array_new.shape[1])*(i+1)-1),axis=0))
        # faces.append(np.flip(((array_new.shape[1])*(i+1)-1,offset+(array_new.shape[1])*i-1,offset+((array_new.shape[1])*(i+1))-1),axis=0))

    faces = np.array(faces)
    faces = np.flip(faces, axis=1)
    return (vertices, faces, array_new)

def mk_box(xmin, ymin, zmin, xmax, ymax, zmax):

    # Create 3 faces of a cube
    data = np.zeros(12, dtype=mesh.Mesh.dtype)

    # Top of the cube
    data['vectors'][0] = np.array(np.flip(([[xmin, ymax, zmax],
                                            [xmax, ymin, zmax],
                                            [xmin, ymin, zmax]]), axis=0))
    data['vectors'][1] = np.array(np.flip(([[xmax, ymin, zmax],
                                            [xmin, ymax, zmax],
                                            [xmax, ymax, zmax]]), axis=0))
    # Right face
    data['vectors'][2] = np.array(np.flip(([[xmax, ymin, zmin],
                                            [xmax, ymin, zmax],
                                            [xmax, ymax, zmin]]), axis=0))
    data['vectors'][3] = np.array([[xmax, ymax, zmax],
                                   [xmax, ymin, zmax],
                                   [xmax, ymax, zmin]])
    # Left face
    data['vectors'][4] = np.array([[xmin, ymin, zmin],
                                   [xmax, ymin, zmin],
                                   [xmax, ymin, zmax]])
    data['vectors'][5] = np.array(np.flip(([[xmin, ymin, zmin],
                                            [xmin, ymin, zmax],
                                            [xmax, ymin, zmax]]), axis=0))
    # Altleft face
    data['vectors'][6] = np.array([[xmin, ymax, zmax],
                                   [xmax, ymax, zmax],
                                   [xmin, ymax, zmin]])
    data['vectors'][7] = np.array([[xmax, ymax, zmax],
                                   [xmax, ymax, zmin],
                                   [xmin, ymax, zmin]])
    # Altrightface
    data['vectors'][8] = np.array(np.flip(([[xmin, ymax, zmin],
                                            [xmin, ymax, zmax],
                                            [xmin, ymin, zmax]]), axis=0))
    data['vectors'][9] = np.array(np.flip(([[xmin, ymax, zmin],
                                            [xmin, ymin, zmax],
                                            [xmin, ymin, zmin]]), axis=0))
    # Alttop
    data['vectors'][10] = np.array(np.flip(([[xmin, ymax, zmin],
                                             [xmin, ymin, zmin],
                                             [xmax, ymin, zmin]]), axis=0))
    data['vectors'][11] = np.array([[xmax, ymin, zmin],
                                    [xmin, ymax, zmin],
                                    [xmax, ymax, zmin]])
    # rotate
    for i in range(12):
        data["vectors"][i] = np.roll(data["vectors"][i], 1, axis=1)
    # Bottom of the cube
    box = mesh.Mesh(data)
    return (box)

def find_mins_maxs(obj):
    minx = maxx = miny = maxy = minz = maxz = None
    for p in obj.points:
        # p contains (x, y, z)
        if minx is None:
            minx = p[stl.Dimension.X]
            maxx = p[stl.Dimension.X]
            miny = p[stl.Dimension.Y]
            maxy = p[stl.Dimension.Y]
            minz = p[stl.Dimension.Z]
            maxz = p[stl.Dimension.Z]
        else:
            maxx = max(p[stl.Dimension.X], maxx)
            minx = min(p[stl.Dimension.X], minx)
            maxy = max(p[stl.Dimension.Y], maxy)
            miny = min(p[stl.Dimension.Y], miny)
            maxz = max(p[stl.Dimension.Z], maxz)
            minz = min(p[stl.Dimension.Z], minz)
    return minx, maxx, miny, maxy, minz, maxz


def translate(_solid, step, padding, multiplier, axis):
    if axis == 'x':
        items = [0, 3, 6]
    elif axis == 'y':
        items = [1, 4, 7]
    elif axis == 'z':
        items = [2, 5, 8]
    for p in _solid.points:
        # point items are ((x, y, z), (x, y, z), (x, y, z))
        for i in range(3):
            p[items[i]] += (step * multiplier) + (padding * multiplier)


def copy_obj(obj, dims, num_rows, num_cols, num_layers):
    w, l, h = dims
    copies = []
    for layer in range(num_layers):
        for row in range(num_rows):
            for col in range(num_cols):
                # skip the position where original being copied is
                if row == 0 and col == 0 and layer == 0:
                    continue
                _copy = mesh.Mesh(obj.data.copy())
                # pad the space between objects by 10% of the dimension being
                # translated
                if col != 0:
                    translate(_copy, w, 0., col, 'x')
                if row != 0:
                    translate(_copy, l, l / 10., row, 'y')
                if layer != 0:
                    translate(_copy, h, h / 10., layer, 'z')
                copies.append(_copy)
    return copies