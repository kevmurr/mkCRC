import numpy as np
from stl import mesh
import stl
from utils import get_total_correction2
from utils import mk_mesh
from utils import mk_vertices_faces
from utils import mk_box
from utils import find_mins_maxs
from utils import translate
from utils import copy_obj

##########################
def run_mk_crc(name,dir,size_mlls=50, coeff_x=0.00015, f_x=250, coeff_y=0.00015, f_y=350, N_lensletts=10, delta_PMMAmod_17_5=1.72e-6 , pillarwidth=10,scale_spacer=0, postheight=20, px_size=1, factor_size_x=1, factor_size_y=1, factor_phase=1.0, bool_box=True, bool_post=True):
    print("Box:",bool_box)
    print("Post:",bool_post)
    flip_y = False
    xdim_o=size_mlls # in µm
    ydim_o=size_mlls # in µm
    # gets the total correction without knowing the lenslet number. Unit is in µm.
    (pupil_size_arr_x,total_corr_arr_x)=get_total_correction2(c=coeff_x,size_sfl_mll=xdim_o,factor_size_sfl=factor_size_x,f_sfl=f_x,delta=delta_PMMAmod_17_5,pxsize=px_size)
    (pupil_size_arr_y,total_corr_arr_y)=get_total_correction2(c=coeff_y,size_sfl_mll=ydim_o,factor_size_sfl=factor_size_y,f_sfl=f_y,delta=delta_PMMAmod_17_5,pxsize=px_size)
    array=np.zeros((pupil_size_arr_y.shape[0],pupil_size_arr_x.shape[0]))
    for i1 in range(0,array.shape[0]):
        for i2 in range(0,array.shape[1]):
            array[i1,i2]=total_corr_arr_x[i2]+total_corr_arr_y[i1]
    #now dividing by the number of refracting surfaces to get each surface
    array=(array-np.amin(array))/(2*N_lensletts)
    if factor_phase!=None:
        array=array*factor_phase
    if flip_y==True:
        array=np.flip(array,axis=0)
    maxheight=np.amax(array)
    #make_box_template
    xdim=np.float(xdim_o*factor_size_x)
    xdimwall=pillarwidth
    ydim=np.float(ydim_o*factor_size_y)
    ydimwall=ydim
    #total z height of each lenslet
    zdim=np.float(2*maxheight+scale_spacer)


    pxsize=xdim/array.shape[1]
    sampling_z=100
    N_px_long=np.amax(array.shape)
    N_px_short=np.amin(array.shape)
    print("Generating lenslett")
    z_mid=scale_spacer/2#(scale_spacer+maxheight)/2
    front_vertices,front_faces,array_used=mk_vertices_faces(array,xdim=xdim,ydim=ydim,z_mid=z_mid,N_px_long=N_px_long)
    print(len(front_vertices),len(front_faces))
    #This is new with numpy-stl

    lenslett_pre=mk_mesh(vertices=front_vertices,faces=front_faces)
    #lenslett_pre=mesh.form_mesh(front_vertices,front_faces)
    print("made mesh")
    boxoverlap = 0
    if bool_box == True:
        # overlap side boxes, factor compared to xdim
        yoffset = -ydim * (array_used.shape[0] - 1) / (array_used.shape[0])
        xmin1 = -pillarwidth - boxoverlap * xdim
        ymin1 = yoffset + 0.001
        zmin1 = -z_mid - maxheight
        xmax1 = boxoverlap * xdim
        ymax1 = 0
        zmax1 = z_mid + maxheight

        xmin2 = (array_used.shape[1] - 1) * xdim / array_used.shape[1] - boxoverlap * xdim
        ymin2 = yoffset + 0.001
        zmin2 = -z_mid - maxheight
        xmax2 = xdim + pillarwidth
        ymax2 = 0
        zmax2 = z_mid + maxheight

        box1 = mk_box(xmin1, ymin1, zmin1, xmax1, ymax1, zmax1)
        box2 = mk_box(xmin2, ymin2, zmin2, xmax2, ymax2, zmax2)

    if bool_post == True:
        yoffset = -ydim * (array_used.shape[0] - 1) / (array_used.shape[0])
        xminp = -pillarwidth - boxoverlap * xdim
        yminp = yoffset - postheight
        zminp = -z_mid - maxheight
        xmaxp = xdim + pillarwidth
        ymaxp = yoffset
        zmaxp = z_mid + maxheight
        post = mk_box(xminp, yminp, zminp, xmaxp, ymaxp, zmaxp)

    # Using an existing stl file:
    main_body = lenslett_pre

    # rotate along Y
    #main_body.rotate([0.0, 0.5, 0.0], math.radians(90))



    if bool_box==False and bool_post==False:
        combined1=main_body
    if bool_box==True and bool_post==False:
        combined1 = mesh.Mesh(np.concatenate([main_body.data, box1.data, box2.data]))
    if bool_box==False and bool_post==True:
        combined1 = mesh.Mesh(np.concatenate([main_body.data, post.data]))
    if bool_box==True and bool_post==True:
        combined1 = mesh.Mesh(np.concatenate([main_body.data, box1.data, box2.data , post.data] ))

    #Now making N lensletts
    tile_offset=2*maxheight+2*z_mid
    minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(combined1)
    w1 = maxx - minx
    l1 = maxy - miny
    h1 = maxz - minz
    print(tile_offset,w1,l1,h1)
    copies = copy_obj(combined1, (tile_offset, l1, h1), 1, N_lensletts, 1)
    CRC=mesh.Mesh(np.concatenate([copy.data for copy in copies]))
    CRC.save('{0}{1}.stl'.format(dir,name), mode=stl.Mode.ASCII)  # save as ASCII
    width_total=2*pillarwidth+xdim
    height_total=ydim
    depth_total=N_lensletts*(2*maxheight+2*z_mid)
    print("expected size (x,y,z): ",(width_total,height_total,depth_total))
    print("Saved {0}.stl".format(name))
    print("Done")
    print("zmid",z_mid)
    return(CRC)