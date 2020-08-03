from PyQt5 import QtWidgets,uic,QtGui
import sys
import mk_crc
import numpy as np
from gui import parameters
from utils import mk_box
from mpl_toolkits import mplot3d
import time

class App(QtWidgets.QMainWindow,parameters.Ui_MainWindow):
    def __init__(self):
        super(App,self).__init__()
        self.setupUi(self)

        #uic.loadUi("gui/parameters.ui",self)
        self.RunButton=self.findChild(QtWidgets.QPushButton,"RunButton")
        self.RunButton.clicked.connect(self.run_stuff)
        self.Button_SetDir=self.findChild(QtWidgets.QPushButton,"Button_SetDir")
        self.Button_SetDir.clicked.connect(self.set_dir)
        self._directory=""
        self._CRC_mesh=None
        self.ButtonPreviewMesh=self.findChild(QtWidgets.QPushButton,"ButtonPreviewMesh")
        self.ButtonPreviewMesh.clicked.connect(self.preview_mesh)
        self.show()
        self.first_mesh_plot()

    def first_mesh_plot(self):
        #remaining bug: Plot doesnt really work
        mpl=self.previewWidget.canvas
        mpl.ax.clear()
        mpl.ax.mouse_init()
        #mpl.ax.imshow(zz)
        self.box=self.create_default_rect()
        polys=mplot3d.art3d.Poly3DCollection(self.box.vectors)
        polys.set_edgecolor("k")
        polys.set_facecolor("c")
        collection=mpl.ax.add_collection3d(polys)
        scale=self.box.points.flatten(-1)
        mpl.ax.set_xlabel("Length (µm)")
        mpl.ax.set_ylabel("Width (µm)")
        mpl.ax.set_zlabel("Height (µm)")
        mpl.ax.auto_scale_xyz(scale,scale,scale)
        mpl.draw()

    def preview_mesh(self):
        if self._CRC_mesh==None:
            print("Run the program first to display mesh!")
            self.error_dialog_run=QtWidgets.QErrorMessage()
            self.error_dialog_run.showMessage("Run the program first to generate mesh of CRC!")
        else:
            time1=time.time()
            mpl=self.previewWidget.canvas
            mpl.ax.clear()
            mpl.ax.mouse_init()
            polys=mplot3d.art3d.Poly3DCollection(self._CRC_mesh.vectors)
            polys.set_edgecolor("k")
            polys.set_facecolor("c")
            collection = mpl.ax.add_collection3d(polys)
            scale = self._CRC_mesh.points.flatten(-1)
            mpl.ax.auto_scale_xyz(scale, scale, scale)
            mpl.ax.set_xlabel("Length (µm)")
            mpl.ax.set_ylabel("Width (µm)")
            mpl.ax.set_zlabel("Height (µm)")
            mpl.draw()
            time2=time.time()
            if time2-time1>1:
                print("WARNING 3D viewer seems to be too inefficient for this model!")
                self.error_dialog_efficiency=QtWidgets.QErrorMessage()
                self.error_dialog_efficiency.showMessage("3D viewer is inefficient and generated model seems to be large. Using an external viewer is advised.")
    def set_dir(self):
        self._directory = str(QtWidgets.QFileDialog.getExistingDirectory(self,"Select directory"))+"/"

    def create_default_rect(self):
        rect=mk_box(-1,-1,-1,1,1,1)
        return(rect)

    def run_stuff(self):
        #reading the entries
        name=str(self.LineEdit_Name.text())
        N_lensletts=int(self.LineEdit_N_lensletts.text())
        mllsize = float(self.LineEdit_mllsize.text())
        coeff_h = float(self.LineEdit_coeff_h.text())
        coeff_v= float(self.LineEdit_coeff_v.text())
        f_h = float(self.LineEdit_f_h.text())
        f_v = float(self.LineEdit_f_v.text())
        delta_mat=float(self.LineEdit_delta_mat.text())
        scale_spacer = float(self.LineEdit_scale_spacer.text())
        pillarwidth= float(self.LineEdit_pillarwidth.text())
        postheight= float(self.LineEdit_postheight.text())
        factor_size_h = 1#float(self.LineEdit_factor_size_h.text())
        factor_size_v = 1#float(self.LineEdit_factor_size_v.text())
        factor_phase = float(self.LineEdit_factor_phase.text())
        pxsize = float(self.LineEdit_pxsize.text())

        bool_post=bool(self.Box_P.checkState())
        bool_box = bool(self.Box_BB.checkState())
        print("This is the directory",self._directory)
        self._CRC_mesh=mk_crc.run_mk_crc(name=name,dir=self._directory,N_lensletts=N_lensletts,size_mlls=mllsize, coeff_x=coeff_h, coeff_y=coeff_v, f_x=f_h, f_y=f_v, bool_post=bool_post, bool_box=bool_box, delta_PMMAmod_17_5=delta_mat, scale_spacer=scale_spacer,pillarwidth=pillarwidth, postheight=postheight, factor_size_x=factor_size_h, factor_size_y=factor_size_v, factor_phase=factor_phase, px_size=pxsize)
app=QtWidgets.QApplication(sys.argv)
window=App()
app.exec_()