# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'parameters.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets
from mplwidget import MPL_WIDGET

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 595)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label_15 = QtWidgets.QLabel(self.centralwidget)
        self.label_15.setGeometry(QtCore.QRect(15, 330, 91, 16))
        self.label_15.setObjectName("label_15")
        self.RunButton = QtWidgets.QPushButton(self.centralwidget)
        self.RunButton.setGeometry(QtCore.QRect(25, 550, 75, 23))
        self.RunButton.setObjectName("RunButton")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(15, 208, 91, 16))
        self.label_9.setObjectName("label_9")
        self.LineEdit_mllsize = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_mllsize.setGeometry(QtCore.QRect(115, 104, 113, 20))
        self.LineEdit_mllsize.setObjectName("LineEdit_mllsize")
        self.label_19 = QtWidgets.QLabel(self.centralwidget)
        self.label_19.setGeometry(QtCore.QRect(15, 390, 91, 16))
        self.label_19.setObjectName("label_19")
        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(15, 300, 91, 16))
        self.label_12.setObjectName("label_12")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(15, 156, 91, 16))
        self.label_4.setObjectName("label_4")
        self.label_16 = QtWidgets.QLabel(self.centralwidget)
        self.label_16.setGeometry(QtCore.QRect(235, 330, 47, 13))
        self.label_16.setObjectName("label_16")
        self.frame = QtWidgets.QFrame(self.centralwidget)
        self.frame.setGeometry(QtCore.QRect(10, 94, 281, 141))
        self.frame.setObjectName("frame")
        self.label_25 = QtWidgets.QLabel(self.centralwidget)
        self.label_25.setGeometry(QtCore.QRect(235, 510, 47, 13))
        self.label_25.setObjectName("label_25")
        self.label_23 = QtWidgets.QLabel(self.centralwidget)
        self.label_23.setGeometry(QtCore.QRect(15, 480, 91, 16))
        self.label_23.setObjectName("label_23")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(15, 130, 91, 16))
        self.label_3.setObjectName("label_3")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(15, 182, 91, 16))
        self.label_7.setObjectName("label_7")
        self.frame_3 = QtWidgets.QFrame(self.centralwidget)
        self.frame_3.setGeometry(QtCore.QRect(10, 280, 281, 261))
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.LineEdit_pillarwidth = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_pillarwidth.setGeometry(QtCore.QRect(115, 360, 113, 20))
        self.LineEdit_pillarwidth.setObjectName("LineEdit_pillarwidth")
        self.label_21 = QtWidgets.QLabel(self.centralwidget)
        self.label_21.setGeometry(QtCore.QRect(15, 420, 91, 16))
        self.label_21.setObjectName("label_21")
        self.label_18 = QtWidgets.QLabel(self.centralwidget)
        self.label_18.setGeometry(QtCore.QRect(235, 360, 47, 13))
        self.label_18.setObjectName("label_18")
        self.frame_2 = QtWidgets.QFrame(self.centralwidget)
        self.frame_2.setGeometry(QtCore.QRect(10, 238, 281, 38))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.Box_BB = QtWidgets.QCheckBox(self.frame_2)
        self.Box_BB.setGeometry(QtCore.QRect(110, 0, 111, 17))
        self.Box_BB.setObjectName("Box_BB")
        self.Box_P = QtWidgets.QCheckBox(self.frame_2)
        self.Box_P.setGeometry(QtCore.QRect(110, 20, 111, 17))
        self.Box_P.setObjectName("Box_P")
        self.LineEdit_f_v = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_f_v.setGeometry(QtCore.QRect(115, 208, 113, 20))
        self.LineEdit_f_v.setObjectName("LineEdit_f_v")
        self.label_14 = QtWidgets.QLabel(self.centralwidget)
        self.label_14.setGeometry(QtCore.QRect(15, 50, 71, 16))
        self.label_14.setObjectName("label_14")
        self.LineEdit_postheight = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_postheight.setGeometry(QtCore.QRect(115, 390, 113, 20))
        self.LineEdit_postheight.setObjectName("LineEdit_postheight")
        self.LineEdit_pxsize = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_pxsize.setGeometry(QtCore.QRect(115, 510, 113, 20))
        self.LineEdit_pxsize.setObjectName("LineEdit_pxsize")
        self.LineEdit_coeff_h = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_coeff_h.setGeometry(QtCore.QRect(115, 130, 113, 20))
        self.LineEdit_coeff_h.setObjectName("LineEdit_coeff_h")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(15, 104, 47, 13))
        self.label.setObjectName("label")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(235, 182, 47, 13))
        self.label_8.setObjectName("label_8")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(235, 104, 47, 13))
        self.label_2.setObjectName("label_2")
        self.LineEdit_delta_mat = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_delta_mat.setGeometry(QtCore.QRect(115, 300, 113, 20))
        self.LineEdit_delta_mat.setObjectName("LineEdit_delta_mat")
        self.label_20 = QtWidgets.QLabel(self.centralwidget)
        self.label_20.setGeometry(QtCore.QRect(235, 390, 47, 13))
        self.label_20.setObjectName("label_20")
        self.LineEdit_N_lensletts = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_N_lensletts.setGeometry(QtCore.QRect(115, 72, 113, 20))
        self.LineEdit_N_lensletts.setObjectName("LineEdit_N_lensletts")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(235, 156, 71, 13))
        self.label_6.setObjectName("label_6")
        self.LineEdit_Name = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_Name.setGeometry(QtCore.QRect(115, 50, 113, 20))
        self.LineEdit_Name.setObjectName("LineEdit_Name")
        self.LineEdit_factor_phase = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_factor_phase.setGeometry(QtCore.QRect(115, 480, 113, 20))
        self.LineEdit_factor_phase.setObjectName("LineEdit_factor_phase")
        self.LineEdit_scale_spacer = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_scale_spacer.setGeometry(QtCore.QRect(115, 330, 113, 20))
        self.LineEdit_scale_spacer.setObjectName("LineEdit_scale_spacer")
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(15, 280, 91, 16))
        self.label_13.setObjectName("label_13")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(235, 130, 71, 13))
        self.label_5.setObjectName("label_5")
        self.LineEdit_coeff_v = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_coeff_v.setGeometry(QtCore.QRect(115, 156, 113, 20))
        self.LineEdit_coeff_v.setObjectName("LineEdit_coeff_v")
        self.label_24 = QtWidgets.QLabel(self.centralwidget)
        self.label_24.setGeometry(QtCore.QRect(15, 510, 91, 16))
        self.label_24.setObjectName("label_24")
        self.LineEdit_f_h = QtWidgets.QLineEdit(self.centralwidget)
        self.LineEdit_f_h.setGeometry(QtCore.QRect(115, 182, 113, 20))
        self.LineEdit_f_h.setObjectName("LineEdit_f_h")
        self.label_17 = QtWidgets.QLabel(self.centralwidget)
        self.label_17.setGeometry(QtCore.QRect(15, 360, 91, 16))
        self.label_17.setObjectName("label_17")
        self.label_22 = QtWidgets.QLabel(self.centralwidget)
        self.label_22.setGeometry(QtCore.QRect(15, 450, 91, 16))
        self.label_22.setObjectName("label_22")
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(15, 72, 71, 16))
        self.label_11.setObjectName("label_11")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(235, 208, 47, 13))
        self.label_10.setObjectName("label_10")
        self.Button_SetDir = QtWidgets.QPushButton(self.centralwidget)
        self.Button_SetDir.setGeometry(QtCore.QRect(20, 20, 111, 23))
        self.Button_SetDir.setObjectName("Button_SetDir")
        self.ButtonPreviewMesh = QtWidgets.QPushButton(self.centralwidget)
        self.ButtonPreviewMesh.setGeometry(QtCore.QRect(310, 20, 81, 23))
        self.ButtonPreviewMesh.setCheckable(False)
        self.ButtonPreviewMesh.setFlat(False)
        self.ButtonPreviewMesh.setObjectName("ButtonPreviewMesh")
        self.previewWidget = MPL_WIDGET(self.centralwidget)
        self.previewWidget.setGeometry(QtCore.QRect(310, 50, 471, 491))
        self.previewWidget.setObjectName("previewWidget")
        self.frame.raise_()
        self.label_15.raise_()
        self.RunButton.raise_()
        self.label_9.raise_()
        self.LineEdit_mllsize.raise_()
        self.label_19.raise_()
        self.label_12.raise_()
        self.label_4.raise_()
        self.label_16.raise_()
        self.label_25.raise_()
        self.label_23.raise_()
        self.label_3.raise_()
        self.label_7.raise_()
        self.frame_3.raise_()
        self.LineEdit_pillarwidth.raise_()
        self.label_21.raise_()
        self.label_18.raise_()
        self.frame_2.raise_()
        self.LineEdit_f_v.raise_()
        self.label_14.raise_()
        self.LineEdit_postheight.raise_()
        self.LineEdit_pxsize.raise_()
        self.LineEdit_coeff_h.raise_()
        self.label.raise_()
        self.label_8.raise_()
        self.label_2.raise_()
        self.LineEdit_delta_mat.raise_()
        self.label_20.raise_()
        self.LineEdit_N_lensletts.raise_()
        self.label_6.raise_()
        self.LineEdit_Name.raise_()
        self.LineEdit_factor_phase.raise_()
        self.LineEdit_scale_spacer.raise_()
        self.label_13.raise_()
        self.label_5.raise_()
        self.LineEdit_coeff_v.raise_()
        self.label_24.raise_()
        self.LineEdit_f_h.raise_()
        self.label_17.raise_()
        self.label_22.raise_()
        self.label_11.raise_()
        self.label_10.raise_()
        self.Button_SetDir.raise_()
        self.ButtonPreviewMesh.raise_()
        self.previewWidget.raise_()
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label_15.setText(_translate("MainWindow", "Size spacer"))
        self.RunButton.setText(_translate("MainWindow", "Run"))
        self.label_9.setText(_translate("MainWindow", "focal length (vert)"))
        self.LineEdit_mllsize.setToolTip(_translate("MainWindow", "<html><head/><body><p>Size of the two MLLs in µm. This should be set sufficiently larger than the actual size.</p></body></html>"))
        self.LineEdit_mllsize.setText(_translate("MainWindow", "50"))
        self.label_19.setText(_translate("MainWindow", "Post height"))
        self.label_12.setText(_translate("MainWindow", "Delta (opt const.)"))
        self.label_4.setText(_translate("MainWindow", "Coefficient (vert)"))
        self.label_16.setText(_translate("MainWindow", "µm"))
        self.label_25.setText(_translate("MainWindow", "µm"))
        self.label_23.setText(_translate("MainWindow", "Factor_Phase"))
        self.label_3.setText(_translate("MainWindow", "Coefficient (horz)"))
        self.label_7.setText(_translate("MainWindow", "focal length (horz)"))
        self.LineEdit_pillarwidth.setToolTip(_translate("MainWindow", "<html><head/><body><p>Size of the stabilizing spacer of the fish.</p></body></html>"))
        self.LineEdit_pillarwidth.setText(_translate("MainWindow", "10"))
        self.label_18.setText(_translate("MainWindow", "µm"))
        self.Box_BB.setToolTip(_translate("MainWindow", "<html><head/><body><p>Specify if two stablizing boxes should be added to the model.</p></body></html>"))
        self.Box_BB.setText(_translate("MainWindow", "Bounding Box"))
        self.Box_P.setToolTip(_translate("MainWindow", "<html><head/><body><p>Specify if an underlying post should be added for easier alignment.</p></body></html>"))
        self.Box_P.setText(_translate("MainWindow", "Post"))
        self.LineEdit_f_v.setToolTip(_translate("MainWindow", "<html><head/><body><p>Focal length of the vertically focusing MLL (in µm).</p></body></html>"))
        self.LineEdit_f_v.setText(_translate("MainWindow", "350"))
        self.label_14.setText(_translate("MainWindow", "Name"))
        self.LineEdit_postheight.setToolTip(_translate("MainWindow", "<html><head/><body><p>Size of the stabilizing spacer of the fish.</p></body></html>"))
        self.LineEdit_postheight.setText(_translate("MainWindow", "20"))
        self.LineEdit_pxsize.setToolTip(_translate("MainWindow", "<html><head/><body><p>The sampling size of the model (x,y).</p></body></html>"))
        self.LineEdit_pxsize.setText(_translate("MainWindow", "1"))
        self.LineEdit_coeff_h.setToolTip(_translate("MainWindow", "<html><head/><body><p>Third order coefficient of the wavefront of the horizontally focusing MLL (in rad/mrad³).</p></body></html>"))
        self.LineEdit_coeff_h.setText(_translate("MainWindow", "0.005"))
        self.label.setText(_translate("MainWindow", "Size of MLLs"))
        self.label_8.setText(_translate("MainWindow", "µm"))
        self.label_2.setText(_translate("MainWindow", "µm"))
        self.LineEdit_delta_mat.setToolTip(_translate("MainWindow", "<html><head/><body><p>Decrement of the refractive index of the used material at the given energy.</p></body></html>"))
        self.LineEdit_delta_mat.setText(_translate("MainWindow", "1.72e-6"))
        self.label_20.setText(_translate("MainWindow", "µm"))
        self.LineEdit_N_lensletts.setToolTip(_translate("MainWindow", "<html><head/><body><p>Number of Lensletts. Each lenslett is fish-shaped and therefore has two refracting surfaces</p></body></html>"))
        self.LineEdit_N_lensletts.setText(_translate("MainWindow", "10"))
        self.label_6.setText(_translate("MainWindow", "rad/mrad³"))
        self.LineEdit_Name.setToolTip(_translate("MainWindow", "<html><head/><body><p>Name of the model. Will create Name.stl and Name.log file.</p></body></html>"))
        self.LineEdit_Name.setText(_translate("MainWindow", "CRC"))
        self.LineEdit_factor_phase.setToolTip(_translate("MainWindow", "<html><head/><body><p>An additional factor to scale the magnitude of the correction.</p></body></html>"))
        self.LineEdit_factor_phase.setText(_translate("MainWindow", "1"))
        self.LineEdit_scale_spacer.setToolTip(_translate("MainWindow", "<html><head/><body><p>Size of the stabilizing spacer of the fish (um).</p></body></html>"))
        self.LineEdit_scale_spacer.setText(_translate("MainWindow", "0"))
        self.label_13.setText(_translate("MainWindow", "Advanced:"))
        self.label_5.setText(_translate("MainWindow", "rad/mrad³"))
        self.LineEdit_coeff_v.setToolTip(_translate("MainWindow", "<html><head/><body><p>Third order coefficient of the wavefront of the vertically focusing MLL (in rad/mrad³).</p></body></html>"))
        self.LineEdit_coeff_v.setText(_translate("MainWindow", "0.005"))
        self.label_24.setText(_translate("MainWindow", "pxsize"))
        self.LineEdit_f_h.setToolTip(_translate("MainWindow", "<html><head/><body><p>Focal length of the horizontally focusing MLL (in µm).</p></body></html>"))
        self.LineEdit_f_h.setText(_translate("MainWindow", "250"))
        self.label_17.setText(_translate("MainWindow", "Bounding Box size"))
        self.label_11.setText(_translate("MainWindow", "N Lensletts"))
        self.label_10.setText(_translate("MainWindow", "µm"))
        self.Button_SetDir.setText(_translate("MainWindow", "Set Save Directory"))
        self.ButtonPreviewMesh.setText(_translate("MainWindow", "Preview Mesh"))

