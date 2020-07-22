from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import QSizePolicy,QWidget,QVBoxLayout
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


class MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig=Figure()
        #self.ax=self.fig.add_subplot(111)
        self.ax=mplot3d.Axes3D(self.fig)
        FigureCanvas.__init__(self,self.fig)
        FigureCanvas.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

class MPL_WIDGET(QWidget):
    def __init__(self,parent=None):
        QWidget.__init__(self,parent)
        self.canvas=MplCanvas()
        #self.navi_toolbar=NavigationToolbar(self.canvas,self)
        self.vbl=QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        #self.vbl.addWidget(self.navi_toolbar)
        self.setLayout(self.vbl)


