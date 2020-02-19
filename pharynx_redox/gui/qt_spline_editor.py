import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import pyqtgraph as pg
from pyqtgraph import Point
from pyqtgraph.graphicsItems.ROI import ROI
from scipy import interpolate
from dataclasses import dataclass, astuple
from pharynx_redox.gui.qt_py_files.spline_editor import Ui_Form
import numpy as np


class Spline:
    def __init__(self, ctrl_pts=[], n_points=100):
        """
        Initialize the Spline object
        
        Parameters
        ----------
        ctrl_pts : list, optional
            a , by default []
        n_points : int, optional
            the number of points under which to evaluate the spline, by default 100
        """
        self._spl = None

        self._xs = np.linspace(0, 1, n_points)
        self.ctrl_pts = []

        self.add_ctrl_pts(ctrl_pts)
        self.update_spline()

    def add_ctrl_pts(self, ctrl_pts):
        for (x, y) in ctrl_pts:
            self.ctrl_pts.append(Point(x, y))

    def update_spline(self):
        xys = np.asarray(list(map(list, self.ctrl_pts)))
        if len(xys) > 0:
            self._spl = interpolate.CubicSpline(np.linspace(0, 1, len(xys)), xys)

    def set_ctrl_pts(self, pos):
        self.ctrl_pts = pos
        self.update_spline()

    def __call__(self, xs=None):
        """
        evaluate the spline
        
        Parameters
        ----------
        xs : list of numbers, optional
            the list of numbers (or numpy array) to evaluate the spline at. must be [0,1]
        
        Returns
        -------
        np.ndarray
            a Mx2 numpy array where 
                M[:, 0] = x
                M[:, 1] = y
        """
        if xs is None:
            return self._spl(self._xs)
        else:
            return self._spl(xs)


class SplineROI(pg.GraphItem):
    def __init__(self, *args, **kwargs):
        self.dragPoint = None
        self.dragOffset = None
        super(SplineROI, self).__init__(*args, **kwargs)

    def setData(self, **kwds):
        self.data = kwds
        if "pos" in self.data:
            npts = self.data["pos"].shape[0]
            self.data["adj"] = np.column_stack(
                (np.arange(0, npts - 1), np.arange(1, npts))
            )
            self.data["data"] = np.empty(npts, dtype=[("index", int)])
            self.data["data"]["index"] = np.arange(npts)
        self.updateGraph()

    def updateGraph(self):
        pg.GraphItem.setData(self, **self.data)

    def mouseDragEvent(self, ev):
        if ev.button() != QtCore.Qt.LeftButton:
            ev.ignore()
            return

        if ev.isStart():
            pos = ev.buttonDownPos()
            pts = self.scatter.pointsAt(pos)
            if len(pts) == 0:
                ev.ignore()
                return
            self.dragPoint = pts[0]
            ind = pts[0].data()[0]
            self.dragOffset = self.data["pos"][ind][1] - pos[1]
        elif ev.isFinish():
            self.dragPoint = None
            return
        else:
            if self.dragPoint is None:
                ev.ignore()
                return

        ind = self.dragPoint.data()[0]
        self.data["pos"][ind][1] = ev.pos()[1] + self.dragOffset
        self.updateGraph()
        ev.accept()


class Graph(pg.GraphItem):
    def __init__(self):
        self.dragPoint = None
        self.dragOffset = None
        self.textItems = []
        self.crv = pg.PlotCurveItem()
        pg.GraphItem.__init__(self)
        self.scatter.sigClicked.connect(self.clicked)

    def setData(self, **kwds):
        self.text = kwds.pop("text", [])
        self.data = kwds
        if "pos" in self.data:
            npts = self.data["pos"].shape[0]
            self.data["data"] = np.empty(npts, dtype=[("index", int)])
            self.data["data"]["index"] = np.arange(npts)
        self.setTexts(self.text)
        self.updateGraph()

    def setTexts(self, text):
        for i in self.textItems:
            i.scene().removeItem(i)
        self.textItems = []
        for t in text:
            item = pg.TextItem(t)
            self.textItems.append(item)
            item.setParentItem(self)

    def updateGraph(self):
        if "pos" in self.data:
            print(self.data["pos"])
            crv_data = Spline(ctrl_pts=self.data["pos"])()
            self.crv.setData(pos=crv_data)
        pg.GraphItem.setData(self, **self.data)
        for i, item in enumerate(self.textItems):
            item.setPos(*self.data["pos"][i])

    def mouseDragEvent(self, ev):
        if ev.button() != QtCore.Qt.LeftButton:
            ev.ignore()
            return

        if ev.isStart():
            # We are already one step into the drag.
            # Find the point(s) at the mouse cursor when the button was first
            # pressed:
            pos = ev.buttonDownPos()
            pts = self.scatter.pointsAt(pos)
            if len(pts) == 0:
                ev.ignore()
                return
            self.dragPoint = pts[0]
            ind = pts[0].data()[0]
            self.dragOffset = self.data["pos"][ind] - pos
        elif ev.isFinish():
            self.dragPoint = None
            return
        else:
            if self.dragPoint is None:
                ev.ignore()
                return

        ind = self.dragPoint.data()[0]
        self.data["pos"][ind] = ev.pos() + self.dragOffset
        self.setData(pos=self.pos)
        ev.accept()

    def clicked(self, pts):
        print("clicked: %s" % pts)

    def paint(self, p, *args):
        p.setPen(self.currentPen)

        # super().paint(p, *args)


class OldSplineROI(pg.GraphicsObject):
    def __init__(self, pos, *args, pen="default", **kwargs):
        super(OldSplineROI, self).__init__(*args, **kwargs)

        self.dragPoint = None
        self.dragOffset = None

        self.spl = Spline(ctrl_pts=pos)

        self.scatter = pg.ScatterPlotItem(pxMode=True)
        self.scatter.setParentItem(self)
        self.scatter.sigClicked.connect(self.clicked)

        self.picture = None
        self.pen = pen
        self.setData(pos=pos)

    def _update(self):
        self.picture = None
        self.prepareGeometryChange()
        self.update()

    def setData(self, **kwds):
        self.data = kwds
        # if "pos" in self.data:
        # npts = len(self.data["pos"])
        # self.data["index"] = np.arange(npts)

        self._update()
        self.scatter.setData(**kwds)
        self.spl.set_ctrl_pts(pos=self.data["pos"])
        self.informViewBoundsChanged()
        self._update()

    def mouseDragEvent(self, ev):
        if ev.button() != QtCore.Qt.LeftButton:
            ev.ignore()
            return

        if ev.isStart():
            # We are already one step into the drag.
            # Find the point(s) at the mouse cursor when the button was first
            # pressed:
            pos = ev.buttonDownPos()
            pts = self.scatter.pointsAt(pos)
            if len(pts) == 0:
                ev.ignore()
                return
            self.dragPoint = pts[0]
            ind = pts[0]._index
            self.dragOffset = self.data["pos"][ind] - pos
        elif ev.isFinish():
            self.dragPoint = None
            return
        else:
            if self.dragPoint is None:
                ev.ignore()
                return

        # ind = self.dragPoint.data()[0]
        ind = self.dragPoint._index
        self.data["pos"][ind] = ev.pos() + self.dragOffset
        self.setData(pos=self.data["pos"])
        ev.accept()

    def clicked(self, pts):
        print("clicked: %s" % pts)

    def generatePicture(self):
        self.picture = QtGui.QPicture()
        if self.pen is None or self.pos is None:
            return

        p = QtGui.QPainter(self.picture)
        try:
            xs = self.spl()[:, 0]
            ys = self.spl()[:, 1]
            pen = self.pen
            if pen == "default":
                pen = pg.getConfigOption("foreground")
            p.setPen(pg.functions.mkPen(pen))
            path = pg.functions.arrayToQPath(x=xs, y=ys, connect="all")
            p.drawPath(path)
        finally:
            p.end()

    def paint(self, p, *args, **kwargs):
        if self.picture == None:
            self.generatePicture()
        if pg.getConfigOption("antialias") is True:
            print("setting antialias")
        p.setRenderHint(True)
        self.picture.play(p)

    def boundingRect(self):
        xmn, xmx = float(np.min(self.spl()[:, 0])), float(np.max(self.spl()[:, 0]))
        ymn, ymx = float(np.min(self.spl()[:, 1])), float(np.max(self.spl()[:, 1]))

        px = py = 0.0
        pxPad = self.pixelPadding()
        if pxPad > 0:
            # determine length of pixel in local x, y directions
            px, py = self.pixelVectors()
            try:
                px = 0 if px is None else px.length()
            except OverflowError:
                px = 0
            try:
                py = 0 if py is None else py.length()
            except OverflowError:
                py = 0

            # return bounds expanded by pixel size
            px *= pxPad
            py *= pxPad

        return QtCore.QRectF(
            xmn - px, ymn - py, (2 * px) + xmx - xmn, (2 * py) + ymx - ymn
        )

    def dataBounds(self, *args, **kwds):
        # return pg.PlotCurveItem(pos=self.data["pos"]).dataBounds(*args, **kwds)
        return self.scatter.dataBounds(*args, **kwds)

    def pixelPadding(self):
        return self.scatter.pixelPadding()


class SplineEditorWidget(QtGui.QWidget):
    """Draw and edit a spline"""

    def __init__(self):
        super(SplineEditorWidget, self).__init__()

        # self.spl = Spline()
        self.spl_roi = OldSplineROI(pos=[(10, 10), (0, 0), (1, 2)])
        # self.spl_roi = Graph()
        # self.spl_roi.setData(pos=np.asarray([(10, 10), (0, 0), (1, 2)]))

        self.setGeometry(300, 300, 450, 450)
        self.setWindowTitle("Bezier Curves")

        ###############################
        # Set up UI
        ###############################
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.vb = self.ui.canvas.addViewBox()
        self.vb.addItem(self.spl_roi)

        ################################
        # Set up State Machine
        ################################
        self.state_machine = QtCore.QStateMachine()

        # Add states
        ############
        self.adding_ctrl_pts_state = QtCore.QState()
        self.not_editing_state = QtCore.QState()
        self.removing_ctrl_pts_state = QtCore.QState()

        # Enter/Exit State functions
        self.adding_ctrl_pts_state.entered.connect(self.enter_adding_pts_state)
        self.removing_ctrl_pts_state.entered.connect(self.enter_removing_pts_state)
        self.not_editing_state.entered.connect(self.enter_default_state)

        # Transitions
        self.not_editing_state.addTransition(
            self.ui.pushButton.pressed, self.adding_ctrl_pts_state,
        )

        self.adding_ctrl_pts_state.addTransition(
            self.ui.pushButton.pressed, self.removing_ctrl_pts_state
        )
        self.removing_ctrl_pts_state.addTransition(
            self.ui.pushButton.pressed, self.not_editing_state
        )

        self.state_machine.addState(self.adding_ctrl_pts_state)
        self.state_machine.addState(self.removing_ctrl_pts_state)
        self.state_machine.addState(self.not_editing_state)

        self.state_machine.setInitialState(self.not_editing_state)

        self.state_machine.start()

    def enter_adding_pts_state(self):
        self.ui.pushButton.setText("remove points")

    def enter_removing_pts_state(self):
        self.ui.pushButton.setText("exit point editor")

    def enter_default_state(self):
        self.ui.pushButton.setText("add points")


if __name__ == "__main__":
    qapp = QtWidgets.QApplication(sys.argv)
    # s = SplineROI([(0, 1), (1, 2)])
    se = SplineEditorWidget()
    se.show()
    qapp.exec_()
