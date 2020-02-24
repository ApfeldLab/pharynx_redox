import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import pyqtgraph as pg
from pyqtgraph import Point
from pyqtgraph.graphicsItems.ROI import ROI
from scipy import interpolate
from scipy.spatial.distance import cdist
from dataclasses import dataclass, astuple
from pharynx_redox.gui.qt_py_files.spline_editor import Ui_Form
import numpy as np


class Spline:
    def __init__(self, ctrl_pts=[], n_points=500):
        """
        Initialize the Spline object
        
        Parameters
        ----------
        ctrl_pts : list, optional
            the list of points that this spline will interpolate between, by default []
        n_points : int, optional
            the number of points under which to evaluate the spline, by default 100
        """
        self._spl = None
        self.ctrl_pts = []

        self.add_ctrl_pts(ctrl_pts)

    def add_ctrl_pts(self, ctrl_pts):
        for (x, y) in ctrl_pts:
            self.ctrl_pts.append(Point(x, y))
        self._update_spline()

    def set_ctrl_pts(self, pos):
        self.ctrl_pts = []
        self.add_ctrl_pts(pos)

    def _update_spline(self):
        xys = np.asarray(list(map(list, self.ctrl_pts)))
        if len(xys) > 0:
            self._spl = interpolate.Akima1DInterpolator(
                np.linspace(0, 1, len(xys)), xys
            )

    def __call__(self, n=500):
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
        xs = np.linspace(0, 1, n)
        return self._spl(xs)


class SplineROI(pg.GraphicsObject):

    sigClicked = QtCore.Signal(object)

    def __init__(self, pos, *args, pen="default", **kwargs):
        super(SplineROI, self).__init__(*args, **kwargs)

        self.dragPoint = None
        self.dragOffset = None

        self.ctrl_pts = pos

        self.spl = Spline(ctrl_pts=pos)

        self.scatter = pg.ScatterPlotItem(pxMode=True)
        self.scatter.setParentItem(self)
        self.scatter.sigClicked.connect(self.clicked)

        self.picture = None
        self.pen = pen
        self.setData(ctrl_pts=pos)

        self._mouseShape = None

    def _update(self):
        self.picture = None
        self.prepareGeometryChange()
        self.mouseShape()
        self.update()

    def setData(self, ctrl_pts):
        self.ctrl_pts = ctrl_pts
        self.scatter.setData(pos=ctrl_pts)
        self.spl.set_ctrl_pts(pos=ctrl_pts)
        self.informViewBoundsChanged()
        self._update()

    def mouseDoubleClickEvent(self, ev):
        if self.mouseShape().contains(ev.pos()):
            ev.accept()
            self.sigClicked.emit(self)

            click_pos = np.asarray([ev.pos().x(), ev.pos().y()])
            crv_data = self.spl(n=1000)

            ctrl_pts_passed = 0
            epsilon = 1

            # make a copy so we can mutate it without losing ctrl pts in the app
            _ctrl_pts = self.ctrl_pts.copy()

            # We find the closest point on the curve to where we clicked
            closest_crv_pt_idx = np.argmin(cdist([click_pos,], crv_data))

            # and move along the curve until we get to that point
            # all the while, we are interested in how many control points we pass
            for crv_pt in crv_data[:closest_crv_pt_idx, :]:
                ctrl_pt = _ctrl_pts[0]
                if cdist([crv_pt,], [ctrl_pt])[0] <= epsilon:
                    # we have reached a control point
                    ctrl_pts_passed += 1

                    # remove that ctrl point from future considerations (because there
                    # are likely many curve points that are close to the curve point we
                    # just passed)
                    _ctrl_pts.pop(0)

            self.ctrl_pts.insert(ctrl_pts_passed, tuple(click_pos))
            self.setData(self.ctrl_pts)

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
            self.dragOffset = self.ctrl_pts[ind] - pos
        elif ev.isFinish():
            self.dragPoint = None
            return
        else:
            if self.dragPoint is None:
                ev.ignore()
                return

        ind = self.dragPoint._index
        self.ctrl_pts[ind] = tuple(ev.pos() + self.dragOffset)
        self.setData(ctrl_pts=self.ctrl_pts)
        ev.accept()

    def remove_point(self, idx):
        self.ctrl_pts.pop(idx)
        self.setData(self.ctrl_pts)

    def clicked(self, pts):
        modifiers = QtWidgets.QApplication.keyboardModifiers()
        if modifiers == QtCore.Qt.ControlModifier:
            idx_clicked = pts.ptsClicked[0]._index
            if len(self.ctrl_pts) > 2:
                self.remove_point(idx_clicked)

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
            p.setPen(pg.functions.mkPen(color="w", width=3))
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
        # return pg.PlotCurveItem(pos=self.pos).dataBounds(*args, **kwds)
        return self.scatter.dataBounds(*args, **kwds)

    def pixelPadding(self):
        return self.scatter.pixelPadding()

    def getPath(self):
        crv_data = self.spl(n=1000)
        x, y = crv_data[:, 0], crv_data[:, 1]
        return pg.functions.arrayToQPath(x, y)

    def mouseShape(self):
        """
        Return a QPainterPath representing the clickable shape of the curve

        """
        # if self._mouseShape is None:
        view = self.getViewBox()
        if view is None:
            return QtGui.QPainterPath()
        stroker = QtGui.QPainterPathStroker()
        path = self.getPath()
        path = self.mapToItem(view, path)
        stroker.setWidth(10)
        mousePath = stroker.createStroke(path)
        self._mouseShape = self.mapFromItem(view, mousePath)
        return self._mouseShape


class SplineEditorTestWidget(QtGui.QWidget):
    """Draw and edit a spline"""

    def __init__(self):
        super(SplineEditorTestWidget, self).__init__()

        self.spl_roi = SplineROI(pos=[(0, 0), (1, 1), (10, 10), (10, 20)])

        self.setGeometry(300, 300, 450, 450)
        self.setWindowTitle("Bezier Curves")

        ###############################
        # Set up UI
        ###############################
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.vb = self.ui.canvas.addViewBox()
        self.vb.addItem(self.spl_roi)
        self.vb.addItem(pg.GridItem())


if __name__ == "__main__":
    qapp = QtWidgets.QApplication(sys.argv)
    se = SplineEditorTestWidget()
    se.show()
    qapp.exec_()
