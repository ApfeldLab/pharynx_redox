import matplotlib.axis as maxis
import numpy as np
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.projections import register_projection
from matplotlib.ticker import NullLocator
from scipy.interpolate import UnivariateSpline


class MidlineSplineEditorPlot(object):
    show_vertices = True
    epsilon = 5  # max pixel distance to count as a vertex hit

    def __init__(self, canvas, ax, splines, spline_xs, n_spline_control_points=5, frame=0):
        """
        x = vector of x positions of midline
        y = vector of y positions of midline
        """
        self.frame = frame
        self.canvas = canvas
        self.ax = ax
        self.splines = splines
        self.spline_xs = spline_xs
        self.current_spline = splines[self.frame]
        self.spl_ctrl_xs = np.linspace(np.min(spline_xs), np.max(spline_xs), n_spline_control_points)
        self.spl_ctrl_ys = self.splines[self.frame](self.spl_ctrl_xs)
        self.spline_resolution = 100

        self.pts = Line2D(self.spl_ctrl_xs, self.splines[self.frame](self.spl_ctrl_xs),
                          zorder=100, color='white', ls='', marker='o', animated=True)
        self.line = Line2D(self.spline_xs, self.splines[self.frame](self.spline_xs), color='orange', animated=True,
                           zorder=99)

        self.active_vertex_idx = None

        self.ax.add_artist(self.line)
        self.ax.add_artist(self.pts)

        self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)

        self.background = self.canvas.copy_from_bbox(self.ax.bbox)

        self.canvas.draw()
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.pts)

    def draw_callback(self, _):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.pts)
        # do not need to blit here, this will fire before the screen is updated

    def get_ind_under_point(self, event):
        d = np.hypot(self.spl_ctrl_xs - event.xdata, self.spl_ctrl_ys - event.ydata)
        index_sequence, = np.nonzero(d == d.min())
        ind = index_sequence[0]
        if d[ind] >= self.epsilon:
            ind = None
        return ind

    def button_press_callback(self, event):
        """Called when a mouse button is pressed"""
        if event.inaxes is None:
            return
        if event.button != 1:  # check that it is left mouse
            return
        self.active_vertex_idx = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        """Called when a mouse button is released"""
        if event.button != 1:
            return
        self.active_vertex_idx = None

    def make_spline(self, xs, ys):
        return UnivariateSpline(xs, ys, ext=0, s=0, k=3)

    def motion_notify_callback(self, event):
        if self.active_vertex_idx is None:
            return
        if event.inaxes is None:
            return
        mouse_x, mouse_y = event.xdata, event.ydata

        old_x = self.spl_ctrl_xs
        old_y = self.spl_ctrl_ys
        try:
            self.spl_ctrl_xs[self.active_vertex_idx] = mouse_x
            self.spl_ctrl_ys[self.active_vertex_idx] = mouse_y

            self.splines[self.frame] = self.make_spline(self.spline_xs, self.splines[self.frame](self.spl_ctrl_xs))
            self.spline_xs = np.linspace(np.min(self.spline_xs), np.max(self.spline_xs), self.spline_resolution)
            self.line.set_data(self.spline_xs, self.splines[self.frame](self.spline_xs))
            self.pts.set_data(self.spl_ctrl_xs, self.spl_ctrl_ys)
        except ValueError:
            self.spl_ctrl_xs = old_x
            self.spl_ctrl_ys = old_y

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.pts)
        self.canvas.blit(self.ax.bbox)


class NoTicksXAxis(maxis.XAxis):
    def reset_ticks(self):
        self._lastNumMajorTicks = 1
        self._lastNumMinorTicks = 1

    def set_clip_path(self, clippath, transform=None):
        pass


class NoTicksYAxis(maxis.YAxis):
    def reset_ticks(self):
        self._lastNumMajorTicks = 1
        self._lastNumMinorTicks = 1

    def set_clip_path(self, clippath, transform=None):
        pass


class ThinAxes(Axes):
    """Thin axes without spines and ticks to accelerate axes creation"""

    name = 'thin'

    def _init_axis(self):
        self.xaxis = NoTicksXAxis(self)
        self.yaxis = NoTicksYAxis(self)

    def cla(self):
        """
        Override to set up some reasonable defaults.
        """
        Axes.cla(self)
        self.xaxis.set_minor_locator(NullLocator())
        self.yaxis.set_minor_locator(NullLocator())
        self.xaxis.set_major_locator(NullLocator())
        self.yaxis.set_major_locator(NullLocator())

    def _gen_axes_spines(self):
        return {}


# Now register the projection with matplotlib so the user can select
# it.
register_projection(ThinAxes)
