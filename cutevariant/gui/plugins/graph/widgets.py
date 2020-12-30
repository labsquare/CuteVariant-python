from cutevariant.gui import plugin, FIcon
from cutevariant.core import sql
from cutevariant.core import querybuilder
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *

import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt

matplotlib.use("Qt5Agg")

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureWidget
from matplotlib.figure import Figure


class GraphWidget(plugin.PluginWidget):
    """Display all fields according categorie

    Usage:

     view = FieldsWidget
     (conn)
     view.columns = ["chr","pos"]

    """

    ENABLE = True

    def __init__(self, conn=None, parent=None):
        super().__init__(parent)

        self.setWindowTitle(self.tr("Columns"))

        self.widget = FigureWidget(Figure(figsize=(5, 5), dpi=100))
        self.widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.ax = self.widget.figure.add_subplot(111)  # create an axis

        vlayout = QVBoxLayout()
        vlayout.addWidget(self.widget)

        self.setLayout(vlayout)

    def on_open_project(self, conn):
        """ Overrided from PluginWidget """
        self.conn = conn
        # self.model.conn = conn
        # self.model.load()
        # self.on_refresh()

    def on_refresh(self):
        """ overrided from PluginWidget """

        q = querybuilder.build_sql_query(
            fields=["dp"],
            source=self.mainwindow.state.source,
            filters=self.mainwindow.state.filters,
        )

        DP = [i["dp"] for i in self.conn.execute(q).fetchall()]

        print(DP)
        # self.ax.plot(range(5), range(5))
        self.ax.clear()
        sns_plot = sns.histplot(DP, ax=self.ax, kde=True)

        self.widget.draw()
        self.widget.updateGeometry()
