from cutevariant.gui import plugin, FIcon
from cutevariant.core import sql
from cutevariant.core import querybuilder
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *

import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from cutevariant.gui.sql_thread import SqlThread
from functools import partial

matplotlib.use("Qt5Agg")

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureWidget
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationBar
from matplotlib.figure import Figure


class GraphWidget(plugin.PluginWidget):
    """Display all fields according categorie

    Usage:

     view = FieldsWidget
     (conn)
     view.columns = ["chr","pos"]

    """

    LOCATION = plugin.CENTRAL_LOCATION
    ENABLE = True

    def __init__(self, conn=None, parent=None):
        super().__init__(parent)

        self.setWindowTitle(self.tr("Columns"))

        self.stack = QStackedWidget()

        self.label = QLabel("Loading")
        self.label.setAlignment(Qt.AlignCenter)
        self.widget = FigureWidget(Figure(figsize=(5, 5), dpi=100))
        self.widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.bar = NavigationBar(self.widget, self)
        self.ax = self.widget.figure.add_subplot(111)  # create an axis

        self.thread = SqlThread()
        vlayout = QVBoxLayout()

        self.stack.addWidget(self.label)
        self.stack.addWidget(self.widget)

        vlayout.addWidget(self.bar)
        vlayout.addWidget(self.stack)

        self.setLayout(vlayout)

        self.thread.result_ready.connect(self.on_result_received)
        self.thread.error.connect(self.on_error)

    def on_result_received(self):

        self.ax.clear()
        sns_plot = sns.histplot(self.thread.results, ax=self.ax)

        self.widget.draw()
        self.widget.updateGeometry()
        self.stack.setCurrentIndex(1)

    def on_error(self, message):
        print("error", message)
        self.stack.setCurrentIndex(1)

    def on_open_project(self, conn):
        """ Overrided from PluginWidget """
        self.conn = conn
        self.thread.conn = conn
        # self.model.conn = conn
        # self.model.load()
        # self.on_refresh()

    def on_refresh(self):
        """ overrided from PluginWidget """

        q = querybuilder.build_full_sql_query(
            self.conn,
            fields=["mq", "qual"],
            source=self.mainwindow.state.source,
            filters=self.mainwindow.state.filters,
            limit=None,
        )

        q = q.replace("DISTINCT", "")  # TODO ... params !

        func = functools.partial(
            lambda conn: [i["qual"] for i in conn.execute(q).fetchall()]
        )

        self.thread.start_function(func)
        self.stack.setCurrentIndex(0)
