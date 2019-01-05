from PySide2.QtCore import *
from PySide2.QtWidgets import *
import sqlite3

from cutevariant.gui.viewquerywidget import ViewQueryWidget
from cutevariant.gui.columnquerywidget import ColumnQueryWidget
from cutevariant.gui.queryrouter import QueryRouter

from cutevariant.core.importer import import_file
from cutevariant.core import Query


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__()
        self.toolbar = self.addToolBar("test")
        self.conn   = None
        self.view_widget   = ViewQueryWidget()
        self.column_widget = ColumnQueryWidget()
        
        # Init router 
        self.router = QueryRouter()
        self.router.addWidget(self.view_widget)
        self.router.addWidget(self.column_widget)

        # Init panel 
        self.addPanel(self.column_widget)


        self.setCentralWidget(self.view_widget)


        self.open("/tmp/qt_cutevariant.db")

    def open(self, filename):

        import_file("exemples/test.csv", filename)
        self.conn = sqlite3.connect(filename)
        self.router.setQuery(Query(self.conn))
        


    def addPanel(self, widget, area = Qt.LeftDockWidgetArea):
        dock = QDockWidget()
        dock.setWindowTitle(widget.windowTitle())
        dock.setWidget(widget)
        self.addDockWidget(area, dock)




