# Standard imports
import typing
import glob
import json
import os
import gzip
import sys
import sqlite3
import copy

# Qt imports
from PySide2.QtCore import *
from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtSql import *

# Custom imports
from cutevariant.gui import style, plugin, FIcon, MainWindow
from cutevariant.core.querybuilder import build_vql_query, build_sql_query
from cutevariant.commons import logger

from cutevariant.core import sql

from cutevariant.gui.widgets import VqlSyntaxHighlighter


LOGGER = logger()


def refGene_to_sqlite(ref_filename: str, db_filename: str):

    # Create databases
    conn = sqlite3.connect(db_filename)

    conn.execute(
        """
        CREATE TABLE refGene(  
        id INTEGER PRIMARY KEY,
        transcript TEXT, 
        txStart INTEGER, 
        txEnd INTEGER,
        cdsStart INTEGER,
        cdsEnd INTEGER,
        exonStarts TEXT,
        exonEnds TEXT,
        gene TEXT
        )

    """
    )

    data = []
    with gzip.open(ref_filename, "rb") as file:
        for index, line in enumerate(file):
            if line:
                line = line.decode("utf-8").strip().split("\t")

                transcript = line[1]
                txStart = line[4]
                txEnd = line[5]
                cdsStart = line[6]
                cdsEnd = line[7]
                exonStarts = line[9]
                exonEnds = line[10]
                gene = line[12]

                data.append(
                    (
                        None,
                        transcript,
                        txStart,
                        txEnd,
                        cdsStart,
                        cdsEnd,
                        exonStarts,
                        exonEnds,
                        gene,
                    )
                )

    conn.executemany("INSERT INTO refGene VALUES(?,?,?,?,?,?,?,?,?);", data)
    conn.commit()


class GeneView(QAbstractScrollArea):
    """docstring for ClassName"""

    def __init__(self, parent=None):
        super().__init__(parent)

        self.name = "NM_000492"
        self.chrom = "chr7"
        self.strand = "+"
        self.tx_start = 117120078
        self.tx_end = 117308718

        self.cds_start = 117120148
        self.cds_end = 117307162

        self.exon_starts = [
            117120078,
            117144306,
            117149087,
            117170952,
            117174329,
            117175301,
            117176601,
            117180153,
            117182069,
            117188694,
            117199517,
            117227792,
            117230406,
            117231987,
            117234983,
            117242879,
            117243585,
            117246727,
            117250572,
            117251634,
            117254666,
            117267575,
            117282491,
            117292895,
            117304741,
            117305512,
            117306961,
        ]
        self.exon_ends = [
            117120201,
            117144417,
            117149196,
            117171168,
            117174419,
            117175465,
            117176727,
            117180400,
            117182162,
            117188877,
            117199709,
            117227887,
            117230493,
            117232711,
            117235112,
            117242917,
            117243836,
            117246807,
            117250723,
            117251862,
            117254767,
            117267824,
            117282647,
            117292985,
            117304914,
            117305618,
            117308718,
        ]

        self.variants = [(117227892, "red"), (117243866, "red")]

        self.exon_count = len(self.exon_starts)

        # style
        self.exon_height = 30
        self.intron_height = 20

        # self.showMaximized()

        self._area = None

        self.scale_factor = 1
        self.translation = 0

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.horizontalScrollBar().setRange(0, 0)

        # QScroller.grabGesture(self, QScroller.LeftMouseButtonGesture);

        self.horizontalScrollBar().valueChanged.connect(self.set_translation)

        self.resize(640, 200)
        QScroller.grabGesture(self.viewport(), QScroller.LeftMouseButtonGesture)

    def paintEvent(self, event: QPaintEvent):

        painter = QPainter()
        painter.begin(self.viewport())
        painter.setBrush(QBrush(QColor("white")))
        painter.drawRect(self.rect())

        # draw guide
        area = self.draw_area()
        painter.drawRect(area.adjusted(-2, -2, 2, 2))

        self.marks = []

        # Draw variants

        painter.save()
        for variant in self.variants:

            pos, color = variant

            x = self._pixel_to_scroll(self._dna_to_pixel(pos)) + area.left()
            pen = QPen()
            pen.setColor(QColor(color))
            pen.setCapStyle(Qt.RoundCap)
            pen.setWidth(10)
            painter.setPen(pen)
            mark = QPoint(x, self.viewport().height() / 4)
            base = QPoint(x, self.viewport().height() / 2)

            painter.drawPoint(mark)

            pen.setWidth(1)
            painter.setPen(pen)
            painter.drawLine(mark, base)

        painter.restore()

        # draw rule
        painter.drawLine(
            area.left(), area.center().y(), area.right(), area.center().y()
        )

        # Draw intron Background
        intron_rect = QRect(area)
        intron_rect.setHeight(self.intron_height)
        intron_rect.moveCenter(QPoint(area.center().x(), area.center().y()))
        linearGrad = QLinearGradient(
            QPoint(0, intron_rect.top()), QPoint(0, intron_rect.bottom())
        )
        linearGrad.setColorAt(0, QColor("#FDFD97"))
        linearGrad.setColorAt(1, QColor("#FDFD97").darker())
        brush = QBrush(linearGrad)
        painter.setBrush(brush)
        painter.drawRect(intron_rect)

        # draw exons

        painter.setClipRect(area)

        charles_cds_start = (
            self._pixel_to_scroll(self._dna_to_pixel(self.cds_start)) + area.left()
        )
        charles_cds_end = (
            self._pixel_to_scroll(self._dna_to_pixel(self.cds_end)) + area.left()
        )

        rect = QRect()
        rect.setLeft(charles_cds_start)
        rect.setRight(charles_cds_end)
        rect.setHeight(20)
        painter.drawRect(rect)

        for i in range(self.exon_count):

            start = self._dna_to_pixel(self.exon_starts[i])
            end = self._dna_to_pixel(self.exon_ends[i])

            start = self._pixel_to_scroll(start)
            end = self._pixel_to_scroll(end)

            for m in self.marks:
                print("mark", m)
                painter.drawLine(m, 0, m, area.bottom())

            # draw exons
            exon_rect = QRect(0, 0, end - start, self.exon_height)
            exon_rect.moveTo(
                start + area.left(), area.center().y() - self.exon_height / 2
            )

            painter.drawText(exon_rect, Qt.AlignCenter, str(i))
            linearGrad = QLinearGradient(
                QPoint(0, exon_rect.top()), QPoint(0, exon_rect.bottom())
            )
            linearGrad.setColorAt(0, QColor("#789FCC"))
            linearGrad.setColorAt(1, QColor("#789FCC").darker())
            brush = QBrush(linearGrad)
            painter.setBrush(brush)
            painter.drawRect(exon_rect)

            # Draw CDS

        painter.end()

    def _pixel_to_dna(self, pixel: int):

        tx_size = self.tx_end - self.tx_start
        scale = tx_size / self.draw_area().width()
        return pixel * scale + self.tx_start

    def _dna_to_pixel(self, dna: int):

        # normalize dna
        dna = dna - self.tx_start
        tx_size = self.tx_end - self.tx_start
        scale = self.draw_area().width() / tx_size
        return dna * scale

    def _pixel_to_scroll(self, pixel):
        return pixel * self.scale_factor - self.translation

    def _scroll_to_pixel(self, pos):
        return (pos + self.translation) / self.scale_factor

    @Slot(int)
    def set_scale(self, x):

        self.scale_factor = x

        min_scroll = 0
        max_scroll = (
            self.draw_area().width() * self.scale_factor
        ) - self.draw_area().width()

        previous = self.horizontalScrollBar().value()
        previous_max = self.horizontalScrollBar().maximum()

        self.horizontalScrollBar().setRange(min_scroll, max_scroll)

        if previous_max > 1:
            new = previous * self.horizontalScrollBar().maximum() / previous_max
        else:
            new = self.horizontalScrollBar().maximum() / 2

        self.horizontalScrollBar().setValue(new)

    def set_translation(self, x):
        self.translation = x
        self.viewport().update()

    def update_scroll(self):
        pass

    def draw_area(self):
        return self.viewport().rect().adjusted(10, 10, -10, -10)

    def wheelEvent(self, event):

        if event.delta() > 0:
            self.set_scale(self.scale_factor + 0.5)
        else:
            if self.scale_factor > 1:
                self.set_scale(self.scale_factor - 0.5)

    def mousePressEvent(self, event):

        super().mousePressEvent(event)
        # print("mark")
        # self.marks.append(self.horizontalScrollBar().value())
        # self.viewport().update()


class GeneViewerWidget(plugin.PluginWidget):
    """Exposed class to manage VQL/SQL queries from the mainwindow"""

    ENABLE = True
    REFRESH_ONLY_VISIBLE = False

    def __init__(self, parent=None):
        super().__init__(parent)

        self.filename = "/home/sacha/refGene.db"

        self.setWindowTitle(self.tr("Gene Viewer"))

        self.view = GeneView()
        self.combo = QComboBox()
        self.combo.setEditable(True)
        self.toolbar = QToolBar()
        self.toolbar.addWidget(self.combo)

        self.vlayout = QVBoxLayout(self)
        self.vlayout.addWidget(self.toolbar)
        self.vlayout.addWidget(self.view)

        self.model = QSqlQueryModel()
        self.combo.setModel(self.model)
        self.load_combo()

        self.combo.currentIndexChanged.connect(self.on_combo_changed)

    def on_open_project(self, conn):
        self.conn = conn

    def on_refresh(self):

        gene = self.combo.currentText()

        filters = copy.deepcopy(self.mainwindow.state.filters)
        fields = ["pos"]

        filters["$and"].append({"ann.gene": gene})

        variants = sql.get_variants(
            self.conn,
            fields,
            self.mainwindow.state.source,
            filters,
        )

        self.view.variants = [(variant["pos"], "red") for variant in variants]

        self.view.viewport().update()

    def load_combo(self):
        db = QSqlDatabase("QSQLITE")
        db.setDatabaseName(self.filename)
        if db.open():
            self.model.setQuery("SELECT gene FROM refGene", db)

    def on_combo_changed(self):

        gene = self.combo.currentText()

        conn = sqlite3.connect(self.filename)

        results = next(conn.execute(f"SELECT * FROM refGene WHERE gene = '{gene}'"))

        (
            _,
            transcript,
            txStart,
            txEnd,
            cdsStart,
            cdsEnd,
            exonsStarts,
            exonEnds,
            gene,
        ) = results

        self.view.tx_start = txStart
        self.view.tx_end = txEnd

        self.view.cds_start = cdsStart
        self.view.cds_end = cdsEnd

        self.view.exon_starts = [int(i) for i in exonsStarts.split(",")[:-1]]
        self.view.exon_ends = [int(i) for i in exonEnds.split(",")[:-1]]

        print(self.view.exon_starts)

        self.view.exon_count = len(self.view.exon_starts)

        self.view.viewport().update()


if __name__ == "__main__":

    pass

    # import sys
    # import sqlite3

    # import os

    # try:
    #     os.remove("/home/sacha/refGene.db")
    # except:
    #     pass

    # refGene_to_sqlite("/home/sacha/refGene.txt.gz", "/home/sacha/refGene.db")

    app = QApplication(sys.argv)

    view = GeneView()
    view.show()

    app.exec_()
