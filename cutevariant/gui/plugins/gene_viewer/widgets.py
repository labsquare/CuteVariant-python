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


# Each variant will have a color index to set a representation color
# Below is the definition of color mapping for these lollipops
VARIANT_LOLLIPOP_COLOR_MAP = {
    0: QColor("#FF0000"),
    1: QColor("#00FF00"),
    2: QColor("#0000FF"),
}


class Gene(QObject):
    """Class to hold a representation of a gene, with structural data and variant annotations.
    Structural data include coding sequence (start, end), exon list (starts, ends), exon count and variants found on the gene.
    """

    def __init__(self, parent: QObject = None):
        super().__init__(parent)
        self.cds_start = None
        self.cds_end = None
        self.exon_starts = None
        self.exon_ends = None
        self.variants = None
        self.tx_start = None
        self.tx_end = None

    @property
    def tx_start(self) -> int:
        return self._tx_start

    @tx_start.setter
    def tx_start(self, value: int):
        self._tx_start = value

    @property
    def tx_end(self) -> int:
        return self._tx_end

    @tx_end.setter
    def tx_end(self, value: int):
        self._tx_end = value

    @property
    def cds_start(self) -> int:
        return self._cds_start

    @cds_start.setter
    def cds_start(self, value: int):
        self._cds_start = value

    @property
    def cds_end(self) -> int:
        return self._cds_end

    @cds_end.setter
    def cds_end(self, value: int):
        self._cds_end = value

    @property
    def exon_starts(self) -> typing.List[int]:
        return self._exon_starts

    @exon_starts.setter
    def exon_starts(self, value: typing.List[int]):
        self._exon_count = len(value) if value else 0
        self._exon_starts = value

    @property
    def exon_ends(self) -> typing.List[int]:
        return self._exon_ends

    @exon_ends.setter
    def exon_ends(self, value: typing.List[int]):
        self._exon_count = len(value) if value else 0
        self._exon_ends = value

    @property
    def variants(self) -> typing.List[typing.Tuple[int, int]]:
        return self._variants

    @variants.setter
    def variants(self, value: typing.List[typing.Tuple[int, int]]):
        self._variants = value

    @property
    def exon_count(self) -> int:
        return self._exon_count


class GeneView(QAbstractScrollArea):
    """A class to visualize variants on a gene diagram, showing introns, exons, and coding sequence."""

    def __init__(self, parent=None):
        super().__init__(parent)

        self.gene = Gene(self)

        # self.variants = [(117227892, "red"), (117243866, "red")]

        # style
        self.cds_height = 40
        self.exon_height = 30
        self.intron_height = 20

        # self.showMaximized()

        self.scale_factor = 1
        self.translation = 0

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.horizontalScrollBar().setRange(0, 0)

        self.horizontalScrollBar().valueChanged.connect(self.set_translation)

        self.resize(640, 200)
        QScroller.grabGesture(self.viewport(), QScroller.LeftMouseButtonGesture)

    def paintEvent(self, event: QPaintEvent):

        painter = QPainter()
        painter.begin(self.viewport())
        painter.setBrush(QBrush(QColor("white")))
        painter.drawRect(self.rect())

        # draw guide
        self.area = self.draw_area()
        painter.drawRect(self.area.adjusted(-2, -2, 2, 2))

        # self.marks = []

        # draw rule
        painter.drawLine(
            self.area.left(),
            self.area.center().y(),
            self.area.right(),
            self.area.center().y(),
        )

        # Draw intron Background
        self._draw_introns(painter)

        # Draw exons
        self._draw_exons(painter)

        # Draw CDS

        # Draw variants

        self._draw_variants(painter)

        # charles_cds_start = (
        #     self._pixel_to_scroll(self._dna_to_pixel(self.gene.cds_start)) + self.area.left()
        # )
        # charles_cds_end = (
        #     self._pixel_to_scroll(self._dna_to_pixel(self.gene.cds_end)) + self.area.left()
        # )

        # rect = QRect()
        # rect.setLeft(charles_cds_start)
        # rect.setRight(charles_cds_end)
        # rect.setHeight(20)
        # painter.drawRect(rect)

        # Draw CDS

        painter.end()

    def _draw_introns(self, painter: QPainter):
        intron_rect = QRect(self.area)
        intron_rect.setHeight(self.intron_height)
        intron_rect.moveCenter(QPoint(self.area.center().x(), self.area.center().y()))
        linearGrad = QLinearGradient(
            QPoint(0, intron_rect.top()), QPoint(0, intron_rect.bottom())
        )
        linearGrad.setColorAt(0, QColor("#FDFD97"))
        linearGrad.setColorAt(1, QColor("#FDFD97").darker())
        brush = QBrush(linearGrad)
        painter.setBrush(brush)
        painter.drawRect(intron_rect)

    def _draw_variants(self, painter: QPainter):
        if self.gene.variants:
            painter.save()

            for variant in self.gene.variants:

                pos, color = variant
                color = VARIANT_LOLLIPOP_COLOR_MAP.get(color, QColor("#FF00FF"))

                x = self._pixel_to_scroll(self._dna_to_pixel(pos)) + self.area.left()
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

    def _draw_exons(self, painter: QPainter):
        if self.gene.exon_count:
            painter.setClipRect(self.area)
            for i in range(self.gene.exon_count):

                start = self._dna_to_pixel(self.gene.exon_starts[i])
                end = self._dna_to_pixel(self.gene.exon_ends[i])

                start = self._pixel_to_scroll(start)
                end = self._pixel_to_scroll(end)

                # for m in self.marks:
                #     print("mark", m)
                #     painter.drawLine(m, 0, m, self.area.bottom())

                # draw exons
                exon_rect = QRect(0, 0, end - start, self.exon_height)
                exon_rect.moveTo(
                    start + self.area.left(),
                    self.area.center().y() - self.exon_height / 2,
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

    def _pixel_to_dna(self, pixel: int):

        tx_size = self.gene.tx_end - self.gene.tx_start
        scale = tx_size / self.area.width()
        return pixel * scale + self.tx_start

    def _dna_to_pixel(self, dna: int):

        # normalize dna
        dna = dna - self.gene.tx_start
        tx_size = self.gene.tx_end - self.gene.tx_start
        scale = self.area.width() / tx_size
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

        self.filename = "/home/charles/refGene.db"

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

        # TODO: What if the filters have an '$or' operator as a root ?
        if "$and" not in filters:
            filters = {"$and": []}

        filters["$and"].append({"ann.gene": gene})

        variants = sql.get_variants(
            self.conn,
            fields,
            self.mainwindow.state.source,
            filters,
        )

        self.view.gene.variants = [(variant["pos"], "red") for variant in variants]

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

        self.view.gene.tx_start = txStart
        self.view.gene.tx_end = txEnd

        self.view.gene.cds_start = cdsStart
        self.view.gene.cds_end = cdsEnd

        self.view.gene.exon_starts = [int(i) for i in exonsStarts.split(",")[:-1]]
        self.view.gene.exon_ends = [int(i) for i in exonEnds.split(",")[:-1]]

        print("EXON STARTS :", self.view.gene.exon_starts)

        self.view.viewport().update()


if __name__ == "__main__":

    pass

    import sys
    import sqlite3

    import os

    # try:
    #     os.remove("/home/charles/refGene.db")
    # except:
    #     pass

    # refGene_to_sqlite("/home/charles/refGene.txt.gz", "/home/charles/refGene.db")

    app = QApplication(sys.argv)

    view = GeneView()
    view.show()

    app.exec_()
