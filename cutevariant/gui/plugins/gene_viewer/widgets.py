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


def sql_to_vector_converter(data):
    return [int(i) for i in str(data, encoding="utf8").split(",") if i.isnumeric()]


def list_to_sql_adapter(data: list):
    return ",".join([d for d in data if d])


# This registration of types in sqlite3 is so important it must be set in this module before anything else
sqlite3.register_converter("VECTOR", sql_to_vector_converter)
sqlite3.register_adapter(list, list_to_sql_adapter)


def refGene_to_sqlite(ref_filename: str, db_filename: str):

    # Create databases
    conn = sqlite3.connect(db_filename)

    conn.execute(
        """
        CREATE TABLE refGene(  
        id INTEGER PRIMARY KEY,
        transcript TEXT, 
        tx_start INTEGER, 
        tx_end INTEGER,
        cds_start INTEGER,
        cds_end INTEGER,
        exon_starts VECTOR,
        exon_ends VECTOR,
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
                exonStarts = line[9].split(",")
                exonEnds = line[10].split(",")
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


class Gene:
    """Class to hold a representation of a gene, with structural data and variant annotations.
    Structural data include coding sequence (start, end), exon list (starts, ends), exon count and variants found on the gene.
    """

    def __init__(self):
        self.cds_start = None
        self.cds_end = None
        self.exon_starts = None
        self.exon_ends = None
        self.variants = []
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


# Defines available mouse modes:
MOUSE_SELECT_MODE = 0  # Mouse clicking and dragging causes rectangle selection
MOUSE_PAN_MODE = 1  # Mouse clicking and dragging causes view panning


class GeneView(QAbstractScrollArea):
    """A class to visualize variants on a gene diagram, showing introns, exons, and coding sequence."""

    def __init__(self, parent=None):
        super().__init__(parent)

        self.gene = Gene()

        # self.variants = [(117227892, "red"), (117243866, "red")]

        # style
        self.cds_height = 40
        self.exon_height = 30
        self.intron_height = 20

        # self.showMaximized()

        self.scale_factor = 1
        self.translation = 0

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.horizontalScrollBar().setRange(0, 0)

        self.horizontalScrollBar().valueChanged.connect(self.set_translation)

        self.resize(640, 200)
        # QScroller.grabGesture(self.viewport(), QScroller.LeftMouseButtonGesture)

        self.region = None

        self.region_brush = QBrush(QColor("#5094CB80"))
        # self.region_brush.setStyle(Qt.Dense4Pattern)
        self.region_pen = QPen(QColor("#1E3252"))

        self.viewport().setMouseTracking(True)

        self.mouse_mode = MOUSE_SELECT_MODE
        self.cursor_selects = False

    @property
    def mouse_mode(self) -> int:
        return self._mouse_mode

    @mouse_mode.setter
    def mouse_mode(self, value: int):

        if value == MOUSE_SELECT_MODE:
            self._mouse_mode = value
            QScroller.ungrabGesture(self.viewport())
            self.setCursor(Qt.SplitHCursor)
        elif value == MOUSE_PAN_MODE:
            self._mouse_mode = value
            QScroller.grabGesture(self.viewport(), QScroller.LeftMouseButtonGesture)
            self.setCursor(Qt.OpenHandCursor)
        else:
            raise ValueError(
                "Cannot set mouse mode to %s (accepted values are MOUSE_PAN_MODE,MOUSE_SELECT_MODE",
                str(value),
            )

    @property
    def cursor_selects(self) -> bool:
        return self._cursor_selects

    @cursor_selects.setter
    def cursor_selects(self, value: bool):
        self._cursor_selects = value

    def paintEvent(self, event: QPaintEvent):

        painter = QPainter()
        painter.begin(self.viewport())
        painter.setBrush(QBrush(QColor("white")))
        painter.drawRect(self.rect())

        # draw guide
        self.area = self.draw_area()
        painter.drawRect(self.area.adjusted(-2, -2, 2, 2))
        painter.setClipRect(self.area)

        if not (self.gene.tx_start and self.gene.tx_end):
            painter.drawText(
                self.area,
                Qt.AlignCenter,
                self.tr("No gene selected. Please chose one in the combobox"),
            )
            painter.end()
            return

        # self.marks = []

        # Draw variants
        self._draw_variants(painter)

        # Draw lollipops to highlight exon selection
        self._draw_exon_handles(painter)

        # Draw intron Background
        self._draw_introns(painter)

        # Draw exons
        self._draw_exons(painter)

        # Draw CDS
        self._draw_cds(painter)

        # Draw rect selection region
        self._draw_region(painter)

        # Draw cursor line
        # if self.mouse_mode == MOUSE_SELECT_MODE:
        #     line_x = self.mapFromGlobal(QCursor.pos()).x()
        #     painter.drawLine(line_x, 0, line_x, self.rect().height())

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
            painter.save()
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

                linearGrad = QLinearGradient(
                    QPoint(0, exon_rect.top()), QPoint(0, exon_rect.bottom())
                )
                linearGrad.setColorAt(0, QColor("#789FCC"))
                linearGrad.setColorAt(1, QColor("#789FCC").darker())
                brush = QBrush(linearGrad)
                painter.setBrush(brush)
                painter.drawRect(exon_rect)
            painter.restore()

    def _draw_exon_handles(self, painter: QPainter):
        if self.gene.exon_count:
            painter.save()
            painter.setClipRect(self.area)
            painter.setRenderHint(QPainter.Antialiasing)

            no_pen = QPen(Qt.NoPen)
            white_pen = QPen(QColor("white"))

            # Schedule rects for drawing, so we know wich one was selected before drawing it on top
            rects_to_draw = []
            mouse_hovered_handle = False
            self.selected_exon = None

            for i in range(self.gene.exon_count):

                start = self._dna_to_scroll(self.gene.exon_starts[i])
                end = self._dna_to_scroll(self.gene.exon_ends[i])

                base = (end + start) / 2

                head_width = 34
                head_height = self.exon_height

                head_rect = QRect(0, 0, head_width, head_height)

                head_rect.moveCenter(
                    QPoint(
                        base + self.area.left(),
                        self.area.center().y() + self.cds_height + 15,
                    )
                )

                # If the head rect is outside the boundaries, do nothing
                if self.area.intersected(head_rect) != head_rect:
                    continue

                # Now that the head rect is at its definitive place, let's check if the cursor hovers it (Thanks @dridk!)
                if head_rect.contains(self.mapFromGlobal(QCursor.pos())):
                    self.selected_exon = i
                    rects_to_draw.insert(0, (i, head_rect))

                else:
                    rects_to_draw.append((i, head_rect))

            for index, rect in rects_to_draw[::-1]:
                if index == self.selected_exon:
                    painter.setBrush(QBrush(self.palette().color(QPalette.Highlight)))

                else:
                    painter.setBrush(QBrush(QColor("white")))

                pen = QPen(QColor("darkgray"))
                pen.setWidth(2)
                painter.setPen(pen)
                painter.drawEllipse(rect)
                font = QFont()
                font.setBold(True)
                painter.setFont(font)
                painter.drawLine(
                    QPoint(rect.center().x(), rect.top()),
                    QPoint(rect.center().x(), self.draw_area().center().y()),
                )

                painter.drawText(rect, Qt.AlignCenter, str(index))

            if self.selected_exon != None:
                self.setCursor(Qt.PointingHandCursor)
            else:
                # Wired.. .TODO refactor
                self.mouse_mode = self.mouse_mode

            # selected_handle = None
            # for exon_index, handle in rects_to_draw:
            #     linearGrad = QLinearGradient(
            #         QPoint(0, handle.top()), QPoint(0, handle.bottom())
            #     )
            #     if exon_index == self.selected_exon:
            #         selected_handle = handle
            #         continue
            #     else:
            #         linearGrad.setColorAt(0, QColor("#789FCC"))
            #         linearGrad.setColorAt(1, QColor("#789FCC"))
            #     brush = QBrush(linearGrad)
            #     painter.setBrush(brush)
            #     painter.setPen(no_pen)
            #     painter.drawEllipse(handle)
            #     painter.setPen(white_pen)
            #     painter.drawText(handle, Qt.AlignCenter, str(exon_index))

            # # Always draw selected handle last, no matter what
            # if selected_handle:
            #     selected_handle.setWidth(selected_handle.width() * 1.2)
            #     selected_handle.setHeight(selected_handle.height() * 1.2)
            #     linearGrad = QLinearGradient(
            #         QPoint(0, selected_handle.top()),
            #         QPoint(0, selected_handle.bottom()),
            #     )
            #     linearGrad.setColorAt(0, QColor("#789FCC"))
            #     linearGrad.setColorAt(1, QColor("#789FCC").darker())
            #     brush = QBrush(linearGrad)
            #     painter.setBrush(brush)
            #     painter.setPen(no_pen)
            #     painter.drawEllipse(selected_handle)
            #     painter.drawPolygon(
            #         QPolygon(
            #             [
            #                 selected_handle.center() + QPoint(5, 0),
            #                 QPoint(
            #                     selected_handle.center().x(), self.area.height() / 2
            #                 ),
            #                 selected_handle.center() + QPoint(-5, 0),
            #             ]
            #         )
            #     )
            #     painter.setPen(white_pen)
            #     painter.drawText(
            #         selected_handle, Qt.AlignCenter, str(self.selected_exon)
            #     )

            painter.restore()

    def _draw_cds(self, painter: QPainter):
        painter.save()
        if self.gene.cds_start and self.gene.cds_end:
            painter.setClipRect(self.area)

            # We draw, on every exon, the CDS rectangle (if existing)
            for i in range(self.gene.exon_count):

                def overlap(interval1, interval2):
                    """
                    Given [0, 4] and [1, 10] returns (True,1, 4)
                    Given [0,10] and [15,20] return (False,0,0)
                    Thanks to https://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
                    for saving me time !
                    """
                    _overlaps = True
                    if interval2[0] <= interval1[0] and interval1[0] <= interval2[1]:
                        start = interval1[0]
                    elif interval1[0] <= interval2[0] and interval2[0] <= interval1[1]:
                        start = interval2[0]
                    else:
                        _overlaps = False
                        start, end = 0, 0

                    if interval2[0] <= interval1[1] <= interval2[1]:
                        end = interval1[1]
                    elif interval1[0] <= interval2[1] <= interval1[1]:
                        end = interval2[1]
                    else:
                        _overlaps = False
                        start, end = 0, 0

                    return (_overlaps, start, end)

                overlaps, start, end = overlap(
                    [self.gene.cds_start, self.gene.cds_end],
                    [self.gene.exon_starts[i], self.gene.exon_ends[i]],
                )

                if not overlaps:
                    continue

                start = self._dna_to_pixel(start)
                end = self._dna_to_pixel(end)

                start = self._pixel_to_scroll(start)
                end = self._pixel_to_scroll(end)

                cds_rect = QRect(0, 0, end - start, self.cds_height)
                cds_rect.moveTo(
                    start + self.area.left(),
                    self.area.center().y() - self.cds_height / 2,
                )

                linearGrad = QLinearGradient(
                    QPoint(0, cds_rect.top()), QPoint(0, cds_rect.bottom())
                )
                linearGrad.setColorAt(0, QColor("#194980"))
                linearGrad.setColorAt(1, QColor("#194980").darker(400))
                brush = QBrush(linearGrad)
                painter.setBrush(brush)
                painter.drawRect(cds_rect)
                painter.drawText(cds_rect, Qt.AlignCenter, str(i))

        painter.restore()

    def _draw_region(self, painter: QPainter):
        if self.region:
            painter.save()
            painter.setBrush(self.region_brush)
            painter.setPen(self.region_pen)
            painter.drawRect(self.region)
            painter.restore()

    def draw_area(self):
        return self.viewport().rect().adjusted(10, 10, -10, -10)

    def _pixel_to_dna(self, pixel: int) -> int:
        """Convert pixel coordinate to dna

        Args:
            pixel (int): coordinate in pixel

        Returns:
            int: coordinate in dna
        """
        tx_size = self.gene.tx_end - self.gene.tx_start
        scale = tx_size / self.area.width()
        return pixel * scale + self.gene.tx_start

    def _dna_to_pixel(self, dna: int) -> int:
        """Convert dna coordinate to pixel

        Args:
            dna (int): coordinate in dna

        Returns:
            int: coordinate in pixel
        """

        # normalize dna
        dna = dna - self.gene.tx_start
        tx_size = self.gene.tx_end - self.gene.tx_start
        scale = self.area.width() / tx_size
        return dna * scale

    def _dna_to_scroll(self, dna: int) -> int:
        return self._pixel_to_scroll(self._dna_to_pixel(dna))

    def _scroll_to_dna(self, pixel: int) -> int:
        return self._pixel_to_dna(self._scroll_to_pixel(pixel))

    def _pixel_to_scroll(self, pixel: int) -> int:
        """Convert pixel coordinate to scroll area

        Args:
            pixel (int): Coordinate in pixel

        Returns:
            int: Return coordinate in scroll area
        """
        return pixel * self.scale_factor - self.translation

    def _scroll_to_pixel(self, pos: int) -> int:
        """Convert scroll coordinate to pixel

        Args:
            pos (int): Coordinate from scrollarea

        Returns:
            int: Coordinate in pixel
        """
        return (pos + self.translation) / self.scale_factor

    @Slot(int)
    def set_scale(self, x: int):
        """Set view scale.

        This method rescale the view . It makes possible to zoom in or zoom out

        Args:
            x (int): scale factor.
        """
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

    def set_translation(self, x: int):
        """Set Translation

        This method translate the view.
        This method is called by scrollbar

        Args:
            x (int): translation factor between 0 and (transcript size in pixel * the scale factor)
        """
        self.translation = x
        self.viewport().update()

    def wheelEvent(self, event):
        """override

        Zoom in or zoom out with the mouse wheel

        """
        if event.delta() > 0:
            self.set_scale(self.scale_factor + 0.5)
        else:
            if self.scale_factor > 1:
                self.set_scale(self.scale_factor - 0.5)

    def zoom_to_dna_interval(self, start: int, end: int):
        """Sets the current view bounds to a rect around DNA sequence spanning from start to end

        Args:
            start (int): DNA start position in the current transcript
            end (int): DNA end position in the current transcript
        """
        dna_pixel_size = end - start

        tx_pixel_size = self.gene.tx_end - self.gene.tx_start
        scale = tx_pixel_size / dna_pixel_size
        self.set_scale(scale)
        self.horizontalScrollBar().setValue(self._dna_to_pixel(start) * scale)
        self.viewport().update()

    def keyPressEvent(self, event: QKeyEvent):
        # If any key is pressed, switch to mouse pan mode

        super().keyPressEvent(event)

        # if event.key() == Qt.Key_Control:
        #     self.setCursor(Qt.PointingHandCursor)
        #     self.cursor_selects = True
        #     self.mouse_mode = MOUSE_PAN_MODE

        # if event.key() == Qt.Key_Shift:
        #     self.setCursor(Qt.OpenHandCursor)
        #     self.cursor_selects = False
        #     self.mouse_mode = MOUSE_PAN_MODE

    def event(self, event):

        # Change cursor when gesture in panning ...
        if event.type() == QEvent.Gesture:
            for g in event.gestures():
                if g.state() == Qt.GestureUpdated:
                    self.setCursor(Qt.ClosedHandCursor)
                else:
                    self.setCursor(Qt.OpenHandCursor)

        super().event(event)

    def keyReleaseEvent(self, event: QKeyEvent):

        super().keyReleaseEvent(event)
        # self.setCursor(Qt.CrossCursor)
        # # Reset mouse mode to rect select (the default)
        # self.mouse_mode = MOUSE_SELECT_MODE
        # self.cursor_selects = False

    def mousePressEvent(self, event: QMouseEvent):

        if self.selected_exon != None:
            self.zoom_to_dna_interval(
                self.gene.exon_starts[self.selected_exon],
                self.gene.exon_ends[self.selected_exon],
            )

        if self.mouse_mode == MOUSE_SELECT_MODE:

            if event.button() == Qt.LeftButton:
                self.region = QRect(0, 0, 0, 0)
                self.region.setHeight(self.viewport().height())
                self.region.setLeft(event.pos().x())
                self.region.setRight(event.pos().x())

            if event.button() == Qt.RightButton:
                self.reset_zoom()

        super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent):

        if self.region:
            self.region.setRight(event.pos().x())

        # # Check if the mouse hovers an exon and, if so, selects it
        # if self.gene.tx_start and self.gene.tx_end:
        #     pass
        # else:
        #     self.selected_exon = -1

        self.viewport().update()

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent):

        if self.region:
            if self.region.normalized().width() < 5:
                self.region = None
                return

            self.region = self.region.normalized()

            if self.region.isEmpty():
                # There is no selection, dna_start - dna_end is too small, zooming will make no sense
                self.reset_zoom()
            else:
                dna_start = self._scroll_to_dna(
                    self.region.left() - self.draw_area().left()
                )

                dna_end = self._scroll_to_dna(
                    self.region.right() - self.draw_area().left()
                )

                print(self.region.width())
                self.zoom_to_dna_interval(dna_start, dna_end)
                self.region = None

        super().mouseReleaseEvent(event)

    def reset_zoom(self):
        """Resets viewer zoom."""
        if self.gene:
            self.zoom_to_dna_interval(self.gene.tx_start, self.gene.tx_end)
        else:
            self.set_scale(1)
            self.set_translation(0)


class GeneViewerWidget(plugin.PluginWidget):
    """Exposed class to manage VQL/SQL queries from the mainwindow"""

    # LOCATION = plugin.FOOTER_LOCATION
    ENABLE = True
    REFRESH_ONLY_VISIBLE = False

    def __init__(self, parent=None):
        super().__init__(parent)

        self.annotations_file_name = "/home/sacha/refGene.db"

        self.setWindowTitle(self.tr("Gene Viewer"))

        self.view = GeneView()
        self.combo = QComboBox()
        self.combo.setEditable(True)
        self.toolbar = QToolBar()
        # self.toolbar.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)

        # mouse_mode_selection_group = QActionGroup(self)

        # # ------------------------------------ Begin setting up set mouse grab action
        # self.mouse_grab_act = QAction(
        #     FIcon(0xF01BE),
        #     self.tr("Mouse grab"),
        # )
        # self.mouse_grab_act.setCheckable(True)
        # mouse_mode_selection_group.addAction(self.mouse_grab_act)
        # self.toolbar.addAction(self.mouse_grab_act)
        # # ------------------------------------- End setting up set mouse grab action

        # # ------------------------------------ Begin setting up set mouse exon select action
        # self.mouse_select_act = QAction(
        #     FIcon(0xF01BD),
        #     self.tr("Exon select"),
        # )
        # self.mouse_select_act.triggered.connect(
        #     lambda: self.set_view_mouse_mode(MOUSE_EXON_SELECTING)
        # )
        # self.mouse_select_act.setCheckable(True)
        # self.mouse_select_act.setChecked(True)
        # mouse_mode_selection_group.addAction(self.mouse_select_act)
        # self.toolbar.addAction(self.mouse_select_act)
        # # ------------------------------------- End setting up set mouse exon select action

        self.toolbar.addAction(
            FIcon(0xF1276), self.tr("Reset zoom"), lambda: self.view.reset_zoom()
        )
        self.toolbar.addWidget(self.combo)

        self.vlayout = QVBoxLayout(self)
        self.vlayout.addWidget(self.toolbar)
        self.vlayout.addWidget(self.view)

        self.annotations_conn = None

        self.model = QStringListModel()
        self.combo.setModel(self.model)
        self.load_combo()

        self.combo.currentIndexChanged.connect(self.on_combo_changed)

        self.current_variant = None

    def on_open_project(self, conn):
        self.conn = conn

    def on_refresh(self):

        # gene = self.combo.currentText()

        self.current_variant = sql.get_one_variant(
            self.conn,
            self.mainwindow.state.current_variant["id"],
            with_annotations=True,
        )

        gene = self.current_variant["annotations"][0]["gene"]

        self.combo.blockSignals(True)
        self.combo.setCurrentText(gene)
        self.combo.blockSignals(False)

        filters = copy.deepcopy(self.mainwindow.state.filters)

        fields = ["pos", "ann.gene"]

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
        # TODO: move this connection, shouldn't be here...
        self.annotations_conn = sqlite3.connect(
            self.annotations_file_name, detect_types=sqlite3.PARSE_DECLTYPES
        )
        print(self.annotations_conn)

        self.annotations_conn.row_factory = sqlite3.Row
        gene_names = [
            s["gene"] for s in self.annotations_conn.execute("SELECT gene FROM refGene")
        ]
        self.model.setStringList(gene_names)

    def on_combo_changed(self):

        if self.annotations_conn:

            gene = self.combo.currentText()

            gene_data = dict(
                self.annotations_conn.execute(
                    f"SELECT transcript,tx_start,tx_end,cds_start,cds_end,exon_starts,exon_ends,gene FROM refGene WHERE gene = '{gene}'"
                ).fetchone()
            )
            # self.view.gene.tx_start = gene_data["tx_start"]
            # self.view.gene.tx_end = gene_data["tx_end"]
            # self.view.gene.cds_start = gene_data["cds_start"]
            # self.view.gene.cds_end = gene_data["cds_end"]
            # self.view.gene.exon_starts = gene_data["exon_starts"]
            # self.view.gene.exon_ends = gene_data["exon_ends"]

            # Set all attributes of our gene from the query
            [
                setattr(self.view.gene, attr_name, gene_data[attr_name])
                for attr_name in gene_data
                if hasattr(self.view.gene, attr_name)
            ]

            self.view.viewport().update()


if __name__ == "__main__":

    pass

    import sys
    import sqlite3

    import os

    # try:
    #     os.remove("/home/sacha/refGene.db")
    # except:
    #     pass

    # refGene_to_sqlite("/home/sacha/refGene.txt.gz", "/home/sacha/refGene.db")

    app = QApplication(sys.argv)

    view = GeneView()
    view.show()
    view.resize(600, 500)

    annotations_conn = sqlite3.connect(
        "/home/sacha/refGene.db", detect_types=sqlite3.PARSE_DECLTYPES
    )
    annotations_conn.row_factory = sqlite3.Row
    gene_data = dict(
        annotations_conn.execute(
            "SELECT transcript,tx_start,tx_end,cds_start,cds_end,exon_starts,exon_ends,gene FROM refGene WHERE gene = 'GJB2'"
        ).fetchone()
    )

    # self.view.gene.tx_start = gene_data["tx_start"]
    # self.view.gene.tx_end = gene_data["tx_end"]
    # self.view.gene.cds_start = gene_data["cds_start"]
    # self.view.gene.cds_end = gene_data["cds_end"]
    # self.view.gene.exon_starts = gene_data["exon_starts"]
    # self.view.gene.exon_ends = gene_data["exon_ends"]

    # Set all attributes of our gene from the query
    [
        setattr(view.gene, attr_name, gene_data[attr_name])
        for attr_name in gene_data
        if hasattr(view.gene, attr_name)
    ]

    variants = [
        (view.gene.tx_start + 1000, 0),
        (view.gene.tx_start + 2000, 0),
        (view.gene.tx_start + 3000, 1),
    ]

    view.gene.variants = variants

    view.viewport().update()

    app.exec_()
