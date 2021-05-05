# Standard imports
import functools
import math
import csv
import io
import time
import itertools as it
from collections import defaultdict
import copy
import sys
import string
import urllib.request  # STRANGE: CANNOT IMPORT URLLIB ALONE
from logging import DEBUG

# dependency
import cachetools

# Qt imports
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *

# Custom imports
from cutevariant.core.querybuilder import build_sql_query
from cutevariant.core import sql

from cutevariant.core import command as cmd
from cutevariant.gui import plugin, FIcon, formatter, style
from cutevariant.gui.sql_thread import SqlThread
from cutevariant.gui.widgets import MarkdownEditor
import cutevariant.commons as cm

LOGGER = cm.logger()


class VariantModel(QAbstractTableModel):
    """VariantModel is a Qt model class which contains variant data from SQL DB.

    It loads paginated data and create an interface for Qt views and controllers.
    The model can group variants by (chr,pos,ref,alt).

    See Qt model/view programming for more information
    https://doc.qt.io/qt-5/model-view-programming.html

    Variants are stored internally as a list of variants.
    By default, there is only one variant per row until a user selects a field
    from annotations or from multiple samples.
    Duplicated variants will be displayed in this case. It is advised to use the
    button "group" in the GUI to easily split between specific fields and common
    fields.

    Signals:
        variant_loaded(bool): Emit when variant are loaded
        count_loaded(bool): Emit when total count are loaded
        error_raised(str): Emit message when threads or something else encounter errors
    """

    # emit when variant results is loaded
    variant_loaded = Signal()
    variant_is_loading = Signal(bool)

    # emit when toutal count is loaded
    count_loaded = Signal()
    count_is_loading = Signal(bool)

    # Emit when all load has started
    load_started = Signal()

    # Emit when all data ( count + variant) has finished
    load_finished = Signal()

    error_raised = Signal(str)
    interrupted = Signal()

    DEFAUT_CACHE_SIZE = 1_048_576 * 32  # Default cache size of 32 Mo

    def __init__(self, conn=None, parent=None):
        super().__init__()
        self.limit = 50
        self.page = 1  #
        self.total = 0
        self.variants = []
        self.headers = []

        # Cache all database fields and their descriptions for tooltips
        # Field names as keys, descriptions as values
        self.fields_descriptions = None

        self.fields = ["chr", "pos", "ref", "alt", "ann.gene"]

        self.filters = dict()
        self.source = "variants"
        self.group_by = []
        self.having = {}
        self.order_by = None
        self.order_desc = False
        self.formatter = None
        self.debug_sql = None
        # Keep after all initialization
        self.conn = conn

        # Thread (1 for getting variant, 1 for getting count variant )
        self._load_variant_thread = SqlThread(self.conn)
        self._load_count_thread = SqlThread(self.conn)

        self._load_variant_thread.started.connect(
            lambda: self.variant_is_loading.emit(True)
        )
        self._load_variant_thread.finished.connect(
            lambda: self.variant_is_loading.emit(False)
        )
        self._load_variant_thread.result_ready.connect(self._on_variant_loaded)
        self._load_variant_thread.error.connect(self.error_raised)

        self._load_count_thread.started.connect(
            lambda: self.count_is_loading.emit(True)
        )
        self._load_count_thread.finished.connect(
            lambda: self.count_is_loading.emit(False)
        )
        self._load_count_thread.result_ready.connect(self._on_count_loaded)

        self._finished_thread_count = 0
        self._user_has_interrupt = False

        # Create results cache because Thread doesn't use the memoization cache from command.py.
        # This is because Thread create a new connection and change the function signature used by the cache.

        self._load_variant_cache = cachetools.LFUCache(
            maxsize=self.DEFAUT_CACHE_SIZE, getsizeof=sys.getsizeof
        )
        self._load_count_cache = cachetools.LFUCache(maxsize=1000)

    @property
    def conn(self):
        """ Return sqlite connection """
        return self._conn

    @conn.setter
    def conn(self, conn):
        """ Set sqlite connection """
        self._conn = conn
        if conn:
            # Note: model is initialized with None connection during start
            # Cache DB fields descriptions
            self.fields_descriptions = {
                field["name"]: field["description"]
                for field in sql.get_fields(self.conn)
            }
            LOGGER.debug("Init async thread")

            # Clear cache with new connection
            self.clear_all_cache()

            # Init Runnables (1 for each query type)
            self._load_variant_thread.conn = conn
            self._load_count_thread.conn = conn

            # self._load_count_thread.error.connect(self._on_error)

    def rowCount(self, parent=QModelIndex()):
        """Overrided : Return children count of index"""
        # If parent is root

        return len(self.variants)

    def columnCount(self, parent=QModelIndex()):
        """Overrided: Return column count of parent .

        Parent is not used here
        """

        #  Check integrity for unit test
        if parent == QModelIndex():
            return len(self.headers)
        return 0

    def clear_all_cache(self):
        """ clear cache """
        self.clear_count_cache()
        self.clear_variant_cache()

    def clear_variant_cache(self):
        self._load_variant_cache.clear()

    def clear_count_cache(self):
        self._load_count_cache.clear()

    def cache_size(self):
        """ Return total cache size """
        return self._load_variant_cache.currsize

    def max_cache_size(self):
        return self._load_variant_cache.maxsize

    def clear(self):
        """Reset the current model

        - clear variants list
        - total of variants is set to 0
        - emit variant_loaded signal
        """
        self.beginResetModel()
        self.variants.clear()
        self.total = 0
        self.endResetModel()
        self.variant_loaded.emit()

    def data(self, index: QModelIndex(), role=Qt.DisplayRole):
        """Overrided: return index data according role.
        This method is called by the Qt view to get data to display according Qt role.

        Params:
            index (QModelIndex): index from where your want to get data
            role (Qt.ItemDataRole): https://doc.qt.io/qt-5/qt.html#ItemDataRole-enum

        Examples:
            index = model.index(row=10, column = 1)
            chromosome_value = model.data(index)
        """

        # Avoid error
        if not index.isValid():
            return

        if self.variants and self.headers:

            column_name = self.headers[index.column()]

            # ---- Display Role ----
            if role == Qt.DisplayRole:
                value = self.variant(index.row())[column_name]
                if value is None:
                    return "NULL"
                else:
                    return str(self.variant(index.row())[column_name])

            if role == Qt.ToolTipRole:
                value = str(self.variant(index.row())[column_name])
                return value

        return None

    def headerData(self, section, orientation=Qt.Horizontal, role=Qt.DisplayRole):
        """Overrided: Return column name and display tooltips on headers

        This method is called by the Qt view to display vertical or horizontal header data.

        Params:
            section (int): row or column number depending on orientation
            orientation (Qt.Orientation): Qt.Vertical or Qt.Horizontal
            role (Qt.ItemDataRole): https://doc.qt.io/qt-5/qt.html#ItemDataRole-enum

        Examples:
            # return 4th column name
            column_name = model.headerData(4, Qt.Horizontal)

        Warnings:
            Return nothing if orientation is != from Qt.Horizontal,
            and role != Qt.DisplayRole
        """
        # Display columns headers
        if orientation == Qt.Horizontal:
            if role == Qt.DisplayRole:
                return self.headers[section]

            if role == Qt.ToolTipRole:
                # Field descriptions on headers
                # Note: fields are set in load()
                if self.fields_descriptions and section != 0:
                    field_name = self.fields[section - 1]
                    return self.fields_descriptions.get(field_name)

    def update_variant(self, row: int, variant: dict):
        """Update a variant at the given row with given content

        Update the variant in the GUI AND in the DB.

        Args:
            row (int): Row id of the variant that will be modified
            variant (dict): Dict of fields to be updated
        """
        # Update in database
        variant_id = self.variants[row]["id"]

        print("variant", variant)

        # Update all variant with same variant_id
        # Use case : When several transcript are displayed
        for row in self.find_row_id_from_variant_id(variant_id):
            print(row)
            left = self.index(row, 0)
            right = self.index(row, self.columnCount() - 1)

            if left.isValid() and right.isValid():
                # Get database id of the variant to allow its update operation
                variant["id"] = self.variants[row]["id"]
                sql.update_variant(self.conn, variant)
                self.variants[row].update(variant)
                self.dataChanged.emit(left, right)

    def find_row_id_from_variant_id(self, variant_id: int) -> list:
        """Find the ids of all rows with the same given variant_id

        Args:
            variant_id (int): Variant sql record id
        Returns:
            (list[int]): ids of rows
        """
        return [
            row_id
            for row_id, variant in enumerate(self.variants)
            if variant["id"] == variant_id
        ]

    def interrupt(self):
        """Interrupt current query if active

        This is a blocking function...

        call interrupt and wait for the error_raised signals ...
        If nothing happen after 1000 ms, by pass and continue
        If I don't use the dead time, it is waiting for an infinite time
        at startup ... Because at startup, loading is called 2 times.
        One time by the register_plugin and a second time by the plugin.show_event
        """

        interrupted = False

        if self._load_count_thread:
            if self._load_count_thread.isRunning():
                self._user_has_interrupt = True
                self._load_count_thread.interrupt()
                self._load_count_thread.wait(1000)
                interrupted = True

        if self._load_variant_thread:
            if self._load_variant_thread.isRunning():
                self._user_has_interrupt = True
                self._load_variant_thread.interrupt()
                self._load_variant_thread.wait(1000)
                interrupted = True

        if interrupted:
            self.interrupted.emit()

        # # Wait for exception ...
        # if loop:
        #     self.error_raised.connect(loop.quit)
        #     loop.exec_()

    def is_running(self):
        if self._load_variant_thread and self._load_count_thread:
            return (
                self._load_variant_thread.isRunning()
                or self._load_count_thread.isRunning()
            )

        return False

    def load(self):
        """Start async queries to get variants and variant count

        Called by:
            - on_change_query() from the view.
            - sort() and setPage() by the model.

        See Also:
            :meth:`loaded`
        """
        if self.conn is None:
            return

        if self.is_running():
            LOGGER.debug(
                "Cannot load data. Thread is not finished. You can call interrupt() "
            )

        LOGGER.debug("Start loading")

        offset = (self.page - 1) * self.limit

        # Add fields from group by
        # self.clear()  # Assume variant = []
        self.total = 0
        self._finished_thread_count = 0
        # LOGGER.debug("Page queried: %s", self.page)

        # Store SQL query for debugging purpose
        self.debug_sql = build_sql_query(
            self.conn,
            fields=self.fields,
            source=self.source,
            filters=self.filters,
            limit=self.limit,
            offset=offset,
            order_desc=self.order_desc,
            order_by=self.order_by,
            group_by=self.group_by,
            having=self.having,
        )

        print(self.debug_sql)
        # Create load_func to run asynchronously: load variants
        load_func = functools.partial(
            cmd.select_cmd,
            fields=self.fields,
            source=self.source,
            filters=self.filters,
            limit=self.limit,
            offset=offset,
            order_desc=self.order_desc,
            order_by=self.order_by,
            group_by=self.group_by,
            having=self.having,
        )

        # Create count_func to run asynchronously: count variants
        count_function = functools.partial(
            cmd.count_cmd,
            fields=self.fields,
            source=self.source,
            filters=self.filters,
            group_by=self.group_by,
        )

        # Start the run
        self._start_timer = time.perf_counter()

        # Create function HASH for CACHE
        self._count_hash = hash(
            count_function.func.__name__ + str(count_function.keywords)
        )
        self._variant_hash = hash(load_func.func.__name__ + str(load_func.keywords))

        self.load_started.emit()

        # Launch the first thread "count" or by pass it using the cache
        if self._count_hash in self._load_count_cache:
            self._load_count_thread.results = self._load_count_cache[self._count_hash]
            self._on_count_loaded()
        else:
            self._load_count_thread.start_function(count_function)

        # Launch the second thread "count" or by pass it using the cache
        if self._variant_hash in self._load_variant_cache:
            self._load_variant_thread.results = self._load_variant_cache[
                self._variant_hash
            ]
            self._on_variant_loaded()

        else:
            self._load_variant_thread.start_function(lambda conn: list(load_func(conn)))

    def _on_variant_loaded(self):
        """
        Triggered when variant_thread is finished

        """

        #  Compute time elapsed since loading

        self._end_timer = time.perf_counter()
        self.elapsed_time = self._end_timer - self._start_timer

        self.beginResetModel()
        self.variants.clear()

        # Add fields from group by
        for g in self.group_by:
            if g not in self.fields:
                self.fields.append(g)

        # Save cache
        self._load_variant_cache[
            self._variant_hash
        ] = self._load_variant_thread.results.copy()

        # Load variants
        self.variants = self._load_variant_thread.results
        if self.variants:
            # Set headers of the view
            self.headers = list(self.variants[0].keys())

        # self.total = self._load_count_thread.results["count"]

        self.endResetModel()
        self.variant_loaded.emit()

        #  Test if both thread are finished
        self._finished_thread_count += 1
        if self._finished_thread_count == 2:
            self.load_finished.emit()

    def _on_count_loaded(self):
        """
        Triggered when count_threaed is finished
        """

        # Save cache
        self._load_count_cache[
            self._count_hash
        ] = self._load_count_thread.results.copy()

        self.total = self._load_count_thread.results["count"]
        self.count_loaded.emit()

        #  Test if both thread are finished
        self._finished_thread_count += 1
        if self._finished_thread_count == 2:
            self.load_finished.emit()

    def hasPage(self, page: int) -> bool:
        """ Return True if <page> exists otherwise return False """
        return (page - 1) >= 0 and (page - 1) * self.limit < self.total

    def setPage(self, page: int):
        """ set the page of the model """
        if self.hasPage(page):
            self.page = page

    def nextPage(self):
        """ Set model to the next page """
        if self.hasPage(self.page + 1):
            self.setPage(self.page + 1)

    def previousPage(self):
        """ Set model to the previous page """
        if self.hasPage(self.page - 1):
            self.setPage(self.page - 1)

    def firstPage(self):
        """ Set model to the first page """
        self.setPage(1)

    def lastPage(self):
        """ Set model to the last page """

        self.setPage(self.pageCount())

    def pageCount(self):
        """ Return total page count """
        return math.ceil(self.total / self.limit)

    def sort(self, column: int, order):
        """Overrided: Sort data by specified column

        column (int): column id
        order (Qt.SortOrder): Qt.AscendingOrder or Qt.DescendingOrder

        """
        if column < self.columnCount():
            field = self.fields[column - 1]

            self.order_by = [field]
            self.order_desc = order == Qt.DescendingOrder
            self.load()

    def variant(self, row: int) -> dict:
        """Return variant data according index"""
        return self.variants[row]

    def is_variant_loading(self):
        if self._load_variant_thread:
            return self._load_variant_thread.isRunning()
        else:
            return False

    def is_count_loading(self):
        if self._load_count_thread:
            return self._load_count_thread.isRunning()
        else:
            return False


class VariantDelegate(QStyledItemDelegate):
    """Specify the aesthetic (style and color) of variants displayed on a view"""

    def __init__(self, parent=None):
        super().__init__(parent)

        self.formatter = None

    def paint(self, painter, option, index):
        """Paint with formatter if defined"""

        if self.formatter is None:
            return super().paint(painter, option, index)

        # Draw selections
        if option.state & QStyle.State_Enabled:
            bg = (
                QPalette.Normal
                if option.state & QStyle.State_Active
                or option.state & QStyle.State_Selected
                else QPalette.Inactive
            )
        else:
            bg = QPalette.Disabled

        if option.state & QStyle.State_Selected:
            painter.fillRect(option.rect, option.palette.color(bg, QPalette.Highlight))

        # Draw formatters
        option.rect = option.rect.adjusted(
            3, 0, 0, 0
        )  # Don't know why I need to adjust the left margin .. .
        self.formatter.paint(painter, option, index)


class LoadingTableView(QTableView):
    """Movie animation displayed on VariantView for long SQL queries executed
    in background.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self._is_loading = False

    def paintEvent(self, event: QPainter):

        if self.is_loading():
            painter = QPainter(self.viewport())

            painter.drawText(
                self.viewport().rect(), Qt.AlignCenter, self.tr("Loading ...")
            )

        else:
            super().paintEvent(event)

    def start_loading(self):
        self._is_loading = True
        self.viewport().update()

    def stop_loading(self):
        self._is_loading = False
        self.viewport().update()

    def is_loading(self):
        return self._is_loading


class VariantView(QWidget):
    """A Variant view with controller like pagination

    Attributes:

        row_to_be_selected (int/None): (optional) Left groupby pane only:
            At the end of the :meth:`loaded` method, the first line is selected
            in order to refresh the current variant in the other pane.
            TL,DR: Select the first row if in grouped mode => refresh the right pane.

    Signals:
        error_raised (str): Emit message when async runnables encounter errors
        no_variant: Emitted when there is no variant to display.
            Used by left pane to clear right pane.
    """

    error_raised = Signal(str)
    no_variant = Signal()

    def __init__(self, parent=None, show_popup_menu=True):
        """
        Args:
            parent: parent widget
            show_popup_menu (boolean, optional: If False, disable the context menu
                over variants. For example the group pane should be disabled
                in order to avoid partial/false informations to be displayed
                in this menu.
                Hacky Note: also used to rename variants to groups in the page box...
        """
        super().__init__(parent)

        self.parent = parent
        self.show_popup_menu = show_popup_menu
        self.view = LoadingTableView()
        self.bottom_bar = QToolBar()

        # self.view.setFrameStyle(QFrame.NoFrame)
        self.view.setAlternatingRowColors(True)
        self.view.horizontalHeader().setStretchLastSection(True)
        self.view.setSelectionMode(QAbstractItemView.SingleSelection)
        self.view.verticalHeader().hide()

        self.view.setSortingEnabled(True)
        self.view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.view.setSelectionMode(QAbstractItemView.ExtendedSelection)
        ## self.view.setIndentation(0)
        self.view.setIconSize(QSize(22, 22))
        self.view.horizontalHeader().setSectionsMovable(True)

        # Setup model
        self.model = VariantModel()
        self.view.setModel(self.model)

        # Setup delegate
        self.delegate = VariantDelegate()
        self.view.setItemDelegate(self.delegate)

        # setup toolbar
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

        self.page_box = QLineEdit()
        self.page_box.setValidator(QIntValidator())
        self.page_box.setFixedWidth(50)
        self.page_box.setValidator(QIntValidator())
        self.page_box.setFocusPolicy(Qt.NoFocus)

        # Display nb of variants/groups and pages
        self.info_label = QLabel()
        self.time_label = QLabel()
        self.cache_label = QLabel()
        self.loading_label = QLabel()
        self.loading_label.setMovie(QMovie(cm.DIR_ICONS + "loading.gif"))

        self.bottom_bar.addAction(FIcon(0xF0A30), "sql", self.on_show_sql)
        self.bottom_bar.addWidget(self.time_label)
        self.bottom_bar.addWidget(self.cache_label)
        self.bottom_bar.addSeparator()
        self.bottom_bar.addWidget(spacer)

        # Add loading action and store action
        self.bottom_bar.addWidget(self.info_label)
        self.loading_action = self.bottom_bar.addWidget(self.loading_label)

        # self.bottom_bar.setIconSize(QSize(16, 16))
        # self.bottom_bar.setMaximumHeight(30)
        self.bottom_bar.setContentsMargins(0, 0, 0, 0)

        self.pagging_actions = []

        action = self.bottom_bar.addAction(FIcon(0xF0600), "<<", self.on_page_clicked)
        action.setShortcut(Qt.CTRL + Qt.Key_Left)
        action.setAutoRepeat(False)
        action.setToolTip(self.tr("First page (%s)" % action.shortcut().toString()))

        self.pagging_actions.append(action)

        action = self.bottom_bar.addAction(FIcon(0xF0141), "<", self.on_page_clicked)
        action.setShortcut(QKeySequence(Qt.Key_Left))
        action.setAutoRepeat(False)
        action.setToolTip(self.tr("Previous page (%s)" % action.shortcut().toString()))
        self.pagging_actions.append(action)

        self.bottom_bar.addWidget(self.page_box)
        self.page_box.returnPressed.connect(self.on_page_changed)

        action = self.bottom_bar.addAction(FIcon(0xF0142), ">", self.on_page_clicked)
        action.setShortcut(QKeySequence(Qt.Key_Right))
        action.setAutoRepeat(False)
        action.setToolTip(self.tr("Next page (%s)" % action.shortcut().toString()))
        self.pagging_actions.append(action)

        action = self.bottom_bar.addAction(FIcon(0xF0601), ">>", self.on_page_clicked)
        action.setShortcut(Qt.CTRL + Qt.Key_Right)
        action.setAutoRepeat(False)
        action.setToolTip(self.tr("Last page (%s)" % action.shortcut().toString()))
        self.pagging_actions.append(action)

        self.fav_action = QAction(FIcon(0xF00C0), self.tr("Toggle favorite "))
        self.fav_action.triggered.connect(lambda: self.update_favorites())
        self.fav_action.setShortcut(QKeySequence(Qt.Key_Space))
        self.fav_action.setShortcutContext(Qt.WidgetShortcut)
        self.fav_action.setAutoRepeat(False)
        self.fav_action.setToolTip(
            self.tr("Toggle favorite (%s)" % self.fav_action.shortcut().toString())
        )
        self.view.addAction(self.fav_action)

        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(self.view)
        main_layout.addWidget(self.bottom_bar)
        self.setLayout(main_layout)

        # Async stuff
        # broadcast focus signal
        self.row_to_be_selected = None
        # Get the status of async load: started/finished
        # self.model.loading.connect(self._set_loading)
        # Queries are finished (yes its redundant with loading signal...)
        self.model.variant_loaded.connect(self._on_variant_loaded)
        self.model.count_loaded.connect(self._on_count_loaded)

        self.model.count_is_loading.connect(self.set_tool_loading)
        self.model.variant_is_loading.connect(self.set_view_loading)

        # Connect errors from async runnables
        self.model.error_raised.connect(self.error_raised)
        #  connect double clicke
        self.view.doubleClicked.connect(self.on_double_clicked)

    # def _set_loading(self, active=True):
    #     """Slot to obtain the status of async load: started/finished

    #     Start/Stop the movie animation on the view.

    #     Keyword Args:
    #         active (bool): True if loaded, False if finished
    #     """
    #     self.setDisabled(active)
    #     if active:
    #         QTimer.singleShot(1000, self.start_loading_after_delay)
    #     else:
    #         self.view.stop_loading()

    # def start_loading_after_delay(self):
    #     self.setDisabled(self.model.is_loading)

    #     if self.model.is_loading:
    #         self.view.start_loading()
    #     else:
    #         self.view.stop_loading()

    def setModel(self, model: VariantModel):
        self.model = model
        self.view.setModel(model)

    def load(self, reset_page=False):

        if reset_page:
            self.model.page = 1
            self.model.order_by = None

        self.model.interrupt()

        self.model.load()

    def _on_variant_loaded(self):
        """Slot called when async queries from the model are finished
        (especially count of variants for page box).

        Signals:
            Emits no_variant signal.
        """
        # if self.row_to_be_selected is not None:
        #     # Left groupby pane only:
        #     # Select by default the first line in order to refresh the
        #     # current variant in the other pane
        #     if self.model.rowCount():
        #         self.select_row(0)
        #     else:
        #         self.no_variant.emit()

        self.time_label.setText(str(" Executed in %.2gs " % (self.model.elapsed_time)))
        cache = cm.bytes_to_readable(self.model.cache_size())
        max_cache = cm.bytes_to_readable(self.model.max_cache_size())
        self.cache_label.setText(str(" Cache {} of {}".format(cache, max_cache)))
        if LOGGER.getEffectiveLevel() != DEBUG:
            self.view.setColumnHidden(0, True)
        self.view.scrollToTop()

        #  Select first row
        if self.model.rowCount():
            self.select_row(0)

    def _on_count_loaded(self):

        self.page_box.clear()
        if self.model.pageCount() - 1 == 0:
            self.set_pagging_enabled(False)
        else:
            # self.page_box.addItems([str(i) for i in range(self.model.pageCount())])
            self.page_box.validator().setRange(1, self.model.pageCount())
            self.page_box.setText(str(self.model.page))
            self.set_pagging_enabled(True)

        text = self.tr("{} line(s) Page {} on {}")
        text = text.format(self.model.total, self.model.page, self.model.pageCount())
        self.info_label.setText(text)

        #  Set focus to view ! Otherwise it stay on page_box
        self.view.setFocus(Qt.ActiveWindowFocusReason)

    def set_formatter(self, formatter):
        self.delegate.formatter = formatter
        self.view.reset()

    @property
    def conn(self):
        return self.model.conn

    @conn.setter
    def conn(self, _conn):
        self.model.conn = _conn

    @property
    def fields(self):
        return self.model.fields

    @fields.setter
    def fields(self, _fields):
        self.model.fields = _fields

    @property
    def source(self):
        return self.model.source

    @source.setter
    def source(self, _source):
        self.model.source = _source

    @property
    def filters(self):
        return self.model.filters

    @filters.setter
    def filters(self, _filters):
        self.model.filters = _filters

    @property
    def group_by(self):
        return self.model.group_by

    @group_by.setter
    def group_by(self, _group_by):
        self.model.group_by = _group_by

    def on_page_clicked(self):

        action_text = self.sender().text()

        if action_text == "<<":
            fct = self.model.firstPage

        if action_text == ">>":
            fct = self.model.lastPage

        if action_text == "<":
            fct = self.model.previousPage

        if action_text == ">":
            fct = self.model.nextPage

        fct()
        self.load()

    def on_page_changed(self):
        """Slot called when page_box is modified and the user has pressed return key

        The validator is ok if this slot is called.
        """
        if self.page_box.text():
            page = int(self.page_box.text())
            self.model.setPage(page)
            self.load()

    def on_show_sql(self):
        """Display debug sql query"""
        msg_box = QMessageBox()
        msg_box.setWindowTitle("SQL debug")
        msg_box.setText(self.model.debug_sql)
        explain_query = []
        for rec in self.conn.execute("EXPLAIN QUERY PLAN " + self.model.debug_sql):
            explain_query.append(str(dict(rec)))

        msg_box.setDetailedText("\n".join(explain_query))

        msg_box.exec_()

    def set_pagging_enabled(self, active=True):
        self.page_box.setEnabled(active)
        for action in self.pagging_actions:
            action.setEnabled(active)

    def set_view_loading(self, active=True):
        self.view.setDisabled(active)

        def show_loading_if_loading():
            if self.model.is_variant_loading():
                self.view.start_loading()

        if active:
            QTimer.singleShot(2000, show_loading_if_loading)
        else:
            self.view.stop_loading()

    def set_tool_loading(self, active=True):

        if active:
            self.info_label.setText(
                self.tr("Counting all variants. This can take a while ... ")
            )
            self.loading_action.setVisible(True)
            self.loading_label.movie().start()
        else:
            self.loading_label.movie().stop()
            self.loading_action.setVisible(False)
            # self.info_label.setText("")

        self.bottom_bar.setDisabled(active)

    def set_loading(self, active=True):

        self.set_view_loading(active)
        self.set_tool_loading(active)

    def _get_links(self) -> list:
        """Get links from settings

        Return list of links from QSettings

        Exemples:
            {
            "name":"google",
            "url": "http://www.google.fr/q={}",
            "is_browser": True   # Open with browser
            "is_default": True   # is a default action
            }

        """
        links = []
        settings = QSettings()
        settings.beginGroup("variant_view")
        size = settings.beginReadArray("links")
        for index in range(size):
            settings.setArrayIndex(index)

            links.append(
                {
                    "name": settings.value("name"),
                    "url": settings.value("url"),
                    "is_default": settings.value("is_default", 0, type=bool),
                    "is_browser": settings.value("is_browser", 0, type=bool),
                }
            )

        settings.endArray()
        settings.endGroup()

        return links

    def contextMenuEvent(self, event: QContextMenuEvent):
        """Override: Show contextual menu over the current variant"""
        if not self.show_popup_menu:
            return

        menu = QMenu(self)
        pos = self.view.viewport().mapFromGlobal(event.globalPos())
        current_index = self.view.indexAt(pos)

        if not current_index.isValid():
            return

        current_variant = self.model.variant(current_index.row())
        full_variant = sql.get_one_variant(self.conn, current_variant["id"])
        # Update variant with currently displayed fields (not in variants table)
        full_variant.update(current_variant)

        # Copy action: Copy the variant reference ID in to the clipboard
        formatted_variant = "{chr}:{pos}-{ref}-{alt}".format(**full_variant)
        menu.addAction(
            FIcon(0xF014C),
            formatted_variant,
            functools.partial(
                QApplication.instance().clipboard().setText, formatted_variant
            ),
        )

        # Create favorite action
        fav_action = menu.addAction(
            self.tr("&Unmark favorite")
            if bool(full_variant["favorite"])
            else self.tr("&Mark as favorite")
        )
        fav_action.setCheckable(True)
        fav_action.setChecked(bool(full_variant["favorite"]))
        fav_action.toggled.connect(self.update_favorites)

        # Create classication action
        class_menu = menu.addMenu(self.tr("Classification"))
        for key, value in cm.CLASSIFICATION.items():

            action = class_menu.addAction(FIcon(cm.CLASSIFICATION_ICONS[key]), value)
            action.setData(key)
            on_click = functools.partial(self.update_classification, current_index, key)
            action.triggered.connect(on_click)

        # Create external links
        links_menu = menu.addMenu(self.tr("External links"))

        # Display only exec_ternal links with placeholders that can be mapped
        for link in self._get_links():
            url = self._create_url(link["url"], full_variant)
            if url:
                func_slot = functools.partial(self._open_url, url, link["is_browser"])
                action = links_menu.addAction(link["url"], func_slot)

        # Comment action
        on_edit = functools.partial(self.edit_comment, current_index)
        menu.addAction(self.tr("&Edit comment ..."), on_edit)

        # Edit menu
        menu.addSeparator()
        menu.addAction(
            FIcon(0xF018F), self.tr("&Copy"), self.copy_to_clipboard, QKeySequence.Copy
        )
        menu.addAction(
            FIcon(0xF0486),
            self.tr("&Select all"),
            self.select_all,
            QKeySequence.SelectAll,
        )

        # Display
        menu.exec_(event.globalPos())

    def _open_url(self, url: QUrl, in_browser=False):

        if in_browser:
            QDesktopServices.openUrl(url)

        else:
            urllib.request.urlopen(url.toString(), timeout=10)

    def _create_url(self, format_string: str, variant: dict) -> QUrl:
        """Create a link from a format string and a variant data

        Args:
            format_string (str): a string with format group like : http://www.google.fr?q={chr}
            variant (dict): a variant dict returned by sql.get_one_variant.

        Returns:
            QUrl: return url or return None

        """
        field_names = {
            name
            for _, name, _, _ in string.Formatter().parse(format_string)
            if name is not None
        }
        if field_names & variant.keys():
            # Full or partial mapping => accepted link
            return QUrl(format_string.format(**variant), QUrl.TolerantMode)
        else:
            return format_string

    def update_favorites(self, checked: bool = None):
        """Update favorite status of multiple selected variants

        if checked is None, toggle favorite

        Warnings:
            BE CAREFUL with this code, we try to limit useless SQL queries as
            much as possible.
        """

        # Do not update the same variant multiple times
        unique_ids = set()

        for index in self.view.selectionModel().selectedRows():
            if not index.isValid():
                continue

            # Get variant id
            variant = self.model.variants[index.row()]
            variant_id = variant["id"]

            if variant_id in unique_ids:
                continue
            unique_ids.add(variant_id)

            # Update GUI + DB
            if checked is None:
                #  Toggle checked box
                update_data = {"favorite": int(not bool(variant["favorite"]))}
            else:
                update_data = {"favorite": int(checked)}

            self.model.update_variant(index.row(), update_data)

    def update_classification(self, index: QModelIndex, value=3):
        """Update classification level of the variant at the given index"""
        if index.isValid():
            update_data = {"classification": int(value)}
            self.model.update_variant(index.row(), update_data)

    def edit_comment(self, index: QModelIndex):
        """Allow a user to add a comment for the selected variant"""
        if not index.isValid():
            return

        # Get comment from DB
        variant_data = sql.get_one_variant(
            self.model.conn, self.model.variant(index.row())["id"]
        )
        comment = variant_data["comment"] if variant_data["comment"] else ""

        editor = MarkdownEditor(default_text=comment)
        if editor.exec_() == QDialog.Accepted:
            # Save in DB
            self.model.update_variant(index.row(), {"comment": editor.toPlainText()})

            # Request a refresh of the variant_info plugin
            self.parent.mainwindow.refresh_plugin("variant_info")

    def select_all(self):
        """Select all variants in the view"""
        self.view.selectAll()

    def select_row(self, row):
        """Select the row with the given index

        Called for left pane by :meth:`loaded`.
        """
        index = self.view.model().index(row, 0)
        self.view.selectionModel().setCurrentIndex(
            index, QItemSelectionModel.SelectCurrent | QItemSelectionModel.Rows
        )

    def keyPressEvent(self, event: QKeyEvent):
        """
        Handles key press events on the VariantView
        Can be used to filter out unexpected behavior with KeySequence conflicts
        """
        if event.matches(
            QKeySequence.Copy
        ):  # Default behavior from QTableView only copies the index that the mouse hovers
            self.copy_to_clipboard()  # So copy to clipboard to get the expected behavior in contextMenuEvent
            event.accept()  # Accept the event before the QTableView handles it in a terrible way

    def copy_to_clipboard(self):
        """Copy the selected variant(s) into the clipboard

        The output data is formated in CSV (delimiter is `\t`)
        """
        # In memory file
        output = io.StringIO()
        # Use CSV to securely format the data
        writer = csv.DictWriter(output, delimiter="\t", fieldnames=self.model.headers)
        writer.writeheader()
        for index in self.view.selectionModel().selectedRows():
            # id col is not wanted
            variant = dict(self.model.variant(index.row()))
            if "id" in variant:
                del variant["id"]
            writer.writerow(variant)

        QApplication.instance().clipboard().setText(output.getvalue())
        output.close()

    def _open_default_link(self, index: QModelIndex):

        current_variant = self.model.variant(index.row())
        full_variant = sql.get_one_variant(self.conn, current_variant["id"])
        # Update variant with annotation data visible in the view ...
        full_variant.update(current_variant)
        #  get default link
        link = [i for i in self._get_links() if i["is_default"] is True]
        if not link:
            return

        link = link[0]

        url = self._create_url(link["url"], full_variant)
        if url:
            self._open_url(url, link["is_browser"])

    def on_double_clicked(self, index: QModelIndex):
        """
        React on double clicked
        TODO : duplicate code with ContextMenu Event ! Need to refactor a bit
        """

        self._open_default_link(index)


class VariantViewWidget(plugin.PluginWidget):
    """Contains the view of query with several controller"""

    variant_clicked = Signal(dict)
    LOCATION = plugin.CENTRAL_LOCATION

    # Emitted when the variants are both counted and loaded in the view. Provides the variants count and the elapsed time (in seconds) to execute the query
    variants_load_finished = Signal(int, float)

    ENABLE = True

    def __init__(self, parent=None):
        super().__init__(parent)

        # Create 2 Panes
        self.splitter = QSplitter(Qt.Horizontal)
        self.main_right_pane = VariantView(parent=self)

        self.splitter.addWidget(self.main_right_pane)

        # Make resizable TODO : ugly ... Make it nicer
        # def _resize1_section(l, o, n):
        #     self.groupby_left_pane.view.horizontalHeader().resizeSection(l, n)

        # def _resize2_section(l, o, n):
        #     self.main_right_pane.view.horizontalHeader().resizeSection(l, n)

        # self.main_right_pane.view.horizontalHeader().sectionResized.connect(_resize1_section)
        # self.groupby_left_pane.view.horizontalHeader().sectionResized.connect(
        #     _resize2_section
        # )

        # self.groupby_left_pane.view.setHorizontalHeader(self.main_right_pane.view.horizontalHeader())

        # Top toolbar
        self.top_bar = QToolBar()
        # PS: Actions with QAction::LowPriority will not show the text labels
        self.top_bar.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)

        # Save selection
        self.save_action = self.top_bar.addAction(
            FIcon(0xF0F87), self.tr("Save selection"), self.on_save_selection
        )
        self.save_action.setToolTip(
            self.tr("Save the current variants into a new selection")
        )
        # self.save_action.setPriority(QAction.LowPriority)

        # Refresh UI button
        action = self.top_bar.addAction(
            FIcon(0xF0450), self.tr("Refresh"), self.on_refresh
        )
        action.setToolTip(self.tr("Refresh the current list of variants"))
        # action.setPriority(QAction.LowPriority)

        # Interrupt current query
        action = self.top_bar.addAction(
            FIcon(0xF04DB), self.tr("Stop"), self.on_interrupt
        )
        action.setToolTip(self.tr("Stop current query"))

        #  edit action
        # TODO : move into view

        # Formatter tools
        self.top_bar.addSeparator()
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.top_bar.addWidget(spacer)

        # Add formatters to combobox
        self.formatter_combo = QComboBox()
        self.add_available_formatters()

        # Error handling
        self.log_edit = QLabel()
        self.log_edit.setMaximumHeight(30)
        self.log_edit.setStyleSheet(
            "QWidget{{background-color:'{}'; color:'{}'}}".format(
                style.WARNING_BACKGROUND_COLOR, style.WARNING_TEXT_COLOR
            )
        )
        self.log_edit.hide()
        self.log_edit.setFrameStyle(QFrame.StyledPanel | QFrame.Raised)

        # Setup layout
        main_layout = QVBoxLayout()
        main_layout.addWidget(self.top_bar)
        main_layout.addWidget(self.splitter)
        main_layout.addWidget(self.log_edit)

        self.setLayout(main_layout)

        # Make connection
        self.main_right_pane.view.selectionModel().currentRowChanged.connect(
            lambda x, _: self.on_variant_clicked(x)
        )

        self.main_right_pane.error_raised.connect(self.set_message)

        # Connect model's signal load_finished to this signal's load_finished
        self.main_right_pane.model.load_finished.connect(
            lambda: self.variants_load_finished.emit(
                self.main_right_pane.model.total,
                self.main_right_pane.model.elapsed_time,
            )
        )

        # Default group
        self.last_group = ["chr"]

    def add_available_formatters(self):
        """Populate the formatters

        Also recall previously selected formatters via config file.
        Default formatter is "SeqoneFormatter".
        """
        # Get previously selected formatter

        settings = QSettings()
        settings.beginGroup("variant_view")
        formatter_name = settings.value("formatter", "CutestyleFormatter")

        # Add formatters to combobox, a click on it will instantiate the class
        selected_formatter_index = 0
        for index, obj in enumerate(formatter.find_formatters()):
            self.formatter_combo.addItem(FIcon(0xF03D8), obj.DISPLAY_NAME, obj)
            if obj.__name__ == formatter_name:
                selected_formatter_index = index

        self.top_bar.addWidget(self.formatter_combo)
        self.formatter_combo.currentTextChanged.connect(self.on_formatter_changed)
        self.formatter_combo.setToolTip(self.tr("Change current style"))

        # Set the previously used/default formatter
        self.formatter_combo.setCurrentIndex(selected_formatter_index)

    def on_formatter_changed(self):
        """Activate the selected formatter and save it in settings

        Called when the current formatter is modified
        """

        settings = QSettings()
        settings.beginGroup("variant_view")

        formatter_class = self.formatter_combo.currentData()
        self.main_right_pane.set_formatter(formatter_class())
        # Save formatter setting
        settings.setValue("formatter", formatter_class.__name__)

    def on_open_project(self, conn):
        """Overrided from PluginWidget"""
        self.conn = conn
        # Set connections of models
        self.main_right_pane.conn = self.conn

        self.on_refresh()

    def on_refresh(self):
        """Overrided from PluginWidget"""
        # Save default data with current query attributes
        # See load(), we use this attr to restore fields after grouping

        if self.mainwindow:
            self.save_fields = self.mainwindow.state.fields
            self.save_filters = self.mainwindow.state.filters

            self.main_right_pane.fields = self.mainwindow.state.fields
            self.main_right_pane.source = self.mainwindow.state.source
            self.main_right_pane.filters = self.mainwindow.state.filters

        # Set formatter
        formatter_class = self.formatter_combo.currentData()

        #        formatter_class = next(formatter.find_formatters())
        self.main_right_pane.set_formatter(formatter_class())

        # Clear only the variant cache ! Because user can edit data
        self.main_right_pane.model.clear_variant_cache()

        # Load ui
        self.load(reset_page=True)

    def on_interrupt(self):
        self.main_right_pane.model.interrupt()

    def on_save_selection(self):
        """Triggered on 'save_selection' button

        - This slot creates a new selection (aka source) from the current state.
        - `source_editor` plugin will be refreshed.
        """

        selection_name, accept = QInputDialog.getText(
            self,
            self.tr("Create a new selection"),
            self.tr("Name of the new selection:"),
        )

        if not accept or not selection_name:
            return

        try:
            result = cmd.create_cmd(
                self.conn,
                selection_name,
                source=self.mainwindow.state.source,
                filters=self.mainwindow.state.filters,
                count=self.main_right_pane.model.total,
            )
            if result:
                self.mainwindow.refresh_plugin("source_editor")
        except Exception as e:
            LOGGER.exception(e)

            # Name already used
            QMessageBox.critical(
                self,
                self.tr("Error while creating selection"),
                self.tr("Error while creating selection '%s'; Name is already used")
                % selection_name,
            )
            # Retry
            self.on_save_selection()

    def on_no_variant(self):
        """Slot called when left pane hasn't any variant to display

        Clear right pane (nothing to select in it)
        """
        self.main_right_pane.model.clear()

    def load(self, reset_page=False):
        """Load all views

        Called by on_refresh, on_group_changed, and _show_group_dialog
        Display/hide groupby_left_pane on user demand.
        """

        # Hide potential errors
        self.log_edit.hide()

        # Grouped => ungrouped
        # Restore fields
        self.main_right_pane.fields = self.save_fields
        # Restore filters
        self.main_right_pane.filters = self.save_filters
        # Refresh model
        self.main_right_pane.load(reset_page)
        # print("saved right:", self.save_fields)

    def on_variant_clicked(self, index: QModelIndex):
        """React on variant clicked

        Args:
            index (QModelIndex): index into item models derived from
                QAbstractItemModel. U used by item views, delegates, and
                selection models to locate an item in the model.
        """

        if index.model() == self.main_right_pane.view.model():
            # Variant clicked on right pane
            variant = self.main_right_pane.model.variant(index.row())

        if self.mainwindow:
            self.mainwindow.state.current_variant = variant
            # Request a refresh of the variant_info plugin
            self.mainwindow.refresh_plugin("variant_info")
            self.mainwindow.refresh_plugin("gene_viewer")

    def _refresh_vql_editor(self):
        """Force refresh of vql_editor if loaded"""
        if "vql_editor" in self.mainwindow.plugins:
            # Request a refresh of the vql_editor plugin
            self.mainwindow.refresh_plugin("vql_editor")

    def copy(self):
        """Copy the selected variant(s) into the clipboard

        See Also: VariantView.copy_to_clipboard
        """
        self.main_right_pane.copy_to_clipboard()

    def select_all(self):
        """Select all variants in the view

        See Also: VariantView.select_all
        """
        self.main_right_pane.select_all()

    def set_message(self, message: str):
        """Show message error at the bottom of the view

        Args:
            message (str): Error message
        """

        if self.log_edit.isHidden():
            self.log_edit.show()

        icon_64 = FIcon(0xF0027, style.WARNING_TEXT_COLOR).to_base64(18, 18)

        self.log_edit.setText(
            """<div height=100%>
            <img src="data:image/png;base64,{}" align="left"/>
             <span> {} </span>
            </div>""".format(
                icon_64, message
            )
        )


if __name__ == "__main__":
    import sys
    from PySide2.QtWidgets import QApplication
    from cutevariant.core.importer import import_file, import_reader
    from cutevariant.core.reader import FakeReader, VcfReader
    from cutevariant.core import sql

    app = QApplication(sys.argv)

    conn = sql.get_sql_connection("/home/sacha/Dev/cutevariant/corpos3.db")

    w = VariantViewWidget()
    w.conn = conn
    w.main_right_pane.model.conn = conn
    w.main_right_pane.load()
    # w.main_view.model.group_by = ["chr","pos","ref","alt"]
    # w.on_refresh()

    w.show()

    app.exec_()
