"""Plugin to view/edit/remove/do set operations on selections in the database
from the GUI.
SourceEditorWidget class is seen by the user and uses selectionModel class
as a model that handles records from the database.
"""
# Qt imports
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *

# Custom imports
from cutevariant.core import sql, command
from cutevariant.core.reader import BedReader
from cutevariant.gui import plugin, FIcon
from cutevariant.commons import logger, DEFAULT_SELECTION_NAME


LOGGER = logger()

# =================== SELECTION_DIALOG ==========================


class SelectionDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.list_view = QListView(self)
        self.model = None

        self.vlayout = QVBoxLayout(self)

        self.ok_cancel = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self
        )
        self.ok_cancel.accepted.connect(self.accept)
        self.ok_cancel.rejected.connect(self.reject)

        self.vlayout.addWidget(self.list_view)
        self.vlayout.addWidget(self.ok_cancel)

    def set_model(self, model):
        self.model = model
        self.list_view.setModel(self.model)

    def get_chosen_selection(self):
        return self.model.data(self.list_view.currentIndex())


# =================== SELECTION MODEL ===========================
class SourceModel(QAbstractTableModel):
    """Model to store all selections from SQLite `selections` table.
    Usage:
        model = selectionModel()
        model.conn = conn
        model.load()
    """

    def __init__(self, conn=None):
        super().__init__()
        self.conn = conn
        self.records = []

    def rowCount(self, parent=QModelIndex()):
        """Overrided from QAbstractTableModel"""
        return len(self.records)

    def columnCount(self, parent=QModelIndex()):
        """Overrided from QAbstractTableModel"""
        return 2  # value and count

    def data(self, index: QModelIndex(), role=Qt.DisplayRole):
        """Return the data stored under the given role for the item referred
        to by the index (row, col)
        Overrided from QAbstractTableModel
        :param index:
        :param role:
        :type index:
        :type role:
        :return: None if no valid index
        :rtype:
        """

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            if index.column() == 0:
                return self.records[index.row()]["name"]

            if index.column() == 1:
                return self.records[index.row()]["count"]

        if role == Qt.DecorationRole:
            if index.column() == 0:
                return QIcon(FIcon(0xF04F1))

        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        """Return the data for the given role and section in the header with
        the specified orientation
        Overrided from QAbstractTableModel
        """
        if not self.records:
            return

        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return "selection"
            if section == 1:
                return "count"

        if orientation == Qt.Vertical and role == Qt.DisplayRole and section > 0:
            return self.records[section].get(
                "id", None
            )  # For debug purpose . displayed in vertical header

    def record(self, index: QModelIndex()):
        """Return Selection records by index"""
        if not index.isValid():
            return None
        return self.records[index.row()]

    def find_record(self, name: str):
        """Find a record by name
        @see: view.selectionModel()
        :return: QModelIndex
        """
        for idx, record in enumerate(self.records):
            if record["name"] == name:
                return self.index(idx, 0, QModelIndex())
        return QModelIndex()

    def remove_record(self, index: QModelIndex()):
        """Delete the selection with the given id in the database
        :return: True if the deletion has been made, False otherwise.
        :rtype: <boolean>
        """
        # Get selected record
        record = self.record(index)

        self.beginRemoveRows(QModelIndex(), index.row(), index.row())

        if sql.delete_selection(self.conn, record["id"]):
            # Delete in model; triggers currentRowChanged signal
            #  Magic... the record disapear ...  ??
            self.endRemoveRows()
            return True
        return False

    def edit_record(self, index, record: dict):
        """Edit the given selection in the database and emit `dataChanged` signal"""
        if sql.edit_selection(self.conn, record):
            self.dataChanged.emit(index, index)

    def load(self):
        """Load all selections into the model"""

        if self.conn is None:
            return

        self.beginResetModel()
        # Add all selections from the database
        # Dictionnary of all attributes of the table.
        #    :Example: {"name": ..., "count": ..., "query": ...}
        self.records = list(sql.get_selections(self.conn))
        self.endResetModel()


# =================== SELECTION VIEW ===========================


class SourceEditorWidget(plugin.PluginWidget):
    """Widget displaying the list of avaible selections.
    User can select one of them to update Query::selection
    """

    ENABLE = True

    def __init__(self, parent=None, conn=None):
        """
        Args:
            parent (QMainWindow): Mainwindow of Cutevariant, passed during
                plugin initialization.
            conn (sqlite3.connexion): Sqlite3 connexion
        """
        super().__init__(parent)

        self.setWindowIcon(FIcon(0xF10E4))
        # conn is always None here but initialized in on_open_project()
        self.model = SourceModel(conn)
        self.view = QTableView()
        self.view.setModel(self.model)
        self.view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.view.setSelectionMode(QAbstractItemView.SingleSelection)
        self.view.horizontalHeader().show()
        self.view.horizontalHeader().setStretchLastSection(False)
        self.view.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeToContents
        )
        self.view.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.view.horizontalHeader().hide()

        self.toolbar = QToolBar()
        # self.toolbar.setIconSize(QSize(16, 16))

        self.view.verticalHeader().hide()
        self.view.verticalHeader().setDefaultSectionSize(26)
        self.view.setShowGrid(False)
        self.view.setAlternatingRowColors(True)

        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.view)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.setLayout(layout)

        # Used to block signals during the insertions (if set to True)
        self.is_loading = False

        # Map the operations of context menu with an internal id not visible
        # from the user
        # This id is used by _create_set_operation_menu
        # Keys: user text; values: set operators
        # See menu_setup()
        self.set_operations_mapping = dict()

        # call on_current_row_changed when item selection changed
        self.view.selectionModel().currentRowChanged.connect(
            self.on_current_row_changed
        )
        # Keep it consistent with other plugins
        self.toolbar.setIconSize(QSize(16, 16))
        self.toolbar.addAction(FIcon(0xF0453), self.tr("Reload"), self.load)

        self.intersect_act = self.toolbar.addAction(
            FIcon(0xF0779), self.tr("Intersect"), self.on_set_operation_pressed
        )
        self.intersect_act.setData("&")
        self.difference_act = self.toolbar.addAction(
            FIcon(0xF077C), self.tr("Difference"), self.on_set_operation_pressed
        )
        self.difference_act.setData("-")
        self.union_act = self.toolbar.addAction(
            FIcon(0xF0778), self.tr("Union"), self.on_set_operation_pressed
        )
        self.union_act.setData("|")

        self.bed_intersect_button = QPushButton(
            FIcon(0xF0219), self.tr("Intersect with BED file ..."), self
        )
        self.bed_intersect_button.clicked.connect(self.create_selection_from_bed)
        self.bed_intersect_button.setFlat(True)
        self.toolbar.addWidget(self.bed_intersect_button)

        self.toolbar.addActions(
            [self.intersect_act, self.difference_act, self.union_act]
        )

    def on_open_project(self, conn):
        self.model.conn = conn
        self.conn = conn
        self.on_refresh()

    def on_refresh(self):
        self.view.selectionModel().blockSignals(True)
        self.model.load()
        self.source = self.mainwindow.state.source
        model_index = self.model.find_record(self.source)
        self.view.setCurrentIndex(model_index)
        self.view.selectionModel().blockSignals(False)

    @Slot()
    def on_current_row_changed(self):
        """This methods trigger the signal for the view
        Note:
            I don't broadcast the signal rowChanged to selectionChanged directly
            because I need to block signals for the view only
        """

        index = self.view.currentIndex()
        source = self.model.record(index)["name"]

        self.mainwindow.state.source = source
        self.mainwindow.refresh_plugins(sender=self)

    def menu_setup(self, locked_selection=False):
        """Setup popup menu. Used on special selections like the default one (named variants).

        Args:
            locked_selection (bool, optional): Allow to mask edit/remove actions. Defaults to False.

        Returns:
            [QMenu]: [A QMenu where all the actions that can be applied to the selected source]
        """
        menu = QMenu()

        if not locked_selection:
            menu.addAction(FIcon(0xF0900), self.tr("Edit"), self.edit_selection)

        #  Create action for bed
        menu.addAction(
            FIcon(0xF0219),
            self.tr("Intersect with BED file ..."),
            self.create_selection_from_bed,
        )

        # Set operations on selections: create mapping and actions
        set_icons_ids = (0xF0779, 0xF077C, 0xF0778)
        set_texts = (self.tr("Intersect"), self.tr("Difference"), self.tr("Union"))
        set_internal_ids = ("&", "-", "|")
        # Map the operations with an internal id not visible for the user
        # This id is used by _create_set_operation_menu
        # Keys: user text; values: internal ids
        self.set_operations_mapping = dict(zip(set_texts, set_internal_ids))

        # Create actions
        [
            menu.addMenu(self._create_set_operation_menu(FIcon(icon_id), text))
            for icon_id, text in zip(set_icons_ids, set_texts)
        ]

        if not locked_selection:
            menu.addSeparator()
            menu.addAction(FIcon(0xF0A7A), self.tr("Remove"), self.remove_selection)
        return menu

    def load(self):
        """Load selection model and update the view"""
        # Block signals during the insertions
        self.is_loading = True
        self.model.load()

        # Select record according to query.selection
        current_index = 0
        if self.source:
            current_index = self.model.find_record(self.source)
        self.view.setCurrentIndex(current_index)

        self.is_loading = False

    def ask_and_check_selection_name(self):
        """Get from the user a selection name validated or None

        TODO : create a sql.selection_exists(name) to check if selection already exists
        """
        name, success = QInputDialog.getText(
            self,
            self.tr("Type a name for selection"),
            self.tr("Selection name:"),
            QLineEdit.Normal,
        )
        if not success:
            return

        if name == DEFAULT_SELECTION_NAME:
            LOGGER.error(
                "SourceEditorWidget:save_current_query:: '%s' is a reserved name for a selection.",
                name,
            )
            self.mainwindow.status_bar.showMessage(
                self.tr("'%s' is a reserved name for a selection!") % name
            )
            QMessageBox.critical(
                self,
                self.tr("Error while creating the selection"),
                self.tr("'%s' is a reserved name for a selection!") % name,
            )
        elif name in {record["name"] for record in self.model.records}:
            LOGGER.error(
                "SourceEditorWidget:save_current_query:: '%s' is already used!", name
            )
            self.mainwindow.status_bar.showMessage(
                self.tr("'%s' is already used for a selection!") % name
            )
            QMessageBox.critical(
                self,
                self.tr("Error while creating the selection"),
                self.tr("'%s' is already used for a selection!") % name,
            )
        else:
            return name

    def on_set_operation_pressed(self):
        """Triggered when a 'Set operation' action is pressed **in the Toolbar**.
        Due to the availability of the data, this is not the exact same
        implementation as self._make_set_operation (though mostly copy-pasted)
        """
        action = self.sender()
        data = action.data()
        if not data:
            return

        # Following is a copy-paste from self._make_set_operation
        selection_name = self.ask_and_check_selection_name()
        if not selection_name:
            return

        second_operand = self.ask_for_existing_selection()
        if not second_operand:
            return

        # Get the action's internal data to know which set operation to do
        # See setData()
        set_operator = action.data()

        # Get the records and extract their database id to build 2 Selections objects
        # {'id': 2, 'name': 'dqzdezd', 'count': 3, 'query': 'SELECT variants.id ...}
        # {'id': 1, 'name': 'variants', 'count': 11, 'query': ''}
        record_1 = self.model.record(self.view.selectionModel().currentIndex())
        record_2 = self.model.record(self.model.find_record(second_operand))

        ret = command.set_cmd(
            self.model.conn,
            selection_name,
            record_1["name"],
            record_2["name"],
            set_operator,
        )
        if not ret:
            QMessageBox.critical(
                self,
                self.tr("Error while creating the selection"),
                self.tr("Error while creating the selection, please check the logs"),
            )
            self.mainwindow.status_bar.showMessage(
                self.tr("Fail to create the selection!")
            )
            return

        self.load()

    def ask_for_existing_selection(self):
        current_selection_name = self.model.record(self.view.currentIndex())["name"]
        records_names = set(record["name"] for record in self.model.records)
        records_names.remove(current_selection_name)

        dialog = SelectionDialog()
        dialog.set_model(QStringListModel(records_names, self))
        if dialog.exec_() == QDialog.Accepted:
            return dialog.get_chosen_selection()
        else:
            return None

    # def save_current_query(self):
    #     """Open a dialog box to save the current query into a selection
    #     TODO: désactivé?
    #     """
    #     selection_name = self.ask_and_check_selection_name()
    #     if selection_name:
    #         self.model.save_current_query(selection_name)

    def contextMenuEvent(self, event: QContextMenuEvent):
        """Overrided: Show popup menu at the cursor position"""
        if not self.model.records:
            return
        # Detect locked selection to mask edit/remove actions
        locked_selection = False
        current_index = self.view.selectionModel().currentIndex()
        record = self.model.record(current_index)
        if record["name"] == DEFAULT_SELECTION_NAME:
            locked_selection = True
        # Show the menu
        menu = self.menu_setup(locked_selection)
        menu.exec_(event.globalPos())

    def _create_set_operation_menu(self, icon, menu_name):
        """Dynamically add submenu with the given name to popup menu"""
        menu = QMenu(menu_name)
        menu.setIcon(icon)

        # Get all the names of selections except the current one
        current_selection_name = self.model.record(self.view.currentIndex())["name"]
        records_names = set(record["name"] for record in self.model.records)
        records_names.remove(current_selection_name)

        menu.addSeparator()
        # Create an action in the menu for all the remaining elements
        for record_name in records_names:
            action = menu.addAction(record_name)
            # This hidden data is used by _make_set_operation() to detect
            # which operation to do
            action.setData(self.set_operations_mapping[menu_name])
            action.triggered.connect(self._make_set_operation)

        return menu

    def _make_set_operation(self):
        """Do set operation
        .. note:: Called when a "set" submenu of the popup menu is triggered.
        """
        action = self.sender()

        selection_name = self.ask_and_check_selection_name()
        if not selection_name:
            return

        # Get the action's internal data to know which set operation to do
        # See setData()
        set_operator = action.data()

        # Get the records and extract their database id to build 2 Selections objects
        # {'id': 2, 'name': 'dqzdezd', 'count': 3, 'query': 'SELECT variants.id ...}
        # {'id': 1, 'name': 'variants', 'count': 11, 'query': ''}
        record_1 = self.model.record(self.view.selectionModel().currentIndex())
        record_2 = self.model.record(self.model.find_record(action.text()))

        ret = command.set_cmd(
            self.model.conn,
            selection_name,
            record_1["name"],
            record_2["name"],
            set_operator,
        )
        if not ret:
            QMessageBox.critical(
                self,
                self.tr("Error while creating the selection"),
                self.tr("Error while creating the selection, please check the logs"),
            )
            self.mainwindow.status_bar.showMessage(
                self.tr("Fail to create the selection!")
            )
            return

        self.load()

    def remove_selection(self):
        """Remove a selection from the database"""
        msg = QMessageBox(icon=QMessageBox.Warning)
        msg.setText(self.tr("Are you sure you want to remove this selection?"))
        msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)

        if msg.exec_() == QMessageBox.Yes:
            for index in self.view.selectionModel().selectedRows():
                self.model.remove_record(index)

            # Reload the UI, otherwise, the old selection is still in popup menu
            self.load()

    def edit_selection(self):
        """Update a selection in the database
        .. note:: We do not reload all the UI for this
        """
        current_index = self.view.selectionModel().currentIndex()

        selection_name = self.ask_and_check_selection_name()
        if current_index and selection_name:
            old_record = self.model.record(current_index)
            old_record["name"] = selection_name
            self.model.edit_record(current_index, old_record)

    def create_selection_from_bed(self):
        """Ask user for a bed file and create a new selection from it """
        # Reload last directory used
        app_settings = QSettings()
        last_directory = app_settings.value("last_directory", QDir.homePath())

        filepath, _ = QFileDialog.getOpenFileName(
            self,
            self.tr("Open BED file"),
            last_directory,
            self.tr("BED - Browser Extensible Data (*.bed)"),
        )
        if not filepath:
            return

        selection_name = self.ask_and_check_selection_name()
        if not selection_name:
            return

        current_index = self.view.selectionModel().currentIndex()
        current_selection = self.model.record(current_index)
        source = current_selection["name"]

        # Open bed intervals & create selection
        intervals = BedReader(filepath)
        sql.create_selection_from_bed(
            self.model.conn, source, selection_name, intervals
        )
        # Refresh UI
        self.model.load()


if __name__ == "__main__":

    import sqlite3
    import sys
    from PySide2.QtWidgets import QApplication

    app = QApplication(sys.argv)

    view = SourceEditorWidget(conn=sqlite3.connect("examples/test.snpeff.vcf.db"))
    view.show()

    app.exec_()
