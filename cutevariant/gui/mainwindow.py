"""Main window of Cutevariant"""
# Standard imports
import os
import sys
import sqlite3
import json
from pkg_resources import parse_version
from functools import partial
from logging import DEBUG

# Qt imports
from PySide2.QtCore import Qt, QSettings, QByteArray, QDir, QUrl, Signal
from PySide2.QtWidgets import *
from PySide2.QtGui import QIcon, QKeySequence, QDesktopServices


# Custom imports
from cutevariant.core import get_sql_connection, get_metadatas, command
from cutevariant.core.sql import get_database_file_name
from cutevariant.core.writer import CsvWriter, PedWriter
from cutevariant.gui import FIcon
from cutevariant.gui.state import State
from cutevariant.gui.wizards import ProjectWizard
from cutevariant.gui.settings import SettingsDialog
from cutevariant.gui.widgets.aboutcutevariant import AboutCutevariant
from cutevariant import commons as cm
from cutevariant.commons import (
    MAX_RECENT_PROJECTS,
    DIR_ICONS,
    MIN_AUTHORIZED_DB_VERSION,
)

from cutevariant.gui.export import ExportDialogFactory, ExportDialog

# Import plugins
from cutevariant.gui import plugin

LOGGER = cm.logger()


class MainWindow(QMainWindow):
    """Main window of Cutevariant"""

    variants_load_finished = Signal(int, float)
    selected_variant_changed = Signal()

    def __init__(self, parent=None):

        super().__init__(parent)

        self.setWindowTitle("Cutevariant")
        self.toolbar = self.addToolBar("maintoolbar")
        self.toolbar.setObjectName("maintoolbar")  # For window saveState
        self.setWindowIcon(QIcon(DIR_ICONS + "app.png"))
        self.setWindowFlags(Qt.WindowContextHelpButtonHint | self.windowFlags())

        # Keep sqlite connection
        self.conn = None

        # App settings
        self.app_settings = QSettings()

        # State variable of application
        # store fields, source, filters, group_by, having data
        # Often changed by plugins
        self.state = State()

        # Central workspace
        self.central_tab = QTabWidget()
        self.footer_tab = QTabWidget()
        self.vsplit = QSplitter(Qt.Vertical)
        self.vsplit.addWidget(self.central_tab)
        self.vsplit.addWidget(self.footer_tab)
        self.setCentralWidget(self.vsplit)

        # Status Bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)

        # Setup menubar
        self.setup_menubar()
        # Setup toobar under the menu bar
        self.setup_toolbar()

        # Register plugins
        self.plugins = {}  # dict of names (not titles) as keys and widgets as values
        self.dialog_plugins = {}  # dict of actions as keys and classes as values
        self.register_plugins()

        # Window geometry
        self.resize(600, 400)
        self.setGeometry(
            QApplication.instance().desktop().rect().adjusted(100, 100, -100, -100)
        )

        self.setTabPosition(Qt.AllDockWidgetAreas, QTabWidget.North)

        # If True, the GUI settings are deleted when the app is closed
        self.requested_reset_ui = False
        # Restores the state of this mainwindow's toolbars and dockwidgets
        self.read_settings()

        # Auto open recent projects
        recent = self.get_recent_projects()
        if recent and os.path.isfile(recent[0]):
            self.open(recent[0])

    def add_panel(self, widget, area=Qt.LeftDockWidgetArea):
        """Add given widget to a new QDockWidget and to view menu in menubar"""
        dock = QDockWidget()
        dock.setWindowTitle(widget.windowTitle())
        dock.setWidget(widget)
        dock.setStyleSheet("QDockWidget { font: bold }")
        # Keep the attached dock to allow further clean deletion
        widget.dock = dock
        # Set the objectName for a correct restoration after saveState
        dock.setObjectName(str(widget.__class__))

        self.addDockWidget(area, dock)
        action = dock.toggleViewAction()
        action.setIcon(widget.windowIcon())
        action.setText(widget.windowTitle())
        self.view_menu.addAction(action)

        self.toolbar.addAction(action)

    def register_plugins(self):
        """Dynamically load plugins to the window

        Wrapper of :meth:`register_plugin`

        Two types of plugins can be registered:
        - widget: Added to the GUI,
        - dialog: DialogBox accessible from the tool menu.

        `setting` type is handled separately in :meth:`cutevariant.gui.settings.load_plugins`
        """
        LOGGER.info("MainWindow:: Registering plugins...")

        # Get classes of plugins
        # Don't forget to skip disabled plugins
        for extension in plugin.find_plugins():
            self.register_plugin(extension)

    def register_plugin(self, extension):
        """Dynamically load plugins to the window

        Two types of plugins can be registered:
        - widget: Added to the GUI,
        - dialog: DialogBox accessible from the tool menu.

        `setting` type is handled separately in :meth:`cutevariant.gui.settings.load_plugins`

        Args:
            extension (dict): Extension dict. See :meth:`cutevariant.gui.plugin.find_plugins`
        """
        LOGGER.debug("Extension: %s", extension)

        name = extension["name"]
        title = extension["title"]
        displayed_title = title

        if "widget" in extension and extension["widget"].ENABLE:
            # New GUI widget
            plugin_widget_class = extension["widget"]

            # Setup new widget
            widget = plugin_widget_class(parent=self)
            widget.setDisabled(True)

            """Variant view is a special widget that loads variants and counts them.
            Every plugin should be warned whenever the variants get loaded, so we need to connect variant view's signals to mainwindow's signals
            Variant view should also tell whenever a variant has been selected in the view
            """
            if name == "variant_view":

                LOGGER.debug("Loading variant view")
                widget.variants_load_finished.connect(
                    lambda count, time: self.variants_load_finished.emit(count, time)
                )
                widget.variant_changed.connect(
                    lambda: self.selected_variant_changed.emit
                )
                LOGGER.debug("Connected variant view signals to mainwindow")

            if not widget.objectName():
                LOGGER.debug(
                    "widget '%s' has no objectName attribute; "
                    "=> fallback to extension name",
                    displayed_title,
                )
                widget.setObjectName(name)

            # Set title
            widget.setWindowTitle(displayed_title)

            # WhatsThis content
            long_description = extension.get("long_description")
            if not long_description:
                long_description = extension.get("description")
            widget.setWhatsThis(long_description)
            # Register (launch first init on some of them)
            widget.on_register(self)
            if self.conn:
                # If register occurs after a project is loaded we must set its
                # connection attribute
                widget.on_open_project(self.conn)

            # Init mainwindow via the constructor or on_register
            if widget.mainwindow != self:
                LOGGER.error(
                    "Bad plugin implementation, <mainwindow> plugin attribute is not set."
                )
                widget.on_close()
            else:
                # Add new plugin to plugins already registered
                self.plugins[name] = widget

                # Set position on the GUI
                if plugin_widget_class.LOCATION == plugin.DOCK_LOCATION:
                    self.add_panel(widget)

                if plugin_widget_class.LOCATION == plugin.CENTRAL_LOCATION:
                    self.central_tab.addTab(widget, widget.windowTitle())

                if plugin_widget_class.LOCATION == plugin.FOOTER_LOCATION:
                    self.footer_tab.addTab(widget, widget.windowTitle())

        if "dialog" in extension and extension["dialog"].ENABLE:
            # New tool menu entry
            plugin_dialog_class = extension["dialog"]

            # Add plugin to Tools menu
            dialog_action = self.tool_menu.addAction(displayed_title + "...")
            self.dialog_plugins[dialog_action] = plugin_dialog_class
            dialog_action.triggered.connect(self.show_dialog)

    def refresh_plugins(self, sender: plugin.PluginWidget = None):
        """Refresh all widget plugins

        It doesn't refresh the sender plugin, and not visible plugins.

        Args:
            sender (PluginWidget): from a plugin, you can pass "self" as argument
        """
        for plugin_obj in self.plugins.values():
            if plugin_obj is not sender and (
                plugin_obj.isVisible() or plugin_obj.REFRESH_ONLY_VISIBLE is False
            ):
                try:
                    plugin_obj.on_refresh()
                except Exception as e:
                    LOGGER.exception(e)

    def refresh_plugin(self, plugin_name: str):
        """Refresh a widget plugin identified by plugin_name

        It doesn't refresh the sender plugin

        Args:
            plugin_name (str): Name of the plugin to be refreshed
        """
        if plugin_name in self.plugins:
            plugin_obj = self.plugins[plugin_name]
            plugin_obj.on_refresh()

    def setup_menubar(self):
        """Menu bar setup: items and actions

        .. note:: Setup tools menu that could be dynamically augmented by plugins.
        """
        ## File Menu
        self.file_menu = self.menuBar().addMenu(self.tr("&File"))
        self.new_project_action = self.file_menu.addAction(
            FIcon(0xF01BA), self.tr("&New project"), self.new_project, QKeySequence.New
        )
        self.open_project_action = self.file_menu.addAction(
            FIcon(0xF095D),
            self.tr("&Open project..."),
            self.open_project,
            QKeySequence.Open,
        )
        ### Recent projects
        self.recent_files_menu = self.file_menu.addMenu(self.tr("Open recent"))
        self.setup_recent_menu()

        self.file_menu.addAction(QIcon(), self.tr("Export..."), self.on_export_pressed)

        self.export_menu = self.file_menu.addMenu(self.tr("Export as"))

        for export_format_name in ExportDialogFactory.get_supported_formats():

            action = self.export_menu.addAction(
                self.tr(f"Export as {export_format_name}..."), self.on_export_pressed
            )

            # Since there are several actions connected to the same slot, we need to pass the format to the receiver
            action.setData(export_format_name)

        # self.export_ped_action = self.file_menu.addAction(
        #     self.tr("Export pedigree PED/PLINK"), self.export_ped
        # )

        self.file_menu.addSeparator()

        self.file_menu.addAction(
            FIcon(0xF0193), self.tr("Save session ..."), self.save_session
        )
        self.file_menu.addAction(
            FIcon(0xF0770), self.tr("Restore session ..."), self.load_session
        )

        self.file_menu.addSeparator()
        ### Misc
        self.file_menu.addAction(
            FIcon(0xF0493), self.tr("Settings..."), self.show_settings
        )
        self.file_menu.addSeparator()
        self.file_menu.addAction(self.tr("&Quit"), self.close, QKeySequence.Quit)

        ## Edit
        # TODO: if variant_view plugin is not loaded, disable this menu entries...
        self.edit_menu = self.menuBar().addMenu(self.tr("&Edit"))
        self.edit_menu.addAction(
            FIcon(0xF018F),
            self.tr("&Copy selected variants"),
            self.copy_variants_to_clipboard,
            "ctrl+shift+c",
        )
        self.edit_menu.addSeparator()
        self.edit_menu.addAction(
            FIcon(0xF0486),
            self.tr("&Select all variants"),
            self.select_all_variants,
            QKeySequence.SelectAll,
        )

        ## View
        self.view_menu = self.menuBar().addMenu(self.tr("&View"))
        self.view_menu.addAction(self.tr("Reset widgets positions"), self.reset_ui)
        # Set toggle footer visibility action
        show_action = self.view_menu.addAction(
            FIcon(0xF018D), self.tr("Show bottom toolbar")
        )
        show_action.setCheckable(True)
        show_action.setChecked(True)
        show_action.toggled.connect(self.toggle_footer_visibility)
        self.view_menu.addSeparator()

        ## Tools
        self.tool_menu = self.menuBar().addMenu(self.tr("&Tools"))

        ## Help
        self.help_menu = self.menuBar().addMenu(self.tr("Help"))

        self.help_menu.addAction(
            FIcon(0xF02D6),
            self.tr("What's this"),
            QWhatsThis.enterWhatsThisMode,
            QKeySequence.WhatsThis,
        )
        self.help_menu.addAction(
            FIcon(0xF059F),
            self.tr("Documentation..."),
            partial(QDesktopServices.openUrl, QUrl(cm.WIKI_URL, QUrl.TolerantMode)),
        )
        self.help_menu.addAction(
            FIcon(0xF0A30),
            self.tr("Report a bug..."),
            partial(
                QDesktopServices.openUrl, QUrl(cm.REPORT_BUG_URL, QUrl.TolerantMode)
            ),
        )

        self.help_menu.addSeparator()
        self.help_menu.addAction(
            self.tr("About Qt..."), QApplication.instance().aboutQt
        )
        self.help_menu.addAction(
            QIcon(DIR_ICONS + "app.png"),
            self.tr("About Cutevariant..."),
            self.about_cutevariant,
        )

    def setup_toolbar(self):
        """Top tool bar setup: items and actions

        .. note:: Require selection_widget and some actions of Menubar
        """
        # Tool bar
        self.toolbar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        self.toolbar.addAction(self.new_project_action)
        self.toolbar.addAction(self.open_project_action)
        self.toolbar.addAction(
            FIcon(0xF02D7), self.tr("Help"), QWhatsThis.enterWhatsThisMode
        )
        self.toolbar.addSeparator()

    def open(self, filepath):
        """Open the given db/project file

        .. note:: Called at the end of a project creation by the Wizard,
            and by Open/Open recent projects slots.

        :param filepath: Path of project file.
        :type filepath: <str>
        """
        if not os.path.isfile(filepath):
            return

        # Save directory
        self.app_settings.setValue("last_directory", os.path.dirname(filepath))

        # Create connection
        self.conn = get_sql_connection(filepath)

        try:
            # DB version filter
            db_version = get_metadatas(self.conn).get("cutevariant_version")
            if db_version and parse_version(db_version) < parse_version(
                MIN_AUTHORIZED_DB_VERSION
            ):
                # Refuse to open blacklisted DB versions
                # Unversioned files are still accepted
                QMessageBox.critical(
                    self,
                    self.tr("Error while opening project"),
                    self.tr("File: {} is too old; please create a new project.").format(
                        filepath
                    ),
                )
                return

            self.open_database(self.conn)
            self.save_recent_project(filepath)

        except sqlite3.OperationalError as e:
            LOGGER.exception(e)
            QMessageBox.critical(
                self,
                self.tr("Error while opening project"),
                self.tr("File: {}\nThe following exception occurred:\n{}").format(
                    filepath, e
                ),
            )
            return

        # Show the project name in title and in status bar
        self.setWindowTitle("Cutevariant - %s" % os.path.basename(filepath))
        self.status_bar.showMessage(self.tr("{} opened").format(filepath))

    def open_database(self, conn):
        """Open the project file and populate widgets

        Args:
            conn (sqlite3.Connection): Sqlite3 Connection
        """
        self.conn = conn

        # Clear memoization cache for count_cmd
        # Clear State variable of application
        # store fields, source, filters, group_by, having data
        self.state = State()

        # Load previous window state for this project (file_path being the key for the settings)
        file_path = get_database_file_name(conn)

        # self.state = self.app_settings.value(f"{file_path}/last_state", State())

        for plugin_obj in self.plugins.values():
            plugin_obj.on_open_project(self.conn)
            plugin_obj.setEnabled(True)

    def save_recent_project(self, path):
        """Save current project into QSettings

        Args:
            path (str): path of project
        """
        paths = self.get_recent_projects()
        if path in paths:
            paths.remove(path)
        paths.insert(0, path)

        self.app_settings.setValue("recent_projects", paths[:MAX_RECENT_PROJECTS])

    def get_recent_projects(self):
        """Return the list of recent projects stored in settings

        Returns:
            list: Recent project paths
        """
        # Reload last projects opened
        recent_projects = self.app_settings.value("recent_projects", list())

        # Check if recent_projects is a list() (as expected)
        if isinstance(recent_projects, str):
            recent_projects = [recent_projects]

        # Return only existing project files
        return [p for p in recent_projects if os.path.exists(p)]

    def clear_recent_projects(self):
        """Slot to clear the list of recent projects"""
        self.app_settings.remove("recent_projects")
        self.setup_recent_menu()

    def setup_recent_menu(self):
        """Setup recent submenu with previously opened projects"""
        self.recent_files_menu.clear()
        for path in self.get_recent_projects():
            self.recent_files_menu.addAction(path, self.on_recent_project_clicked)

        self.recent_files_menu.addSeparator()
        self.recent_files_menu.addAction(self.tr("Clear"), self.clear_recent_projects)

    def on_recent_project_clicked(self):
        """Slot to load a recent project"""
        action = self.sender()
        LOGGER.debug(action.text())
        self.open(action.text())

    def new_project(self):
        """Slot to allow creation of a project with the Wizard"""
        wizard = ProjectWizard()
        if not wizard.exec_():
            return

        db_filename = (
            wizard.field("project_path")
            + QDir.separator()
            + wizard.field("project_name")
            + ".db"
        )
        try:
            self.open(db_filename)
        except Exception as e:
            self.status_bar.showMessage(e.__class__.__name__ + ": " + str(e))
            raise

    def open_project(self):
        """Slot to open an already existing project from a QFileDialog"""
        # Reload last directory used
        last_directory = self.app_settings.value("last_directory", QDir.homePath())

        filepath, _ = QFileDialog.getOpenFileName(
            self,
            self.tr("Open project"),
            last_directory,
            self.tr("Cutevariant project (*.db)"),
        )
        if filepath:
            try:
                self.open(filepath)
            except Exception as e:
                self.status_bar.showMessage(e.__class__.__name__ + ": " + str(e))
                raise

    def export_as_csv(self):
        """Export variants into CSV file"""
        # Reload last directory used
        last_directory = self.app_settings.value("last_directory", QDir.homePath())

        filepath, _ = QFileDialog.getSaveFileName(
            self, self.tr("Save project"), last_directory, self.tr("(*.csv)")
        )

        if filepath:
            with open(filepath, "w") as file:
                writer = CsvWriter(file)
                writer.save(self.conn)

    def export_ped(self):
        """Export samples into PED/PLINK file"""
        # Reload last directory used
        last_directory = self.app_settings.value("last_directory", QDir.homePath())

        # noinspection PyCallByClass
        filepath, _ = QFileDialog.getSaveFileName(
            self, self.tr("Save project"), last_directory, "(*.tfam)"
        )

        if filepath:
            filepath = filepath if filepath.endswith(".tfam") else filepath + ".tfam"

            with open(filepath, "w") as file:
                writer = PedWriter(file)
                writer.save(self.conn)

    def show_settings(self):
        """Slot to show settings window"""
        widget = SettingsDialog(self)
        widget.exec_()

    def show_dialog(self):
        """Show Plugin dialog box after a click on the tool menu"""
        action = self.sender()
        if action in self.dialog_plugins:
            # Get class object and instantiate it
            dialog_class = self.dialog_plugins[action]
            # Send SQL connection
            dialog = dialog_class(conn=self.conn)
            dialog.mainwindow = self
            dialog.exec_()

    def about_cutevariant(self):
        """Slot to show about window"""
        dialog_window = AboutCutevariant(self)
        dialog_window.exec_()

    def reload_ui(self):
        """Reload *without* reset the positions of the widgets"""
        geometry = self.saveGeometry()
        ui_state = self.saveState()
        self.reset_ui()
        # Reload positions
        self.restoreGeometry(geometry)
        self.restoreState(ui_state)

    def reset_ui(self):
        """Reset the positions of docks (and their widgets) to the default state

        All the widgets are deleted and reinstantiated on the GUI.
        GUI settings via QSettings are also reset.
        """
        # Set reset ui boolean (used by closeEvent)
        self.requested_reset_ui = True
        # Reset settings
        self.write_settings()

        # Remove widgets in QTabWidgets
        for index in range(self.central_tab.count()):
            self.central_tab.removeTab(index)

        for index in range(self.footer_tab.count()):
            self.footer_tab.removeTab(index)

        # Remove tool menu actions linked to the dialog plugins
        for action in self.dialog_plugins:
            # LOGGER.debug("Remove action <%s>", action.text())
            self.tool_menu.removeAction(action)
            action.deleteLater()

        # Purge widgets (central/footer and others) and related docks
        for plugin_obj in self.plugins.values():
            # LOGGER.debug("Remove plugin <%s>", plugin_obj)
            self.removeDockWidget(plugin_obj.dock)
            plugin_obj.on_close()

        # Clean registered plugins
        self.plugins = {}
        self.dialog_plugins = {}
        # Register new plugins
        self.register_plugins()
        # We reload everything, but do not forget the project's file name !
        self.open_database(self.conn)
        # Allow a user to save further modifications
        self.requested_reset_ui = False

    def deregister_plugin(self, extension):
        """Delete plugin from the UI; Called from app settings when a plugin is disabled.

        - dialogs plugins: Remove action from tool menu
        - widgets plugins: Delete widget and its dock if available via its
            on_close() method.

        Args:
            extension (dict): Extension dict. See :meth:`cutevariant.gui.plugin.find_plugins`
        """
        name = extension["name"]
        title = extension["title"]
        displayed_title = name if LOGGER.getEffectiveLevel() == DEBUG else title

        # Remove tool menu actions linked to the dialog plugins
        for action in self.dialog_plugins:
            if action.text() == displayed_title:
                # LOGGER.debug("Remove action <%s>", action.text())
                self.tool_menu.removeAction(action)
                action.deleteLater()
                del self.dialog_plugins[action]
                break

        # Purge widgets and related docks
        if name in self.plugins:
            plugin_obj = self.plugins[name]
            # LOGGER.debug("Remove plugin <%s>", plugin_obj)
            self.removeDockWidget(plugin_obj.dock)
            plugin_obj.on_close()
            del self.plugins[name]

    def copy_variants_to_clipboard(self):
        if "variant_view" in self.plugins:
            self.plugins["variant_view"].copy()

    def select_all_variants(self):
        """Select all elements in the current tab's view"""
        if "variant_view" in self.plugins:
            self.plugins["variant_view"].select_all()

    def closeEvent(self, event):
        """Save the current state of this mainwindow's toolbars and dockwidgets

        .. warning:: Make sure that the property objectName is unique for each
            QToolBar and QDockWidget added to the QMainWindow.

        .. note:: Reset windowState if asked.
        """
        self.write_settings()

        # Don't forget to tell all the plugins that the window is being closed
        for plugin_obj in self.plugins.values():
            plugin_obj.on_close()
        super().closeEvent(event)

    def save_session(self):
        """ save plugin state into a json file """

        filename, _ = QFileDialog.getSaveFileName(
            self,
            self.tr("Save the session"),
            QDir.home().path(),
            "Session file (*.session.json)",
        )

        if os.path.exists(filename):
            ret = QMessageBox.warning(
                self,
                "file already exists",
                "Overwrite?",
                QMessageBox.Yes | QMessageBox.No,
            )

            if ret == QMessageBox.No:
                return

        if not filename.endswith(".session.json"):
            filename = filename + ".session.json"

        #  write sessions
        session = {}
        for name, plugin in self.plugins.items():
            session[name] = plugin.to_json()

        #  write file
        with open(filename, "w") as file:
            json.dump(session, file)

    def load_session(self):
        filename, _ = QFileDialog.getOpenFileName(
            self,
            self.tr("Save the session"),
            QDir.home().path(),
            "Session file (*.session.json)",
        )

        if not os.path.exists(filename):
            return

        # read sessions
        with open(filename) as file:
            state = json.load(file)

        #  set plugins
        for name, plugin in self.plugins.items():
            if name in state:
                plugin.from_json(state[name])

    def write_settings(self):
        """Store the state of this mainwindow.

        .. note:: This methods is called by closeEvent
        """
        if self.requested_reset_ui:
            # Delete window state setting
            self.app_settings.remove("windowState")
        else:
            self.app_settings.setValue("geometry", self.saveGeometry())
            #  TODO: handle UI changes by passing UI_VERSION to saveState()
            self.app_settings.setValue("windowState", self.saveState())

    def read_settings(self):
        """Restore the state of this mainwindow's toolbars and dockwidgets

        .. note:: If windowState is not stored, current state is written.
        """
        # Init reset ui boolean
        self.requested_reset_ui = False

        self.restoreGeometry(QByteArray(self.app_settings.value("geometry")))
        #  TODO: handle UI changes by passing UI_VERSION to saveState()
        window_state = self.app_settings.value("windowState")
        if window_state:
            self.restoreState(QByteArray(window_state))
        else:
            # Setting has been deleted: set the current default state
            #  TODO: handle UI changes by passing UI_VERSION to saveState()
            self.app_settings.setValue("windowState", self.saveState())

    def toggle_footer_visibility(self, visibility):
        """Toggle visibility of the bottom toolbar and all its plugins"""
        self.footer_tab.setVisible(visibility)

    def on_export_pressed(self):
        """
        Slot called by any export action.
        Use QAction.setData() on the sender to set format name.
        Otherwise, this slot will guess the format based on the save file dialog
        """

        if not self.conn:
            QMessageBox.information(
                self,
                self.tr("Info"),
                self.tr("No project opened, nothing to export.\nAborting."),
            )
            return

        settings = QSettings()
        default_save_dir = settings.value("last_save_file_dir", QDir.homePath())

        # Supported export extensions to filter names in the save file dialog (all of them by default)

        factory = ExportDialogFactory()

        exts = factory.get_supported_formats()

        format_name = self.sender().data()

        # Narrow down available extensions based on the action that was triggered
        if format_name:
            exts = [format_name]

        filters_and_exts = {
            f"{ext.upper()} file (*.{ext})": ext for ext in exts
        }  # Hack to get only ext out of savefiledialog result's second element (associates the filter's message with the extension name)

        file_name, chosen_ext = QFileDialog.getSaveFileName(
            self,
            self.tr("Please chose a filename you'd like to save the database to"),
            default_save_dir,
            ";;".join(filters_and_exts.keys()),
        )

        if not file_name:
            QMessageBox.information(
                self,
                self.tr("Info"),
                self.tr("No file name specified, nothing will be written"),
            )
            return

        settings.setValue("last_save_file_dir", os.path.dirname(file_name))

        chosen_ext = filters_and_exts[
            chosen_ext
        ]  # Hacky, extracts extension from second element from getSaveFileName result

        # Automatic extension of file_name
        file_name = (
            file_name if file_name.endswith(chosen_ext) else f"{file_name}.{chosen_ext}"
        )

        export_dialog: ExportDialog = ExportDialogFactory.create_dialog(
            self.conn,
            chosen_ext,
            file_name,
            fields=self.state.fields,
            source=self.state.source,
            filters=self.state.filters,
        )

        # # TODO : refactor self.state
        # export_dialog.state = {
        # "fields" : self.state.fields,
        # "source": self.state.source,
        # "filters": self.state.filters,
        # }

        success = export_dialog.exec_()
        if success == QDialog.Accepted:
            QMessageBox.information(
                self,
                self.tr("Success!"),
                self.tr(f"Successfully saved {os.path.basename(file_name)}"),
            )
        else:
            QMessageBox.critical(
                self,
                self.tr("Error!"),
                self.tr(f"Cannot save file to {os.path.basename(file_name)}"),
            )

    # @Slot()
    # def on_query_model_changed(self):
    #     for name, _plugin in self.plugins.items():
    #         if _plugin.isVisible():
    #             _plugin.on_query_model_changed(self.query_model)

    # @Slot()
    # def on_variant_changed(self, variant):
    #     for name, _plugin in self.plugins.items():
    #         if _plugin.isVisible():
    #             _plugin.on_variant_changed(variant)


if __name__ == "__main__":
    app = QApplication(sys.argv)

    w = MainWindow()

    w.show()

    app.exec_()
