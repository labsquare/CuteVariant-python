# Standard imports
import sqlite3
from typing import Callable

# Qt imports
from PySide2.QtCore import QThread, QObject, Signal

# Custom imports
from cutevariant.core.sql import get_sql_connection
from cutevariant.commons import logger

LOGGER = logger()


class SqlThread(QThread):
    """Used to execute SQL/VQL queries in a separated thread

    Attributes:
        - db_file (str): File path of the database.
        - async_conn (sqlite3.Connection): sqlite3 Connection
        - function (Callable): Function to be executed
        - caching_hash: A user defined hash of the request
        - results: Contain the result of the threaded function.
            `None`, as long as the function has not finished its execution done.



    Signals:
        - result_ready(): Emitted when results is available. If no results or error ,
            This signal is not emitted
        - error(str): Emitted when the function has encountered an error during
            its execution. The message is formatted with the type and the
            message of the exception.
        - finished() : inherits from QThread
        - starte() : inherits from QThread


    Notes:
        AutoDelete flag of such objects is set to False.
        This means that QThreadPool **will not** delete these objects after
        calling run() method.
    """

    error = Signal(str)
    result_ready = Signal()  # TODO rename to result_received

    def __init__(self, conn: sqlite3.Connection = None, function: Callable = None):
        """Init a Thread with sqlite connection and callable

        Notes:
            Since sqlite3 Connection objects are not thread safe, we create a new
            connection based on it. 



        Examples:

            >>> thread = SqlThread(conn)
            >>> thread.start_function(lambda conn: conn.execute("SELECT COUNT(*) FROM variants"))
            >>> thread.result_ready.connect(self.on_received)



        Args:
            conn (sqlite3.Connection): sqlite3 Connexion

        """

        super().__init__()

        self.conn = conn
        self._async_conn = None
        self.results = None
        self.function = function
        self.last_error = None

    @property
    def conn(self) -> sqlite3.Connection:
        """Return the Application thread connection
        
        Returns:
            sqlite3.Connection
        """
        return self._conn

    @property
    def async_conn(self) -> sqlite3.Connection:
        """Return the thread connection 
        It is a clone of application connection for thread safe
        
        Returns:
            TYPE: Description
        """
        return self._async_conn

    @conn.setter
    def conn(self, conn):
        if conn:
            self._conn = conn
            self._async_conn = None
            self.db_file = conn.execute("PRAGMA database_list").fetchone()["file"]

    def run(self):
        """Execute the function in a new thread

        Note:
            All the context of this function is executed into a separated thread.

        """
        self.last_error = None

        if self.function is None:
            LOGGER.exception("no function defined")
            return

        self._async_conn = get_sql_connection(self.db_file)
        assert self.async_conn

        try:
            LOGGER.debug("thread start ")
            self.results = self.function(self.async_conn)
            LOGGER.debug("Thread finished")
        except Exception as e:
            # LOGGER.exception(e)
            self.last_error = "%s: %s" % (e.__class__.__name__, str(e))
            self.async_conn.close()
            self.error.emit(self.last_error)
        else:
            self.async_conn.close()
            self.result_ready.emit()

        return

    def start_function(self, function: Callable, caching_hash=None):
        """Execute a function in the thread
        
        Args:
            function (Callable): Function must take con as arguments 

        Examples:
            >>> thread.exec_function(lambda conn: conn.execute("SELECT ..."))

        """
        assert isinstance(function, Callable)
        self.function = function
        self.start()

    def interrupt(self):
        """Interrupt the thread connection

        Use this method instead of terminate.
        """
        if self.async_conn:
            self.async_conn.interrupt()

    @property
    def function(self):
        return self._function

    @function.setter
    def function(self, value: Callable):
        # A valid function must be set
        self._function = value


if __name__ == "__main__":
    from PySide2.QtWidgets import QApplication
    import sys

    app = QApplication(sys.argv)

    app.exec_()
