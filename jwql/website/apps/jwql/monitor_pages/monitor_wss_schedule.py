"""This module contains code for the WSS Schedule Bokeh plots.

Authors
-------

    - Nicolas Flagey (based on Cosmic Rays script by Bryan Hilbert)

Use
---

    This module is intended to be imported and use as such:
    ::

        from jwql.website.apps.jwql import monitor_pages
        monitor_template = monitor_pages.WSSScheduleMonitor()

Bokeh figures will then be in:
        monitor_template.???
        monitor_template.???
"""

from datetime import datetime, timedelta
import os

from jwql.website.apps.jwql.monitor_models.schedule import WSSScheduleQueryHistory, WSSScheduleStats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


class WSSScheduleMonitor():
    def __init__(self, instrument, aperture):
        """Create instance

        Parameters
        ----------
        instrument : str
            Name of JWST instrument. e.g. 'nircam'

        aperture : str
            Name of aperture. e.g. 'NRCA1_FULL'
        """
        self._instrument = instrument
        self._aperture = aperture
        self.create_figures()

    def create_figures(self):
        """Wrapper function to create both the history and histogram plots
        for a given instrument/aperture.
        """
        # Get the data
        self.load_data()

        # Create the history plot
        self.history_figure = self.history_plot()

        # Create the histogram plot
        self.histogram_figure = self.histogram_plot()

    def get_histogram_data(self):
        """Get data required to create cosmic ray histogram from the
        database query.
        """



    def get_history_data(self):
        """Extract data on the history of cosmic ray numbers from the
        database query result
        """

    def histogram_plot(self):
        """Create the histogram figure of CR magnitudes.
        """

    def history_plot(self):
        """Create the plot of CR rates versus time
        """


    def identify_tables(self):
        """Determine which database tables as associated with
        a given instrument"""

    def load_data(self):
        """Query the database tables to get data"""

        # Determine which database tables are needed based on instrument
