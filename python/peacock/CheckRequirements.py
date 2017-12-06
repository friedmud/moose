def print_error(msg):
    print("\nError starting peacock: %s\n" % msg)
    print("You need to use the miniconda module. Run:")
    print("\tmodule load miniconda")

class ErrorObserver(object):
    """
    Just captures a VTK error event
    """
    def __init__(self):
        self._occurred = False
        self._message = None
        self.CallDataType = "string0"

    def __call__(self, obj, event, message):
        self._message = message
        self._occurred = True

    def errorOccurred(self):
        if self._occurred:
            self._occurred = False
            return True
        return False

    def message(self):
        return self._message

def check_vtk():
    """
    Start up a vtkRenderWindow and see if we get any errors.
    """
    try:
        import vtk
        w = vtk.vtkRenderWindow()
        w.SetSize(1, 1) # make it small so we don't notice it on startup
        e = ErrorObserver()
        w.AddObserver('ErrorEvent', e)
        w.Start()
        if e.errorOccurred():
            print(e.message())
            print("\nCouldn't show a VTK window.")
            print("This is most likely caused by outdated Mesa drivers.")
            print("You can try upgrading to Mesa version 10.6.5 or later or update your video card driver.")
            print("If you are running in a virtual machine, you can try turning OFF 3D acceleration.")
            print("\nYou can also try to use an alternative VTK package in a conda environment.")
            print("You can create a new conda environment by doing (with the miniconda module loaded):")
            print("\n\tconda create --clone root -n <env_name>")
            print("\tsource activate <env_name>")
            print("\tconda install -c idaholab vtk=7.1.1.opengl")
            print("\nThen you can try running peacock again.")
            print("Note that this VTK package may have various rendering problems and some peacock tests will fail.")
            return False
    except ImportError:
        print_error("Could not import vtk")
        return False
    return True

def check_matplotlib():
    """
    Make sure that matplotlib has a good Qt5 backend.
    """
    try:
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
        assert FigureCanvasQTAgg # silence pyflakes warning
    except ImportError:
        print_error("Could not import matplotlib")
        return False
    return True

def check_qt():
    """
    Make sure that we can import PyQt5
    """
    try:
        import PyQt5
        assert PyQt5 # silence pyflakes warning
    except ImportError:
        print_error("Could not import Qt5")
        return False
    return True

def check_pandas():
    """
    Make sure that we can import pandas. This is required for VectorPostprocessorReader
    """
    try:
        import pandas
        assert pandas # silence pyflakes warning
    except ImportError:
        print_error("Could not import pandas")
        return False
    return True

def has_requirements():
    return check_qt() and check_vtk() and check_matplotlib() and check_pandas()
