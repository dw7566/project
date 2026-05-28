from . import MZMfitting

def runfile(path, save, show, csv_save):
    """
    This function is now a wrapper around the main analysis script.
    The original parameters are no longer used as the new script handles
    configuration internally.
    """
    # The path from run.py is a glob pattern, like 'data/**/*.xml'
    # The main function in ivfitting expects the root data directory.
    # The default is 'data', which is correct for this project structure.
    ivfitting.main()
