from subprocess import Popen

def load_jupyter_server_extension(nbapp):
    """serve the design.ipynb directory with bokeh server"""
    Popen(["panel", "serve", "design.ipynb", "--allow-websocket-origin=*"])