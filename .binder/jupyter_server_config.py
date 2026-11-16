# Configuration for Jupyter Server in Binder
try:
    c = get_config()  # type: ignore[name-defined]  # noqa: F821

    # Increase IOPub message and data rate limits for 3D Trame / PyVista interactive visualization
    c.ServerApp.iopub_msg_rate_limit = 100000.0
    c.ServerApp.iopub_data_rate_limit = 100000000.0

    # Backward compatibility for legacy notebook server
    c.NotebookApp.iopub_msg_rate_limit = 100000.0
    c.NotebookApp.iopub_data_rate_limit = 100000000.0
except NameError:
    pass
