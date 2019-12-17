import os
from pathlib import Path


def remake_qt_ui():
    ui_file_dir = Path(
        os.path.join(os.path.dirname(__file__), "qt_ui_files/")
    ).resolve()
    ui_py_file_dir = Path(
        os.path.join(os.path.dirname(__file__), "qt_py_files/")
    ).resolve()
    for ui_file in ui_file_dir.iterdir():
        os.system(
            f"pyuic5 {ui_file.absolute()} -o {ui_py_file_dir.joinpath(ui_file.stem)}.py"
        )


if __name__ == "__main__":
    remake_qt_ui()
