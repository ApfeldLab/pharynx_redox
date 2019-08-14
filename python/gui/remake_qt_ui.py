import os
from pathlib import Path


def remake_qt_ui():
    ui_file_dir = Path("qt_ui_files/")
    ui_py_file_dir = Path("qt_py_files/")
    for ui_file in ui_file_dir.iterdir():
        os.system(
            f"pyuic5 {ui_file.absolute()} -o {ui_py_file_dir.joinpath(ui_file.stem)}.py"
        )


if __name__ == "__main__":
    remake_qt_ui()
