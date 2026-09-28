"""GUI smoke test (headless)."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")


def test_main_window_uses_system_fixed_font():
    from PySide6.QtGui import QFontInfo
    from PySide6.QtWidgets import QApplication

    from tdtk_analyzer.gui import MainWindow

    app = QApplication.instance() or QApplication([])
    w = MainWindow()
    font = w.log_view.font()
    # a hard-coded family such as "Monospace" does not exist on macOS/Windows and makes
    # Qt scan all fonts ("qt.qpa.fonts: Populating font family aliases ...")
    assert font.family() != "Monospace"
    assert QFontInfo(font).fixedPitch()
    assert w.tabs.count() == 4
    w.close()
    app.processEvents()


def test_quit_button_and_menu_close_the_window():
    from PySide6.QtGui import QKeySequence
    from PySide6.QtWidgets import QApplication

    from tdtk_analyzer.gui import MainWindow

    app = QApplication.instance() or QApplication([])
    w = MainWindow()
    w.show()
    app.processEvents()
    w.quit_btn.click()
    app.processEvents()
    assert not w.isVisible()

    w2 = MainWindow()
    w2.show()
    app.processEvents()
    quit_actions = [a for menu in w2.menuBar().actions() for a in menu.menu().actions() if "Quit" in a.text()]
    assert quit_actions and quit_actions[0].shortcut() == QKeySequence("Ctrl+Q")
    quit_actions[0].trigger()
    app.processEvents()
    assert not w2.isVisible()
