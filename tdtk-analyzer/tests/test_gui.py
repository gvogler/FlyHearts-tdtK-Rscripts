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
    assert w.tabs.count() == 5
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


def test_review_tab_adds_and_removes_traces(tmp_path):
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication

    from synthetic import heart_movie, write_cxd
    from tdtk_analyzer.gui import MainWindow
    from tdtk_analyzer.gui_review import C_IN
    from tdtk_analyzer.pipeline import Callbacks, Pipeline, Settings

    movies = tmp_path / "movies"
    movies.mkdir()
    for name, seed, per in [("MAYO0001_1_1wf_a.cxd", 0, 0.25), ("MAYO0002_2_1wf_a.cxd", 3, 0.22)]:
        write_cxd(str(movies / name), heart_movie(seed=seed, period=per))
    out = tmp_path / "out"
    Pipeline(Settings(movie_dir=str(movies), output_dir=str(out), run_analysis=False, workers=1,
                      min_file_size_mb=1), Callbacks(log=lambda m: None)).run()

    app = QApplication.instance() or QApplication([])
    w = MainWindow()
    w.output.setText(str(out))
    w.show()                                                 # shortcuts only work in a visible window
    rp = w.review_panel
    w.tabs.setCurrentWidget(rp)
    rp.refresh()
    rp.show_combo.setCurrentText("Not in the analysis")
    assert rp.table.rowCount() > 0
    rp.table.selectRow(0)
    rp.table.setFocus()
    app.processEvents()
    name = rp.table.item(0, C_IN).data(Qt.UserRole)
    QTest.keyClick(rp.table, Qt.Key_A)                       # keyboard shortcut: add
    app.processEvents()
    assert (out / "balled" / "excellent traces" / name).exists()
    assert rp.stale.isVisibleTo(rp)
    rp.show_combo.setCurrentText("Your decisions")
    assert rp.table.rowCount() == 1
    rp.table.item(0, C_IN).setCheckState(Qt.Unchecked)       # untick: remove again
    assert not (out / "balled" / "excellent traces" / name).exists()
    w.close()
    app.processEvents()
