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


def _contrast(a, b):
    def lum(c):
        ch = [v / 255 for v in (c.red(), c.green(), c.blue())]
        ch = [v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in ch]
        return 0.2126 * ch[0] + 0.7152 * ch[1] + 0.0722 * ch[2]
    hi, lo = sorted((lum(a), lum(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


@pytest.mark.parametrize("base, text", [("#ffffff", "#000000"), ("#1e1e1e", "#e6e6e6"), ("#000000", "#ffffff")])
def test_review_colours_readable_in_light_and_dark_themes(base, text):
    from PySide6.QtGui import QColor, QPalette
    from PySide6.QtWidgets import QApplication

    from tdtk_analyzer.gui_review import theme

    QApplication.instance() or QApplication([])
    pal = QPalette()
    pal.setColor(QPalette.Base, QColor(base))
    pal.setColor(QPalette.Text, QColor(text))
    th = theme(pal)
    assert _contrast(th["in_row"], th["text"]) >= 7          # table rows in the analysis
    assert _contrast(th["in_row"], th["decision"]) >= 4.5    # "added by you"
    assert _contrast(th["base"], th["decision"]) >= 4.5
    assert _contrast(th["base"], th["line"]) >= 3            # diameter plot
    assert _contrast(th["base"], th["axis"]) >= 4.5


def test_app_icon_has_every_size():
    from PySide6.QtWidgets import QApplication

    from tdtk_analyzer.gui import ICON_SIZES, app_icon

    QApplication.instance() or QApplication([])
    icon = app_icon()
    assert not icon.isNull()
    assert {s.width() for s in icon.availableSizes()} == set(ICON_SIZES)
