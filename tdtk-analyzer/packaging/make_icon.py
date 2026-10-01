"""Draws the app icon: a heart M-mode (kymograph) band with traced walls.

    python packaging/make_icon.py

Writes tdtk_analyzer/resources/icon_<size>.png (used by the running app) and
packaging/tdtk-analyzer.ico / .icns (used by PyInstaller for the app bundle).
"""

import math
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QPointF, QRectF, Qt  # noqa: E402
from PySide6.QtGui import (QColor, QGuiApplication, QImage, QLinearGradient, QPainter, QPainterPath,  # noqa: E402
                           QPen)

HERE = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(HERE, "..", "tdtk_analyzer", "resources")
SIZES = (16, 24, 32, 48, 64, 128, 256, 512, 1024)


def _half_width(u: float) -> float:
    """Heart half-width over one beat (u in [0, 1)): fast contraction, slower relaxation."""
    wide, narrow = 0.25, 0.085
    if u < 0.30:                                   # systole
        k = 0.5 - 0.5 * math.cos(math.pi * u / 0.30)
        return wide - (wide - narrow) * k
    if u < 0.75:                                   # relaxation
        k = 0.5 - 0.5 * math.cos(math.pi * (u - 0.30) / 0.45)
        return narrow + (wide - narrow) * k
    return wide                                    # diastole


def draw(size: int) -> QImage:
    img = QImage(size, size, QImage.Format_ARGB32)
    img.fill(Qt.transparent)
    p = QPainter(img)
    p.setRenderHint(QPainter.Antialiasing)
    p.scale(size / 100.0, size / 100.0)            # draw in a 100 x 100 box

    m = 2.0 if size >= 48 else 0.5
    tile = QPainterPath()
    tile.addRoundedRect(QRectF(m, m, 100 - 2 * m, 100 - 2 * m), 22, 22)
    bg = QLinearGradient(0, 0, 0, 100)
    bg.setColorAt(0, QColor("#232a3d"))
    bg.setColorAt(1, QColor("#0d1018"))
    p.fillPath(tile, bg)
    p.setClipPath(tile)

    beats, x0, x1, cy = 2.0, -2.0, 102.0, 50.0
    n = 240
    xs = [x0 + (x1 - x0) * i / n for i in range(n + 1)]
    hw = [100 * _half_width(((x - x0) / (x1 - x0) * beats + 0.62) % 1.0) for x in xs]
    top = [QPointF(x, cy - h) for x, h in zip(xs, hw)]
    bottom = [QPointF(x, cy + h) for x, h in zip(xs, hw)]

    band = QPainterPath(top[0])
    for q in top[1:]:
        band.lineTo(q)
    for q in reversed(bottom):
        band.lineTo(q)
    band.closeSubpath()
    fill = QLinearGradient(0, cy - 25, 0, cy + 25)
    fill.setColorAt(0.0, QColor("#a3135f"))
    fill.setColorAt(0.5, QColor("#ff4fc8"))
    fill.setColorAt(1.0, QColor("#a3135f"))
    p.fillPath(band, fill)

    wall = 4.2 if size >= 48 else 6.5               # thicker walls stay visible at 16-32 px
    pen = QPen(QColor("#4dff8a"), wall)
    pen.setCapStyle(Qt.RoundCap)
    pen.setJoinStyle(Qt.RoundJoin)
    p.setPen(pen)
    for edge in (top, bottom):
        path = QPainterPath(edge[0])
        for q in edge[1:]:
            path.lineTo(q)
        p.drawPath(path)
    p.end()
    return img


def main() -> int:
    QGuiApplication.instance() or QGuiApplication(sys.argv)
    os.makedirs(RESOURCES, exist_ok=True)
    files = {}
    for s in SIZES:
        path = os.path.join(RESOURCES, f"icon_{s}.png")
        draw(s).save(path)
        files[s] = path
    from PIL import Image

    imgs = {s: Image.open(f).convert("RGBA") for s, f in files.items()}
    # every size drawn on its own (small ones with thicker lines), not scaled down from the big one
    ico = [imgs[s] for s in (16, 24, 32, 48, 64, 128, 256)]
    ico[-1].save(os.path.join(HERE, "tdtk-analyzer.ico"), sizes=[im.size for im in ico], append_images=ico[:-1])
    imgs[1024].save(os.path.join(HERE, "tdtk-analyzer.icns"),
                    append_images=[imgs[s] for s in (16, 32, 64, 128, 256, 512)])
    os.remove(files[1024])                         # only needed for the .icns
    print("icons written to", os.path.normpath(RESOURCES), "and", HERE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
