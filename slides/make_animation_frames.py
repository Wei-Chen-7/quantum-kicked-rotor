"""Split figures/04_momentum_distribution.gif into PNG frames for the Beamer deck.

The deck embeds the animation with the LaTeX animate package, which needs one
image file per frame. A GIF frame has at most 256 colors, so each frame is
saved as a palette PNG, which is several times smaller than full RGB; the
script checks that this loses nothing. Run from the repo root after
experiment 04:
    python slides/make_animation_frames.py
"""

import pathlib

import numpy as np
from PIL import Image

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
GIF = REPO_ROOT / "figures" / "04_momentum_distribution.gif"
OUT = REPO_ROOT / "slides" / "anim"


def main():
    OUT.mkdir(exist_ok=True)
    for old in OUT.glob("frame_*.png"):
        old.unlink()
    with Image.open(GIF) as gif:
        for i in range(gif.n_frames):
            gif.seek(i)
            rgb = gif.convert("RGB")
            frame = rgb.convert("P", palette=Image.Palette.ADAPTIVE, colors=256)
            if not np.array_equal(np.asarray(frame.convert("RGB")), np.asarray(rgb)):
                raise RuntimeError(f"frame {i} has more than 256 colors; a palette PNG would change it")
            frame.save(OUT / f"frame_{i:03d}.png", optimize=True)
        n_frames = gif.n_frames
    size = sum(path.stat().st_size for path in OUT.glob("frame_*.png"))
    print(f"Wrote {n_frames} frames to slides/anim/ (frame_000.png to frame_{n_frames - 1:03d}.png),"
          f" {size / 1e6:.1f} MB in total")


if __name__ == "__main__":
    main()
