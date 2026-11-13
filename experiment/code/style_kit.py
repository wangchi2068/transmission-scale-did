"""Shared publication-style helpers for experiment-suite figures."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Bbox
from matplotlib.text import Text
from mpl_toolkits.mplot3d.art3d import Text3D


PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "green_1": "#DDF3DE",
    "green_2": "#AADCA9",
    "green_3": "#8BCF8B",
    "red_1": "#F6CFCB",
    "red_2": "#E9A6A1",
    "red_strong": "#B64342",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
    "teal": "#42949E",
    "violet": "#9A4D8E",
}

DEFAULT_COLORS = [
    PALETTE["blue_main"],
    PALETTE["green_3"],
    PALETTE["red_strong"],
    PALETTE["teal"],
    PALETTE["violet"],
    PALETTE["neutral_mid"],
]

PALETTE_NMI_PASTEL = [
    "#484878",
    "#7884B4",
    "#B4C0E4",
    "#E4E4F0",
    "#E4CCD8",
    "#F0C0CC",
]

CJK_RE = re.compile(r"[\u3400-\u9fff]")
LATIN_OR_DIGIT_RE = re.compile(r"[A-Za-z0-9]")
TIMES_FAMILY = "Times New Roman"
SONGTI_FAMILIES = ("SimSun", "Songti SC", "STSong")
_ACTIVE_LANGUAGE: str | None = None
_ACTIVE_SONGTI_FAMILY: str | None = None


def _normalise_family(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _font_entry(path: str):
    resolved = Path(path).resolve()
    return next(
        (entry for entry in font_manager.fontManager.ttflist
         if Path(entry.fname).resolve() == resolved),
        None,
    )


def _font_weight(entry) -> int:
    if entry is None:
        return 400
    weight = entry.weight
    if isinstance(weight, str):
        return int(font_manager.weight_dict.get(weight.lower(), 400))
    return int(weight)


def _require_font_family(candidates: tuple[str, ...]) -> str:
    for family in candidates:
        properties = FontProperties(family=[family])
        try:
            path = font_manager.findfont(properties, fallback_to_default=False)
        except ValueError:
            continue
        resolved_name = FontProperties(fname=path).get_name()
        entry = _font_entry(path)
        if (
            _normalise_family(resolved_name) == _normalise_family(family)
            and _font_weight(entry) <= 500
        ):
            return family
    joined = ", ".join(candidates)
    raise RuntimeError(
        f"No regular/light publication face is resolvable for: {joined}. "
        "Install an individually addressable font file or use TikZ/PGF; do not "
        "silently substitute DejaVu, Noto, or a Black/Bold TTC face."
    )


def manuscript_font_properties(script: str) -> FontProperties:
    """Return an exact manuscript font for a separate Latin or CJK text artist."""
    key = script.lower()
    if key in {"latin", "english", "en"}:
        family = _require_font_family((TIMES_FAMILY,))
    elif key in {"cjk", "chinese", "zh", "cn"}:
        family = _require_font_family(SONGTI_FAMILIES)
    else:
        raise ValueError("script must be 'latin' or 'cjk'")
    return FontProperties(family=[family])


def apply_publication_style(
    font_size: int = 9,
    axes_linewidth: float = 0.6,
    language: str = "english",
) -> None:
    """Apply a shared journal-style matplotlib configuration."""
    global _ACTIVE_LANGUAGE, _ACTIVE_SONGTI_FAMILY
    if font_size < 7:
        raise ValueError("Publication figure font_size must be at least 7 pt")
    is_chinese = language.lower() in {"chinese", "zh", "zh-cn", "cn"}
    _require_font_family((TIMES_FAMILY,))
    _ACTIVE_LANGUAGE = "chinese" if is_chinese else "english"
    _ACTIVE_SONGTI_FAMILY = _require_font_family(SONGTI_FAMILIES) if is_chinese else None
    font_family = _ACTIVE_SONGTI_FAMILY if is_chinese else TIMES_FAMILY
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": [font_family],
        "font.size": font_size,
        "axes.labelsize": font_size,
        "axes.titlesize": font_size + 1,
        "xtick.labelsize": max(font_size - 1, 7),
        "ytick.labelsize": max(font_size - 1, 7),
        "legend.fontsize": max(font_size - 1, 7),
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": axes_linewidth,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.4,
        "lines.linewidth": 1.4,
        "lines.markersize": 4.5,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def add_panel_label(ax, label: str, xlabel: str | None = None) -> None:
    """Prefix a panel marker to the x-axis label without writing into the plot."""
    base_label = ax.get_xlabel() if xlabel is None else xlabel
    ax.set_xlabel(f"{label} {base_label}".strip())


def add_simulated_watermark(ax, x: float = 0.99, y: float = 0.02) -> None:
    ax.text(
        x, y, "simulated",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color="gray",
        alpha=0.7,
        style="italic",
    )


def _visible_text_artists(fig, include_text3d: bool = False):
    hidden_tick_labels = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            ticks = [*axis.get_major_ticks(), *axis.get_minor_ticks()]
            if not ax.axison or not axis.get_visible():
                for tick in ticks:
                    hidden_tick_labels.update((tick.label1, tick.label2))
                continue

            view_min, view_max = sorted(axis.get_view_interval())
            tolerance = max(abs(view_max - view_min) * 1e-10, 1e-12)
            for tick in ticks:
                if not view_min - tolerance <= tick.get_loc() <= view_max + tolerance:
                    hidden_tick_labels.update((tick.label1, tick.label2))

    for text in fig.findobj(match=Text):
        # Text3D reports unprojected display extents; validate it from the render.
        if isinstance(text, Text3D) and not include_text3d:
            continue
        if text in hidden_tick_labels:
            continue
        if not text.get_visible():
            continue
        if not text.get_text().strip():
            continue
        yield text


def assert_manuscript_fonts(fig) -> None:
    """Fail on font fallback or mixed-script text that cannot be exact-font typeset."""
    if _ACTIVE_LANGUAGE is None:
        raise RuntimeError("Call apply_publication_style() before saving a figure")

    fig.canvas.draw()
    approved_songti = {_normalise_family(name) for name in SONGTI_FAMILIES}
    required_times = _normalise_family(TIMES_FAMILY)
    for text in _visible_text_artists(fig, include_text3d=True):
        value = text.get_text()
        has_cjk = bool(CJK_RE.search(value))
        has_latin = bool(LATIN_OR_DIGIT_RE.search(value))
        if not has_cjk and not has_latin:
            continue
        if text.get_fontsize() < 7:
            raise RuntimeError(
                f"Figure text is below the 7 pt source minimum: {value[:80]!r}"
            )
        if has_cjk and has_latin:
            raise RuntimeError(
                f"Mixed CJK/Latin text must use separate text artists: {value[:80]!r}"
            )
        try:
            path = font_manager.findfont(
                text.get_fontproperties(), fallback_to_default=False
            )
        except ValueError as exc:
            raise RuntimeError(f"Unresolved figure font for {value[:80]!r}") from exc
        entry = _font_entry(path)
        actual = _normalise_family(FontProperties(fname=path).get_name())
        if _font_weight(entry) > 500:
            raise RuntimeError(
                f"Figure text must use a regular/light face, got weight "
                f"{_font_weight(entry)}: {value[:80]!r}"
            )
        if has_cjk and actual not in approved_songti:
            raise RuntimeError(
                f"Chinese figure text must use SimSun/Songti/STSong, got {actual!r}: {value[:80]!r}"
            )
        if has_latin and actual != required_times:
            raise RuntimeError(
                f"Latin figure text must use Times New Roman, got {actual!r}: {value[:80]!r}"
            )


def assert_no_text_collisions(fig, pad_px: float = 1.5) -> None:
    """Fail when visible text boxes overlap or leave the figure canvas."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    figure_box = fig.bbox
    boxes = []

    for text in _visible_text_artists(fig):
        bbox = text.get_window_extent(renderer=renderer)
        if bbox.width <= 0 or bbox.height <= 0:
            continue
        collision_bbox = bbox.expanded(
            (bbox.width + 2 * pad_px) / bbox.width,
            (bbox.height + 2 * pad_px) / bbox.height,
        )
        label = text.get_text().replace("\n", " ")[:80]
        if (
            bbox.x0 < figure_box.x0 - pad_px
            or bbox.y0 < figure_box.y0 - pad_px
            or bbox.x1 > figure_box.x1 + pad_px
            or bbox.y1 > figure_box.y1 + pad_px
        ):
            raise RuntimeError(f"Text clipped outside figure canvas: {label!r}")
        boxes.append((label, collision_bbox))

    for idx, (left_label, left_box) in enumerate(boxes):
        for right_label, right_box in boxes[idx + 1:]:
            if left_box.overlaps(right_box):
                raise RuntimeError(
                    f"Text collision: {left_label!r} overlaps {right_label!r}"
                )


def assert_no_legend_reference_line_collisions(fig, pad_px: float = 2.0) -> None:
    """Fail when an in-axes legend covers a dashed/dotted reference line."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        legend = ax.get_legend()
        if legend is None or not legend.get_visible():
            continue
        legend_box = legend.get_window_extent(renderer=renderer)
        legend_box = Bbox.from_extents(
            legend_box.x0 - pad_px,
            legend_box.y0 - pad_px,
            legend_box.x1 + pad_px,
            legend_box.y1 + pad_px,
        )
        for line in ax.lines:
            if not line.get_visible() or line.get_linestyle() not in {"--", ":", "-."}:
                continue
            display_path = line.get_path().transformed(line.get_transform())
            if display_path.intersects_bbox(legend_box, filled=False):
                label = line.get_label()
                raise RuntimeError(
                    f"Legend overlaps reference line: {label!r}"
                )


def add_centered_colorbar(
    fig,
    ax,
    mappable,
    label: str,
    orientation: str = "vertical",
    length: float = 0.62,
    thickness: float | None = None,
    gap: float | None = None,
    text_clearance_px: float = 14.0,
):
    """Place a colorbar outside and clear of its main axes and axis text."""
    fig.canvas.draw()
    box = ax.get_position()
    if orientation == "horizontal":
        thickness = 0.026 if thickness is None else thickness
        gap = 0.16 if gap is None else gap
        bounds = [
            box.x0 + box.width * (1 - length) / 2,
            box.y0 - gap,
            box.width * length,
            thickness,
        ]
    elif orientation == "vertical":
        thickness = 0.018 if thickness is None else thickness
        gap = 0.065 if gap is None else gap
        renderer = fig.canvas.get_renderer()
        axes = [ax.xaxis, ax.yaxis]
        if hasattr(ax, "zaxis"):
            axes.append(ax.zaxis)
        text_artists = [axis.label for axis in axes]
        for axis in axes:
            for tick in [*axis.get_major_ticks(), *axis.get_minor_ticks()]:
                text_artists.extend((tick.label1, tick.label2))
        right_edges = []
        for text in text_artists:
            if not text.get_visible() or not text.get_text().strip():
                continue
            text_box = text.get_window_extent(renderer=renderer)
            if 0 <= text_box.x1 <= fig.bbox.x1 * 1.5:
                right_edges.append(text_box.x1)
        text_safe_x = (
            (max(right_edges) + text_clearance_px) / fig.bbox.width
            if right_edges
            else 0.0
        )
        bounds = [
            max(box.x1 + gap, text_safe_x),
            box.y0 + box.height * (1 - length) / 2,
            thickness,
            box.height * length,
        ]
    else:
        raise ValueError("orientation must be 'horizontal' or 'vertical'")

    colorbar_ax = fig.add_axes(bounds)
    colorbar = fig.colorbar(mappable, cax=colorbar_ax, orientation=orientation)
    colorbar.set_label(label, labelpad=9 if orientation == "horizontal" else 11)
    return colorbar, colorbar_ax


def place_3d_zaxis_left(
    ax,
    labelpad: float = 10.0,
) -> None:
    """Place 3D z ticks and a closely spaced stable title on the left."""
    if not hasattr(ax, "zaxis") or not hasattr(ax.zaxis, "_axinfo"):
        raise TypeError("place_3d_zaxis_left requires a Matplotlib 3D axes")
    ax.zaxis._axinfo["juggled"] = (1, 2, 0)
    ax.zaxis.labelpad = labelpad
    label = ax.get_zlabel()
    previous = getattr(ax, "_left_zaxis_label_artist", None)
    if previous is not None:
        previous.remove()
    ax.zaxis.label.set_visible(False)
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tick_boxes = []
    for tick in [*ax.zaxis.get_major_ticks(), *ax.zaxis.get_minor_ticks()]:
        for tick_label in (tick.label1, tick.label2):
            if tick_label.get_visible() and tick_label.get_text().strip():
                tick_box = tick_label.get_window_extent(renderer=renderer)
                if tick_box.width > 0 and tick_box.height > 0:
                    tick_boxes.append(tick_box)
    if not tick_boxes:
        raise RuntimeError("Cannot position the 3D z-axis title without visible ticks")

    box = ax.get_position()
    label_artist = fig.text(
        0.0,
        (box.y0 + box.y1) / 2,
        label,
        rotation=90,
        ha="center",
        va="center",
        fontproperties=ax.zaxis.label.get_fontproperties(),
        color=ax.zaxis.label.get_color(),
    )
    fig.canvas.draw()
    label_box = label_artist.get_window_extent(renderer=renderer)
    tick_left = min(tick_box.x0 for tick_box in tick_boxes)
    label_center_x = tick_left - labelpad - label_box.width / 2
    label_center_x = max(label_box.width / 2 + 2.0, label_center_x)
    label_artist.set_x(label_center_x / fig.bbox.width)
    fig.canvas.draw()
    label_box = label_artist.get_window_extent(renderer=renderer)
    ax._left_zaxis_label_artist = label_artist
    ax._left_zaxis_label_gap_px = tick_left - label_box.x1


def assert_3d_zaxis_left(ax, colorbar_ax=None) -> None:
    """Verify the z-axis side contract and separation from a right colorbar."""
    if not hasattr(ax, "zaxis") or not hasattr(ax.zaxis, "_axinfo"):
        raise RuntimeError("3D z-axis QA requires a Matplotlib 3D axes")
    if tuple(ax.zaxis._axinfo.get("juggled", ())) != (1, 2, 0):
        raise RuntimeError("3D z-axis is not configured on the left")
    label_artist = getattr(ax, "_left_zaxis_label_artist", None)
    if label_artist is None or not label_artist.get_visible() or not label_artist.get_text().strip():
        raise RuntimeError("3D z-axis title is missing from the left side")
    if label_artist.get_position()[0] >= ax.get_position().x0:
        raise RuntimeError("3D z-axis title is not left of the 3D axes")
    label_gap = getattr(ax, "_left_zaxis_label_gap_px", None)
    if label_gap is None or label_gap < 4.0:
        raise RuntimeError("3D z-axis title overlaps its tick labels")
    if label_gap > 24.0:
        raise RuntimeError(f"3D z-axis title is too far from its ticks: {label_gap:.1f}px")
    if colorbar_ax is not None:
        main_box = ax.get_position()
        colorbar_box = colorbar_ax.get_position()
        if colorbar_box.x0 <= main_box.x1:
            raise RuntimeError("Right colorbar is not outside the 3D axes")


def assert_colorbar_alignment(
    main_ax,
    colorbar_ax,
    orientation: str,
    tolerance: float = 0.015,
    text_clearance_px: float = 10.0,
) -> None:
    """Verify colorbar centering, external placement, and axis-text clearance."""
    main_box = main_ax.get_position()
    colorbar_box = colorbar_ax.get_position()
    if orientation == "horizontal":
        error = abs(
            (main_box.x0 + main_box.x1) / 2
            - (colorbar_box.x0 + colorbar_box.x1) / 2
        )
        outside = colorbar_box.y1 < main_box.y0
    elif orientation == "vertical":
        error = abs(
            (main_box.y0 + main_box.y1) / 2
            - (colorbar_box.y0 + colorbar_box.y1) / 2
        )
        outside = colorbar_box.x0 > main_box.x1
        fig = main_ax.figure
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        axes = [main_ax.xaxis, main_ax.yaxis]
        if hasattr(main_ax, "zaxis"):
            axes.append(main_ax.zaxis)
        text_artists = [axis.label for axis in axes]
        for axis in axes:
            for tick in [*axis.get_major_ticks(), *axis.get_minor_ticks()]:
                text_artists.extend((tick.label1, tick.label2))
        right_edges = []
        for text in text_artists:
            if not text.get_visible() or not text.get_text().strip():
                continue
            text_box = text.get_window_extent(renderer=renderer)
            if 0 <= text_box.x1 <= fig.bbox.x1 * 1.5:
                right_edges.append(text_box.x1)
        colorbar_left_px = colorbar_box.x0 * fig.bbox.width
        text_clear = (
            not right_edges
            or colorbar_left_px - max(right_edges) >= text_clearance_px
        )
    else:
        raise ValueError("orientation must be 'horizontal' or 'vertical'")
    if orientation == "horizontal":
        text_clear = True
    if not outside or not text_clear or error > tolerance:
        raise RuntimeError(
            f"Misaligned {orientation} colorbar: outside={outside}, "
            f"text_clear={text_clear}, error={error:.4f}"
        )


def save_publication_bundle(
    fig,
    out_base: str,
    dpi: int = 600,
    check_layout: bool = True,
    extra_formats: bool = False,
) -> None:
    """Save the default PDF/PNG bundle, optionally with heavier edit formats."""
    assert_manuscript_fonts(fig)
    if check_layout:
        assert_no_text_collisions(fig)
        assert_no_legend_reference_line_collisions(fig)
    out_path = Path(out_base)
    fig.savefig(out_path.with_suffix(".pdf"))
    fig.savefig(out_path.with_suffix(".png"), dpi=200)
    if extra_formats:
        fig.savefig(out_path.with_suffix(".svg"))
        fig.savefig(out_path.with_suffix(".tiff"), dpi=dpi)
