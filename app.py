import streamlit as st
import matplotlib.pyplot as plt
import math
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages

VERSION = "1.4"

# -----------------------------
# Streamlit UI
# -----------------------------
st.set_page_config(layout="wide")
st.title(f"Розрахунок смуг для резервуара — v{VERSION}")

# -----------------------------
# Helpers
# -----------------------------
def clamp(x, a, b):
    return max(a, min(b, x))

def primitive_F(R: float, y: float) -> float:
    """
    Первісна для 2*sqrt(R^2 - y^2):
    F(y) = y*sqrt(R^2-y^2) + R^2*asin(y/R)
    """
    y = clamp(y, -R, R)
    inside = max(0.0, R * R - y * y)
    return y * math.sqrt(inside) + R * R * math.asin(clamp(y / R, -1.0, 1.0))

def chord_len(R: float, y: float) -> float:
    """Довжина хорди кола радіуса R на рівні y."""
    if y <= -R or y >= R:
        return 0.0
    return 2.0 * math.sqrt(max(0.0, R * R - y * y))

def build_patterns(Wrem: float):
    """
    Повертає список патернів (рядків), де кожен рядок — список сегментів по ширині (м).
    - Для Wrem 4..6 — “шахматка” (як було)
    - Для довільних значень — сегментація на 3 м + залишок
    """
    eps = 1e-6
    if abs(Wrem - 4.0) < eps:
        return [[3.0, 1.0], [1.0, 3.0]]
    if abs(Wrem - 5.0) < eps:
        return [[3.0, 2.0], [2.0, 3.0]]
    if abs(Wrem - 6.0) < eps:
        return [[1.0, 3.0, 2.0], [2.0, 3.0, 1.0]]

    # Довільний Wrem: ріжемо на 3.0 м + залишок
    segs = []
    remain = Wrem
    while remain > 3.0 + 1e-9:
        segs.append(3.0)
        remain -= 3.0
    if remain > 1e-9:
        segs.append(round(remain, 2))
    if not segs:
        segs = [round(Wrem, 2)]
    return [segs]  # без шахматки

def fmt_m(x):
    return f"{x:.2f} м"

def fmt_m2(x):
    return f"{x:.3f} м²"

def fmt_kg(x):
    return f"{x:.2f} кг"

def make_table_figure(title: str, rows, col_labels=None, fontsize=10):
    """
    rows: list[list[str]]
    col_labels: list[str] or None
    """
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    fig.clf()
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.01, 0.98, title, fontsize=14, fontweight="bold", va="top")

    table_ax = fig.add_axes([0.05, 0.05, 0.90, 0.88])
    table_ax.axis("off")

    tbl = table_ax.table(
        cellText=rows,
        colLabels=col_labels,
        cellLoc="left",
        loc="upper left"
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.scale(1.0, 1.25)

    # легкий стиль
    if col_labels:
        for (r, c), cell in tbl.get_celld().items():
            if r == 0:
                cell.set_text_props(weight="bold")
                cell.set_linewidth(1.2)
            else:
                cell.set_linewidth(0.6)
    else:
        for (r, c), cell in tbl.get_celld().items():
            cell.set_linewidth(0.6)

    return fig

# -----------------------------
# Inputs (main)
# -----------------------------
col1, col2, col3, col4 = st.columns(4)
with col1:
    D = st.number_input("Діаметр резервуара, м", value=2.5, min_value=0.1, step=0.1, format="%.1f")
with col2:
    L = st.number_input("Довжина резервуара, м", value=4.3, min_value=0.1, step=0.1, format="%.1f")
with col3:
    Wrem = st.number_input("Ширина ремонтованої ділянки (по колу), м", value=1.0, min_value=0.1, step=0.1, format="%.1f")
with col4:
    overlap_cm = st.number_input("Нахлест, см", value=5.0, min_value=0.0, step=0.5, format="%.1f")

st.markdown("### Зварювання та електроди (для звіту)")
e1, e2, e3, e4 = st.columns(4)
with e1:
    electrode_d_mm = st.number_input("Діаметр електрода, мм", value=3.2, min_value=2.0, step=0.1, format="%.1f")
with e2:
    a_mm = st.number_input("Катет кутового шва, мм", value=4.0, min_value=2.0, step=0.5, format="%.1f")
with e3:
    eta = st.number_input("ККД наплавлення (η) для MMA", value=0.65, min_value=0.40, max_value=0.90, step=0.01, format="%.2f")
with e4:
    k_loss = st.number_input("Коеф. втрат (K)", value=1.10, min_value=1.00, max_value=1.30, step=0.01, format="%.2f")

c5, c6, c7 = st.columns(3)
with c5:
    pack_mass_kg = st.number_input("Маса 1 пачки електродів, кг", value=2.5, min_value=0.5, step=0.5, format="%.1f")
with c6:
    include_shell_to_bottom = st.checkbox("Додати шви приварки циліндричної частини до днищ (2×довжина кола)", value=True)
with c7:
    double_overlap_weld = st.checkbox("Нахлест проварюється двома лініями (технологічно, для електродів)", value=True)

# -----------------------------
# Calculation
# -----------------------------
if st.button("Розрахувати"):
    if D <= 0 or L <= 0 or Wrem <= 0 or overlap_cm < 0:
        st.error("Введіть додатні D, L, Wrem та невід'ємний нахлест.")
        st.stop()

    overlap = overlap_cm / 100.0  # m
    R = D / 2.0
    h_smuha = 0.5
    circumference = 2.0 * math.pi * R

    if Wrem > circumference + 1e-9:
        st.warning(f"Wrem = {Wrem:.2f} м більша за довжину кола {circumference:.2f} м. Розгортка буде умовною.")

    # ---- Hcrit (днище) ----
    alpha = (Wrem / 2.0) / R
    Hcrit = R * (1.0 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))

    # ============================================
    # 1) Bottom (one)
    # ============================================
    widths_bot = []
    used_heights_bot = []
    areas_bot = []
    smuhaDict = {}

    fig_bot, ax_bot = plt.subplots(figsize=(6, 6))
    ax_bot.set_aspect("equal")
    ax_bot.set_title("Днище резервуара (синій — смуга, червоний — зона нахлесту)", fontsize=13, fontweight="bold")

    circle = plt.Circle((0, 0), R, edgecolor="black", facecolor="lightyellow", alpha=0.3, zorder=0)
    ax_bot.add_patch(circle)

    ax_bot.axvline(0, color="red", linestyle="--", linewidth=1, alpha=0.7)
    ax_bot.text(
        R * 0.60, -R + Hcrit + 0.03,
        f"H = {Hcrit:.3f} м",
        fontsize=10,
        bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.2"),
        zorder=10
    )

    y_bot_global = -R
    y_cut = -R + Hcrit

    for j in range(n_bot):
        y_bot = y_bot_global + j * h_smuha
        y_top_nom = y_bot + h_smuha

        # real height for last strip
        real_h = h_smuha
        if j == n_bot - 1:
            real_h = max(0.0, min(y_cut - y_bot, h_smuha))

        # width reference y
        if y_top_nom <= 0:
            y_ref = y_top_nom
        elif y_bot >= 0:
            y_ref = y_bot
        else:
            y_ref = 0.0

        if abs(y_ref) >= R:
            width = 0.0
        else:
            width = 2.0 * math.sqrt(max(0.0, R * R - y_ref * y_ref))

        widths_bot.append(width)
        used_heights_bot.append(real_h)

        y_top_clip = min(y_bot + real_h, y_cut)
        area_strip_exact = primitive_F(R, y_top_clip) - primitive_F(R, y_bot)
        areas_bot.append(area_strip_exact)

        # strip accounting (2 symmetric)
        width_rnd = round(width, 2)
        key_bot = f"{width_rnd:>5.2f}м x {h_smuha:>3.2f}м"
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2

        x_left = -width / 2.0

        ax_bot.add_patch(plt.Rectangle((x_left, y_bot), width, h_smuha,
                                       edgecolor="black", facecolor="skyblue", alpha=0.45, zorder=2))

        ov_h = h_smuha + (overlap / 2.0 if j < n_bot - 1 else 0.0)
        ax_bot.add_patch(plt.Rectangle((x_left, y_bot), width, ov_h,
                                       edgecolor="red", facecolor="none",
                                       linestyle="--", linewidth=1, alpha=0.7, zorder=3))

        if real_h > 0:
            ax_bot.add_patch(plt.Rectangle((x_left, y_bot), width, real_h,
                                           edgecolor="black", facecolor="skyblue", alpha=0.55, zorder=4))

        extra_h = h_smuha - real_h
        if extra_h > 1e-9:
            ax_bot.add_patch(plt.Rectangle((x_left, y_bot + real_h), width, extra_h,
                                           edgecolor="red", facecolor="none", hatch="///", alpha=0.5, zorder=5))

        ax_bot.text(0, y_bot + h_smuha / 2.0, f"S{j+1}\n{areas_bot[-1]:.3f} м²",
                    ha="center", va="center", fontsize=8, zorder=7)

    cum_area_bot = sum(areas_bot)                 # one bottom
    total_area_both_bottoms = 2.0 * cum_area_bot  # two bottoms

    margin = R * 0.1
    ax_bot.set_xlim(-R - margin, R + margin)
    ax_bot.set_ylim(-R - margin, R + margin)
    ax_bot.set_xlabel("x (м)")
    ax_bot.set_ylabel("y (м)")
    ax_bot.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # ============================================
    # 2) Cylinder (unfolded)
    # ============================================
    # number of rows along length with vertical overlap
    full_rows = 1
    covered_height = h_smuha
    while covered_height + 1e-9 < L:
        full_rows += 1
        covered_height += (h_smuha - overlap)

    patterns = build_patterns(Wrem)

    fig_cyl, ax_cyl = plt.subplots(figsize=(8, 6))
    ax_cyl.set_aspect("equal", adjustable="box")  # 1:1 scale
    ax_cyl.set_title("Розгорнута поверхня циліндра (чергування смуг, червоний — зона нахлесту)",
                     fontsize=13, fontweight="bold")

    # tank body background (full circumference)
    ax_cyl.add_patch(plt.Rectangle((-circumference/2.0, 0), circumference, L,
                                   edgecolor="black", facecolor="#d0d0d0", linewidth=1, zorder=1))
    ax_cyl.axvline(0, color="red", linestyle="--", linewidth=1, alpha=0.7)
    ax_cyl.axhline(L, color="red", linestyle="--", linewidth=1, alpha=0.7)

    x_start = -Wrem / 2.0
    площа_cyl_з_нахлестом = 0.0
    max_top_for_ylim = L

    for rowNum in range(full_rows):
        y_off = rowNum * (h_smuha - overlap)
        if y_off >= L + overlap/2.0:
            break

        pat = patterns[rowNum % len(patterns)]
        x_off = x_start

        visible_height = max(0.0, min(h_smuha, L - y_off))
        extra_down = min(overlap/2.0, max(0.0, y_off))
        is_last_row = (rowNum == full_rows - 1)

        top_ov_end = y_off + visible_height + (0.0 if is_last_row else overlap/2.0)
        ov_y_base = max(0.0, y_off - extra_down)
        ov_h = max(0.0, top_ov_end - ov_y_base)

        if is_last_row:
            max_top_for_ylim = max(max_top_for_ylim, y_off + h_smuha)
        else:
            max_top_for_ylim = max(max_top_for_ylim, top_ov_end)

        for seg_idx, seg in enumerate(pat):
            base_h = h_smuha if is_last_row else visible_height

            ax_cyl.add_patch(plt.Rectangle((x_off, y_off), seg, base_h,
                                           edgecolor="black",
                                           facecolor=("orange" if (rowNum % 2 == 0) else "lightgreen"),
                                           alpha=0.7, zorder=2))

            if is_last_row and (y_off + h_smuha > L):
                ax_cyl.add_patch(plt.Rectangle((x_off, L), seg, (y_off + h_smuha - L),
                                               edgecolor="red", facecolor="none", hatch="///", alpha=0.5, zorder=4))

            # side overlap only for 4..6 "checker"
            is_checker = abs(Wrem - 4.0) < 1e-6 or abs(Wrem - 5.0) < 1e-6 or abs(Wrem - 6.0) < 1e-6
            if is_checker and 0 < seg_idx < len(pat) - 1:
                ov_x = x_off - overlap/2.0
                ov_w = seg + overlap
            else:
                ov_x = x_off
                ov_w = seg

            if ov_h > 0:
                ax_cyl.add_patch(plt.Rectangle((ov_x, ov_y_base), ov_w, ov_h,
                                               edgecolor="red", facecolor="none",
                                               linestyle="--", linewidth=1, alpha=0.7, zorder=3))
            площа_cyl_з_нахлестом += ov_w * ov_h

            ax_cyl.text(x_off + seg/2.0, y_off + (base_h/2.0 if base_h > 0 else 0.02),
                        f"{seg:>5.2f}м", ha="center", va="center", fontsize=8, zorder=5)

            key_cyl = f"{round(seg, 2):>5.2f}м x {h_smuha:>3.2f}м"
            smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1

            x_off += seg

    ax_cyl.set_xlim(-circumference/2.0, circumference/2.0)
    ax_cyl.set_ylim(0, max_top_for_ylim)
    ax_cyl.set_xlabel("довжина поверхні (м)")
    ax_cyl.set_ylabel("висота (м)")
    ax_cyl.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # ============================================
    # 3) Areas + sheets
    # ============================================
    площа_cyl = Wrem * L
    реально_ремонт_площа = total_area_both_bottoms + площа_cyl

    площа_днищ_з_нахлестом_факт = 2.0 * sum(
        w * (h + (overlap/2.0 if idx < n_bot - 1 else 0.0))
        for idx, (w, h) in enumerate(zip(widths_bot, used_heights_bot))
    )
    фактична_площа_матеріалів = площа_днищ_з_нахлестом_факт + площа_cyl_з_нахлестом

    sheet_area = 6.0 * 1.5
    total_sheets = math.ceil(фактична_площа_матеріалів / sheet_area)

    # ============================================
    # 4) Welds (GEOMETRY) + weld of shell-to-bottom
    # ============================================
    # Cylinder repair welds (geometry, your previous method)
    weld_h_cyl = (full_rows + 1) * Wrem

    def row_visible_height(i: int) -> float:
        y = i * (h_smuha - overlap)
        if y >= L:
            return 0.0
        return min(h_smuha, L - y)

    weld_v_cyl = 0.0
    for i in range(full_rows):
        h_eff = row_visible_height(i)
        if h_eff <= 0:
            break
        pat = patterns[i % len(patterns)]
        n_vertical = (len(pat) - 1) + 2  # inside + 2 edges
        weld_v_cyl += n_vertical * h_eff

    weld_cyl_repair_geom = weld_h_cyl + weld_v_cyl

    # Bottom repair welds (geometry)
    y0 = -R
    y_top = -R + Hcrit

    levels = []
    k = 0
    while True:
        y_level = y0 + k * h_smuha
        if y_level >= y_top - 1e-9:
            break
        levels.append(y_level)
        k += 1
    levels.append(y_top)

    weld_h_bot_one = sum(chord_len(R, y) for y in levels)

    theta_top = math.asin(clamp(y_top / R, -1.0, 1.0))
    arc_len_one_side = R * (theta_top - (-math.pi / 2.0))
    weld_side_bot_one = 2.0 * arc_len_one_side

    weld_bot_one_geom = weld_h_bot_one + weld_side_bot_one
    weld_both_bottoms_geom = 2.0 * weld_bot_one_geom

    # Shell-to-bottom factory/constructive seams (2 circumferential seams)
    weld_shell_to_bottom = (2.0 * circumference) if include_shell_to_bottom else 0.0

    total_weld_geom = weld_cyl_repair_geom + weld_both_bottoms_geom + weld_shell_to_bottom

    # ============================================
    # 5) Technological "bead length" for electrodes (overlap double line)
    #    Idea: geometry length + extra lines due to double overlap weld.
    # ============================================
    # Cylinder: between rows there are (full_rows-1) overlaps; geometry counts 1 line each.
    # Double overlap => +1 extra line per overlap => add (full_rows-1)*Wrem
    extra_cyl_overlap = 0.0
    if double_overlap_weld and full_rows > 1:
        extra_cyl_overlap = (full_rows - 1) * Wrem

    # Bottom: between strips there are (n_bot-1) internal levels; each has chord_len at that y.
    # Geometry already includes these once; double overlap => add them once more.
    extra_bot_overlap_one = 0.0
    if double_overlap_weld and n_bot > 1:
        # internal levels only (exclude first y0 and final y_top)
        internal_levels = []
        for idx in range(1, len(levels)-1):
            internal_levels.append(levels[idx])
        extra_bot_overlap_one = sum(chord_len(R, y) for y in internal_levels)
    extra_bot_overlap_both = 2.0 * extra_bot_overlap_one

    total_bead_length_tech = total_weld_geom + extra_cyl_overlap + extra_bot_overlap_both

    # ============================================
    # 6) Electrode consumption (professional, by deposited metal mass)
    # ============================================
    # Fillet weld metal area approx: A = a^2/2 (mm^2)
    # Volume per meter: A*1e-6 (m^3)
    # Deposit mass per meter: rho * volume
    rho = 7850.0  # kg/m^3 steel
    area_mm2 = (a_mm * a_mm) / 2.0
    deposit_mass_per_m = rho * area_mm2 * 1e-6  # kg/m (deposited metal)

    # Electrode mass per meter: deposit / eta * K
    electrode_mass_per_m = (deposit_mass_per_m / max(eta, 1e-9)) * k_loss

    electrodes_total_kg = total_bead_length_tech * electrode_mass_per_m
    packs_needed = math.ceil(electrodes_total_kg / max(pack_mass_kg, 1e-9))

    # “Reference phrase” example: on 20 m
    ref_L = 20.0
    ref_kg = ref_L * electrode_mass_per_m
    ref_packs = math.ceil(ref_kg / max(pack_mass_kg, 1e-9))

    # ============================================
    # 7) Strip listing
    # ============================================
    items = sorted(smuhaDict.items())

    # ============================================
    # Output (screen)
    # ============================================
    st.markdown("## Результати розрахунків")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("### Графіки")
        st.pyplot(fig_bot, clear_figure=False)
        st.pyplot(fig_cyl, clear_figure=False)

    with c2:
        st.markdown("### Професійний підсумок (для звіту клієнту)")
        st.markdown("#### Площі та матеріали")
        st.text(f"Площа одного днища: {fmt_m2(cum_area_bot)}")
        st.text(f"Площа обох днищ:    {fmt_m2(total_area_both_bottoms)}")
        st.text(f"Площа циліндричної частини: {fmt_m2(площа_cyl)}")
        st.markdown(f"<p style='color:blue;'>Ремонтована площа (фактична): {fmt_m2(реально_ремонт_площа)}</p>", unsafe_allow_html=True)
        st.markdown(f"<p style='color:red;'>Використана площа матеріалу (з нахлестом {overlap_cm:.1f} см): {fmt_m2(фактична_площа_матеріалів)}</p>", unsafe_allow_html=True)
        st.text(f"Листи 6×1,5 м (орієнтовно): {total_sheets} шт")

        st.markdown("#### Зварювання")
        st.text(f"Довжина швів (геометрична): {fmt_m(total_weld_geom)}")
        st.text(f"  • ремонт циліндра: {fmt_m(weld_cyl_repair_geom)}")
        st.text(f"  • ремонт днищ (2 шт): {fmt_m(weld_both_bottoms_geom)}")
        if include_shell_to_bottom:
            st.text(f"  • приварка циліндра до днищ (2×коло): {fmt_m(weld_shell_to_bottom)}")

        st.markdown("#### Електроди (розрахунок через масу наплавлення)")
        st.text(f"Довжина наплавлення (технологічна): {fmt_m(total_bead_length_tech)}")
        if double_overlap_weld:
            st.text(f"  • додатково по нахлестах: {fmt_m(extra_cyl_overlap + extra_bot_overlap_both)} (враховано)")
        st.text(f"Електроди: Ø {electrode_d_mm:.1f} мм | катет шва a = {a_mm:.1f} мм | η = {eta:.2f} | K = {k_loss:.2f}")
        st.markdown(f"<p style='color:green; font-weight:bold;'>Потрібно електродів: {fmt_kg(electrodes_total_kg)} ≈ {packs_needed} пач(ки) по {pack_mass_kg:.1f} кг</p>", unsafe_allow_html=True)

        st.markdown("#### Приклад формулювання для клієнта (довіра)")
        st.text(f"Наприклад, на {ref_L:.0f} м наплавлення потрібно {fmt_kg(ref_kg)} ≈ {ref_packs} пач(ки) електродів Ø {electrode_d_mm:.1f} мм (по {pack_mass_kg:.1f} кг).")

        st.markdown("#### Смуги матеріалу")
        st.text("Розмір смуги         К-сть")
        for key, cnt in items:
            st.text(f"{key:>14s}   {cnt:>3d} шт")
        st.text(f"Нахлест: {overlap_cm:.1f} см")

    # ============================================
    # PDF export (graphs + tables)
    # ============================================
    st.markdown("---")
    st.markdown("### Збереження результатів (PDF з таблицями)")

    # Build table pages
    inputs_rows = [
        ["Діаметр резервуара D", f"{D:.2f} м"],
        ["Довжина резервуара L", f"{L:.2f} м"],
        ["Ширина ремонтованої ділянки Wrem", f"{Wrem:.2f} м"],
        ["Нахлест", f"{overlap_cm:.1f} см"],
        ["Довжина кола", f"{circumference:.3f} м"],
        ["Висота смуги", f"{h_smuha:.2f} м"],
        ["К-сть рядів на циліндрі", f"{full_rows} шт"],
        ["К-сть смуг на днище (по висоті)", f"{n_bot} шт"],
        ["Діаметр електрода", f"{electrode_d_mm:.1f} мм"],
        ["Катет шва a", f"{a_mm:.1f} мм"],
        ["ККД наплавлення η", f"{eta:.2f}"],
        ["Коеф. втрат K", f"{k_loss:.2f}"],
        ["Маса 1 пачки", f"{pack_mass_kg:.1f} кг"],
        ["Подвійна лінія нахлесту (для електродів)", "Так" if double_overlap_weld else "Ні"],
        ["Додано шви приварки циліндра до днищ", "Так" if include_shell_to_bottom else "Ні"],
    ]

    areas_rows = [
        ["Площа одного днища", fmt_m2(cum_area_bot)],
        ["Площа обох днищ", fmt_m2(total_area_both_bottoms)],
        ["Площа циліндричної частини", fmt_m2(площа_cyl)],
        ["Ремонтована площа (фактична)", fmt_m2(реально_ремонт_площа)],
        [f"Використана площа матеріалу (нахлест {overlap_cm:.1f} см)", fmt_m2(фактична_площа_матеріалів)],
        ["Листи 6×1,5 м (орієнтовно)", f"{total_sheets} шт"],
    ]

    weld_rows = [
        ["Довжина швів (геометрична)", fmt_m(total_weld_geom)],
        ["  • ремонт циліндра", fmt_m(weld_cyl_repair_geom)],
        ["  • ремонт днищ (2 шт)", fmt_m(weld_both_bottoms_geom)],
    ]
    if include_shell_to_bottom:
        weld_rows.append(["  • приварка циліндра до днищ (2×коло)", fmt_m(weld_shell_to_bottom)])

    weld_rows += [
        ["Довжина наплавлення (технологічна)", fmt_m(total_bead_length_tech)],
    ]
    if double_overlap_weld:
        weld_rows.append(["  • додатково по нахлестах (враховано)", fmt_m(extra_cyl_overlap + extra_bot_overlap_both)])

    electrodes_rows = [
        ["Метод", "Розрахунок через масу наплавлення (катет шва, η, K)"],
        ["Електрод", f"Ø {electrode_d_mm:.1f} мм"],
        ["Катет шва", f"a = {a_mm:.1f} мм"],
        ["Витрата електродів на 1 м (оцінка)", fmt_kg(electrode_mass_per_m)],
        ["Потрібно електродів всього", fmt_kg(electrodes_total_kg)],
        ["Пачки електродів", f"{packs_needed} пач(ки) × {pack_mass_kg:.1f} кг"],
        ["Приклад для формулювання", f"На {ref_L:.0f} м — {fmt_kg(ref_kg)} ≈ {ref_packs} пач(ки)"],
    ]

    strips_rows = [["Розмір смуги", "К-сть, шт"]]
    for key, cnt in items:
        strips_rows.append([key, str(cnt)])

    pdf_buffer = BytesIO()
    with PdfPages(pdf_buffer) as pdf:
        # Graph pages
        pdf.savefig(fig_bot, bbox_inches="tight")
        pdf.savefig(fig_cyl, bbox_inches="tight")

        # Tables pages
        fig1 = make_table_figure("1) Вхідні дані", inputs_rows, col_labels=["Параметр", "Значення"], fontsize=10)
        pdf.savefig(fig1, bbox_inches="tight")

        fig2 = make_table_figure("2) Площі та матеріали", areas_rows, col_labels=["Показник", "Значення"], fontsize=11)
        pdf.savefig(fig2, bbox_inches="tight")

        fig3 = make_table_figure("3) Зварювання та електроди (прозорий розрахунок)", weld_rows + [["", "" ]] + electrodes_rows,
                                 col_labels=["Показник", "Значення"], fontsize=10)
        pdf.savefig(fig3, bbox_inches="tight")

        fig4 = make_table_figure("4) Відомість смуг матеріалу", strips_rows[1:], col_labels=strips_rows[0], fontsize=10)
        pdf.savefig(fig4, bbox_inches="tight")

        plt.close(fig1); plt.close(fig2); plt.close(fig3); plt.close(fig4)

    st.download_button(
        label="Зберегти як PDF (графіки + таблиці)",
        data=pdf_buffer.getvalue(),
        file_name="резервуар_звіт.pdf",
        mime="application/pdf"
    )

    # PNG export
    cc1, cc2 = st.columns(2)
    with cc1:
        png_buffer_bot = BytesIO()
        fig_bot.savefig(png_buffer_bot, format="png", dpi=200, bbox_inches="tight")
        st.download_button(
            label="Зберегти днище (PNG)",
            data=png_buffer_bot.getvalue(),
            file_name="днище.png",
            mime="image/png"
        )
    with cc2:
        png_buffer_cyl = BytesIO()
        fig_cyl.savefig(png_buffer_cyl, format="png", dpi=200, bbox_inches="tight")
        st.download_button(
            label="Зберегти циліндр (PNG)",
            data=png_buffer_cyl.getvalue(),
            file_name="циліндр.png",
            mime="image/png"
        )

    # Close main figs to avoid memory leaks
    plt.close(fig_bot)
    plt.close(fig_cyl)
