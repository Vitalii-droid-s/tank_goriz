import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
from datetime import date
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg

VERSION = "1.7"

# =============================
# Matplotlib fonts (стабільно в PDF)
# =============================
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False

# =============================
# Streamlit
# =============================
st.set_page_config(layout="wide")
st.title(f"Розрахунок ремонту резервуара — v{VERSION}")

# =============================
# Helpers
# =============================
def clamp(x, a, b):
    return max(a, min(b, x))

def primitive_F(R: float, y: float) -> float:
    y = clamp(y, -R, R)
    inside = max(0.0, R * R - y * y)
    return y * math.sqrt(inside) + R * R * math.asin(clamp(y / R, -1.0, 1.0))

def chord_len(R: float, y: float) -> float:
    if y <= -R or y >= R:
        return 0.0
    return 2.0 * math.sqrt(max(0.0, R * R - y * y))

def build_patterns(Wrem: float):
    """
    - Wrem 4..6: "шахматка"
    - інакше: сегментація на 3.0 м + залишок
    """
    eps = 1e-6
    if abs(Wrem - 4.0) < eps:
        return [[3.0, 1.0], [1.0, 3.0]]
    if abs(Wrem - 5.0) < eps:
        return [[3.0, 2.0], [2.0, 3.0]]
    if abs(Wrem - 6.0) < eps:
        return [[1.0, 3.0, 2.0], [2.0, 3.0, 1.0]]

    segs = []
    remain = Wrem
    while remain > 3.0 + 1e-9:
        segs.append(3.0)
        remain -= 3.0
    if remain > 1e-9:
        segs.append(round(remain, 2))
    return [segs]

def fig_to_img_array(fig, dpi=160):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return mpimg.imread(buf)

def make_table_ax(
    ax, title, rows,
    col_labels=("Параметр", "Значення"),
    font_size=10,
    title_pad=16,
    scale_y=1.25,
    bbox=(0.0, 0.02, 1.0, 0.86)
):
    ax.axis("off")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=title_pad)

    tbl = ax.table(
        cellText=rows,
        colLabels=list(col_labels),
        colLoc="center",
        cellLoc="left",
        loc="center",
        bbox=list(bbox)
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(font_size)
    tbl.scale(1.0, scale_y)

    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.8)
        if r == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#f2f2f2")
            cell.set_height(cell.get_height() * 1.10)

    return tbl

def fmt_m(x):  return f"{x:.2f} м"
def fmt_m2(x): return f"{x:.3f} м²"
def fmt_kg(x): return f"{x:.2f} кг"

# =============================
# Реквізити звіту
# =============================
meta_c1, meta_c2, meta_c3, meta_c4 = st.columns([2.2, 1.0, 1.0, 1.0])
with meta_c1:
    tank_type = st.text_input("Тип резервуара", value="Паливний резервуар")
with meta_c2:
    tank_id = st.text_input("ID резервуара", value="—")
with meta_c3:
    serial_no = st.text_input("Заводський №", value="—")
with meta_c4:
    repair_date = st.date_input("Дата ремонту", value=date.today())

st.markdown("---")

# =============================
# Режим роботи
# =============================
mode = st.radio(
    "Режим роботи",
    ["Ремонт резервуара (повний)", "Ремонт резервуара (латочний)"],
    horizontal=True
)

# =============================
# Спільні параметри електродів
# =============================
e1, e2, e3 = st.columns(3)
with e1:
    electrode_diam_mm = st.number_input("Діаметр електрода, мм", value=3.2, min_value=1.6, step=0.2, format="%.1f")
with e2:
    spec_consumption = st.number_input("Питома витрата електродів, кг/м наплавлення", value=0.25, min_value=0.01, step=0.01, format="%.2f")
with e3:
    pack_mass = st.number_input("Маса 1 пачки електродів, кг", value=2.5, min_value=0.5, step=0.5, format="%.1f")

# =============================
# Параметри для розрахунку наплавлення
# =============================
st.caption("Примітка: «Довжина наплавлення» = геометрична довжина швів × кількість проходів × кількість ліній шва (де задано).")
p1, p2 = st.columns([1, 2])
with p1:
    passes = st.number_input("К-сть проходів (для наплавлення)", value=1, min_value=1, step=1)
with p2:
    note_weld = st.checkbox("Показувати у звіті геометрію швів та наплавлення окремо", value=True)

st.markdown("---")

# =============================
# 1) ПОВНИЙ РЕМОНТ
# =============================
if mode == "Ремонт резервуара (повний)":
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        D = st.number_input("Діаметр резервуара, м", value=2.5, min_value=0.1, step=0.1, format="%.2f")
    with col2:
        L = st.number_input("Довжина резервуара, м", value=4.3, min_value=0.1, step=0.1, format="%.2f")
    with col3:
        Wrem = st.number_input("Ширина ремонтованої ділянки (по розгортці), м", value=1.0, min_value=0.1, step=0.1, format="%.2f")
    with col4:
        overlap_cm = st.number_input("Нахлест, см", value=5.0, min_value=0.0, step=0.5, format="%.1f")

    include_shell_to_bottom = st.checkbox(
        "Враховувати шви приєднання циліндричної частини до днищ",
        value=True
    )

    if st.button("Розрахувати"):
        if D <= 0 or L <= 0 or Wrem <= 0 or overlap_cm < 0:
            st.error("Введіть додатні D, L, Wrem та невід'ємний нахлест.")
            st.stop()

        overlap = overlap_cm / 100.0
        R = D / 2.0
        h_smuha = 0.5
        circumference = 2.0 * math.pi * R

        # ====== Днище: висота ремонту ======
        alpha = (Wrem / 2.0) / R
        Hcrit = R * (1.0 - math.cos(alpha))
        n_bot = max(1, math.ceil(Hcrit / h_smuha))
        y_cut = -R + Hcrit

        # =============================
        # ДНИЩЕ (графік)
        # =============================
        widths_bot, used_heights_bot, areas_bot = [], [], []
        smuhaDict = {}

        fig_bot, ax_bot = plt.subplots(figsize=(6.2, 6.2))
        ax_bot.set_aspect("equal")
        ax_bot.set_title("Днище резервуара (смуги + зона нахлесту)", fontsize=12, fontweight="bold")

        circle = plt.Circle((0, 0), R, edgecolor="black", facecolor="lightyellow", alpha=0.3, zorder=0)
        ax_bot.add_patch(circle)

        ax_bot.axvline(0, color="red", linestyle="--", linewidth=1, alpha=0.7)
        ax_bot.text(
            R * 0.55, -R + Hcrit + 0.03,
            f"H = {Hcrit:.3f} м",
            fontsize=9,
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.2"),
            zorder=10
        )

        y_bot_global = -R
        for j in range(n_bot):
            y_bot = y_bot_global + j * h_smuha
            y_top_nom = y_bot + h_smuha

            real_h = h_smuha
            if j == n_bot - 1:
                real_h = max(0.0, min(y_cut - y_bot, h_smuha))

            if y_top_nom <= 0:
                y_ref = y_top_nom
            elif y_bot >= 0:
                y_ref = y_bot
            else:
                y_ref = 0.0

            width = 0.0 if abs(y_ref) >= R else 2.0 * math.sqrt(max(0.0, R * R - y_ref * y_ref))
            widths_bot.append(width)
            used_heights_bot.append(real_h)

            y_top_clip = min(y_bot + real_h, y_cut)
            area_strip = primitive_F(R, y_top_clip) - primitive_F(R, y_bot)
            areas_bot.append(area_strip)

            key_bot = f"{round(width, 2):>5.2f}м × {h_smuha:>3.2f}м"
            smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2

            x_left = -width / 2.0
            ax_bot.add_patch(plt.Rectangle((x_left, y_bot), width, h_smuha,
                                           edgecolor="black", facecolor="skyblue", alpha=0.45, zorder=2))

            ov_h = h_smuha + (overlap / 2.0 if j < n_bot - 1 else 0.0)
            ax_bot.add_patch(plt.Rectangle((x_left, y_bot), width, ov_h,
                                           edgecolor="red", facecolor="none", linestyle="--",
                                           linewidth=1, alpha=0.7, zorder=3))

            if real_h > 0:
                ax_bot.add_patch(plt.Rectangle((x_left, y_bot), width, real_h,
                                               edgecolor="black", facecolor="skyblue", alpha=0.55, zorder=4))

            extra_h = h_smuha - real_h
            if extra_h > 1e-9:
                ax_bot.add_patch(plt.Rectangle((x_left, y_bot + real_h), width, extra_h,
                                               edgecolor="red", facecolor="none", hatch="///",
                                               alpha=0.5, zorder=5))

            ax_bot.text(0, y_bot + h_smuha / 2.0, f"S{j + 1}\n{area_strip:.3f} м²",
                        ha="center", va="center", fontsize=8, zorder=7)

        cum_area_bot = sum(areas_bot)
        total_area_both_bottoms = 2.0 * cum_area_bot

        m = R * 0.1
        ax_bot.set_xlim(-R - m, R + m)
        ax_bot.set_ylim(-R - m, R + m)
        ax_bot.set_xlabel("x (м)")
        ax_bot.set_ylabel("y (м)")
        ax_bot.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

        # =============================
        # ЦИЛІНДР (графік 1:1 по зоні ремонту)
        # =============================
        full_rows = 1
        covered_height = h_smuha
        while covered_height + 1e-9 < L:
            full_rows += 1
            covered_height += (h_smuha - overlap)

        patterns = build_patterns(Wrem)

        fig_cyl, ax_cyl = plt.subplots(figsize=(6.2, 6.2))
        ax_cyl.set_aspect("equal", adjustable="box")
        ax_cyl.set_title("Розгорнута поверхня (масштаб 1:1, зона ремонту)", fontsize=12, fontweight="bold")

        ax_cyl.add_patch(plt.Rectangle((-circumference / 2.0, 0), circumference, L,
                                       edgecolor="black", facecolor="#d0d0d0", linewidth=1, zorder=1))

        ax_cyl.axvline(0, color="red", linestyle="--", linewidth=1, alpha=0.7)
        ax_cyl.axhline(L, color="red", linestyle="--", linewidth=1, alpha=0.7)

        x_start = -Wrem / 2.0
        area_cyl_with_overlap = 0.0
        max_top_for_ylim = L

        for rowNum in range(full_rows):
            y_off = rowNum * (h_smuha - overlap)
            if y_off >= L + overlap / 2.0:
                break

            pat = patterns[rowNum % len(patterns)]
            x_off = x_start

            visible_height = max(0.0, min(h_smuha, L - y_off))
            extra_down = min(overlap / 2.0, max(0.0, y_off))
            is_last_row = (rowNum == full_rows - 1)

            top_ov_end = y_off + visible_height + (0.0 if is_last_row else overlap / 2.0)
            ov_y_base = max(0.0, y_off - extra_down)
            ov_h = max(0.0, top_ov_end - ov_y_base)
            max_top_for_ylim = max(max_top_for_ylim, (y_off + h_smuha) if is_last_row else top_ov_end)

            for seg_idx, seg in enumerate(pat):
                base_h = h_smuha if is_last_row else visible_height

                ax_cyl.add_patch(plt.Rectangle((x_off, y_off), seg, base_h,
                                               edgecolor="black",
                                               facecolor=("orange" if (rowNum % 2 == 0) else "lightgreen"),
                                               alpha=0.75, zorder=2))

                # боковий нахлест тільки для шахматки 4..6
                if abs(Wrem - 4.0) < 1e-6 or abs(Wrem - 5.0) < 1e-6 or abs(Wrem - 6.0) < 1e-6:
                    if 0 < seg_idx < len(pat) - 1:
                        ov_x = x_off - overlap / 2.0
                        ov_w = seg + overlap
                    else:
                        ov_x = x_off
                        ov_w = seg
                else:
                    ov_x = x_off
                    ov_w = seg

                if ov_h > 0:
                    ax_cyl.add_patch(plt.Rectangle((ov_x, ov_y_base), ov_w, ov_h,
                                                   edgecolor="red", facecolor="none",
                                                   linestyle="--", linewidth=1, alpha=0.7, zorder=3))

                area_cyl_with_overlap += ov_w * ov_h

                ax_cyl.text(x_off + seg / 2.0, y_off + base_h / 2.0, f"{seg:.2f} м",
                            ha="center", va="center", fontsize=9, zorder=5)

                key_cyl = f"{round(seg, 2):>5.2f}м × {h_smuha:>3.2f}м"
                smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1

                x_off += seg

        ax_cyl.set_xlim(-Wrem / 2.0 - 0.6, Wrem / 2.0 + 0.6)
        ax_cyl.set_ylim(0, max_top_for_ylim)
        ax_cyl.set_xlabel("довжина поверхні (м)")
        ax_cyl.set_ylabel("висота (м)")
        ax_cyl.grid(True, linestyle="--", linewidth=0.5, alpha=0.25)

        # =============================
        # Площі
        # =============================
        area_cyl = Wrem * L
        real_repair_area = total_area_both_bottoms + area_cyl

        area_bottoms_with_overlap_fact = 2.0 * sum(
            w * (h + (overlap / 2.0 if idx < n_bot - 1 else 0.0))
            for idx, (w, h) in enumerate(zip(widths_bot, used_heights_bot))
        )
        used_material_area = area_bottoms_with_overlap_fact + area_cyl_with_overlap
        sheet_area = 6.0 * 1.5
        total_sheets = math.ceil(used_material_area / sheet_area)

        # =============================
        # ШВИ (геометрія)
        # =============================
        weld_cyl_h = (full_rows + 1) * Wrem

        def row_visible_height(i: int) -> float:
            y = i * (h_smuha - overlap)
            if y >= L:
                return 0.0
            return min(h_smuha, L - y)

        weld_cyl_v = 0.0
        for i in range(full_rows):
            h_eff = row_visible_height(i)
            if h_eff <= 0:
                break
            pat = patterns[i % len(patterns)]
            n_vertical = (len(pat) - 1) + 2
            weld_cyl_v += n_vertical * h_eff

        weld_cyl_total = weld_cyl_h + weld_cyl_v

        # днище 1 шт
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

        weld_bot_h_one = sum(chord_len(R, y) for y in levels)

        theta_top = math.asin(clamp(y_top / R, -1.0, 1.0))
        arc_len_one_side = R * (theta_top - (-math.pi / 2.0))
        weld_bot_side_one = 2.0 * arc_len_one_side

        weld_bot_one = weld_bot_h_one + weld_bot_side_one
        weld_bottoms_total = 2.0 * weld_bot_one

        weld_shell_to_bottom = 0.0
        if include_shell_to_bottom:
            weld_shell_to_bottom = 2.0 * circumference

        weld_geom_total = weld_cyl_total + weld_bottoms_total + weld_shell_to_bottom
        weld_deposition_total = weld_geom_total * passes

        # =============================
        # Електроди
        # =============================
        electrodes_mass = weld_deposition_total * spec_consumption
        packs_needed = math.ceil(electrodes_mass / pack_mass)

        # =============================
        # WEB ВИВІД
        # =============================
        left_lines = [
            f"Тип резервуара: {tank_type}",
            f"ID резервуара: {tank_id}",
            f"Заводський №: {serial_no}",
            f"Дата ремонту: {repair_date.strftime('%d.%m.%Y')}",
            f"D: {D:.2f} м  |  R: {R:.2f} м",
            f"L: {L:.2f} м",
            f"Wrem: {Wrem:.2f} м",
            f"Нахлест: {overlap_cm:.1f} см",
            f"Висота ремонту днища H: {Hcrit:.3f} м",
            f"К-сть рядів на циліндрі: {full_rows} шт",
            "",
            f"Площа 1 днища: {cum_area_bot:.3f} м²",
            f"Площа 2 днищ: {total_area_both_bottoms:.3f} м²",
            f"Площа ремонту циліндра: {area_cyl:.3f} м²",
            f"Реально ремонтована площа: {real_repair_area:.3f} м²",
            f"Фактична площа матеріалу (з нахлестом): {used_material_area:.3f} м²",
            f"Оцінка листів 6×1.5 м: {total_sheets} шт",
            "",
        ]

        if note_weld:
            left_lines += [
                f"Шви (геометрія): {weld_geom_total:.2f} м",
                f"Наплавлення (проходи = {passes}): {weld_deposition_total:.2f} м",
                f"  • циліндр: {weld_cyl_total:.2f} м; днища: {weld_bottoms_total:.2f} м; приєднання: {weld_shell_to_bottom:.2f} м",
                "",
            ]
        else:
            left_lines += [
                f"Загальна довжина швів: {weld_geom_total:.2f} м",
                "",
            ]

        left_lines += [
            f"Електроди Ø{electrode_diam_mm:.1f} мм: {electrodes_mass:.2f} кг  (~ {packs_needed} пач. по {pack_mass:.2f} кг)",
        ]

        strips_items = sorted(smuhaDict.items(), key=lambda x: x[0])
        right_lines = ["Розмір смуги           К-сть"]
        for k, v in strips_items:
            right_lines.append(f"{k:>16s}   {v:>3d} шт")

        st.markdown("## Результати розрахунків")
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("### Графіки")
            st.pyplot(fig_bot, clear_figure=False)
            st.pyplot(fig_cyl, clear_figure=False)

        with c2:
            st.markdown("### Підсумкові дані")
            for line in left_lines:
                if line.startswith("Фактична площа"):
                    st.markdown(f"<p style='color:#b00020; font-weight:600;'>{line}</p>", unsafe_allow_html=True)
                elif line.startswith("Наплавлення"):
                    st.markdown(f"<p style='color:#006400; font-weight:700;'>{line}</p>", unsafe_allow_html=True)
                elif line.startswith("Електроди"):
                    st.markdown(f"<p style='color:#0033aa; font-weight:700;'>{line}</p>", unsafe_allow_html=True)
                else:
                    st.text(line)

            st.markdown("#### Відомість смуг")
            for line in right_lines:
                st.text(line)

        # =============================
        # PDF (3 сторінки, landscape)
        # =============================
        pdf_buffer = BytesIO()
        with PdfPages(pdf_buffer) as pdf:
            # Page 1: мета + 2 графіки
            fig1 = plt.figure(figsize=(11.69, 8.27))
            fig1.suptitle("ЗВІТ РОЗРАХУНКУ РЕМОНТУ (ПОВНИЙ)", fontsize=16, fontweight="bold", y=0.965)

            gs1 = gridspec.GridSpec(
                2, 2, figure=fig1,
                height_ratios=[0.33, 0.67],
                wspace=0.10, hspace=0.20,
                left=0.04, right=0.98, top=0.90, bottom=0.06
            )
            ax_meta = fig1.add_subplot(gs1[0, :])
            ax_a = fig1.add_subplot(gs1[1, 0])
            ax_b = fig1.add_subplot(gs1[1, 1])

            ax_meta.axis("off")
            meta_lines = [
                f"Тип: {tank_type}    |    ID: {tank_id}    |    Зав. №: {serial_no}    |    Дата: {repair_date.strftime('%d.%m.%Y')}",
                f"D = {D:.2f} м;  L = {L:.2f} м;  Wrem = {Wrem:.2f} м;  Нахлест = {overlap_cm:.1f} см",
                f"Шви (геометрія): {weld_geom_total:.2f} м    |    Наплавлення (проходи={passes}): {weld_deposition_total:.2f} м",
                f"Матеріал (факт): {used_material_area:.3f} м²    |    Листи 6×1.5 м: {total_sheets} шт",
                f"Електроди Ø{electrode_diam_mm:.1f} мм: {electrodes_mass:.2f} кг  (~ {packs_needed} × {pack_mass:.2f} кг; {spec_consumption:.2f} кг/м)",
            ]
            ax_meta.text(0.01, 0.92, meta_lines[0], fontsize=11, fontweight="bold", va="top")
            ax_meta.text(0.01, 0.66, meta_lines[1], fontsize=10, va="top")
            ax_meta.text(0.01, 0.41, meta_lines[2], fontsize=10, va="top")
            ax_meta.text(0.01, 0.22, meta_lines[3], fontsize=10, va="top")
            ax_meta.text(0.01, 0.05, meta_lines[4], fontsize=10, va="top")

            img_bot = fig_to_img_array(fig_bot, dpi=160)
            img_cyl = fig_to_img_array(fig_cyl, dpi=160)
            ax_a.imshow(img_bot); ax_b.imshow(img_cyl)
            ax_a.set_title("Рис. 1 — Схема смуг на днищі", fontsize=11, fontweight="bold", pad=8)
            ax_b.set_title("Рис. 2 — Розгортка циліндра (масштаб 1:1)", fontsize=11, fontweight="bold", pad=8)
            ax_a.axis("off"); ax_b.axis("off")

            pdf.savefig(fig1)
            plt.close(fig1)

            # Page 2: Ремонт + Шви
            fig2 = plt.figure(figsize=(11.69, 8.27))
            fig2.suptitle("Розрахунок ремонту та зварних швів", fontsize=16, fontweight="bold", y=0.965)

            gs2 = gridspec.GridSpec(1, 2, figure=fig2, wspace=0.12, left=0.04, right=0.98, top=0.90, bottom=0.08)
            ax21 = fig2.add_subplot(gs2[0, 0])
            ax22 = fig2.add_subplot(gs2[0, 1])

            repair_rows = [
                ["Тип резервуара", tank_type],
                ["ID резервуара", tank_id],
                ["Заводський №", serial_no],
                ["Дата ремонту", repair_date.strftime("%d.%m.%Y")],
                ["Діаметр D", fmt_m(D)],
                ["Радіус R", fmt_m(R)],
                ["Довжина L", fmt_m(L)],
                ["Ширина ремонту Wrem", fmt_m(Wrem)],
                ["Нахлест", f"{overlap_cm:.1f} см"],
                ["Висота ремонту днища H", fmt_m(Hcrit)],
                ["Кількість рядів на циліндрі", f"{full_rows} шт"],
            ]

            weld_rows = [
                ["Циліндр: горизонтальні", fmt_m(weld_cyl_h)],
                ["Циліндр: вертикальні", fmt_m(weld_cyl_v)],
                ["Циліндр: разом", fmt_m(weld_cyl_total)],
                ["Днища: разом", fmt_m(weld_bottoms_total)],
                ["Приєднання днищ до циліндра", fmt_m(weld_shell_to_bottom)],
                ["Шви (геометрія) — разом", fmt_m(weld_geom_total)],
                ["Проходи", f"{passes}"],
                ["Наплавлення — разом", fmt_m(weld_deposition_total)],
            ]

            make_table_ax(ax21, "1) Дані по ремонту", repair_rows, font_size=10, title_pad=16)
            make_table_ax(ax22, "2) Довжина швів / наплавлення", weld_rows, font_size=10, title_pad=16)

            pdf.savefig(fig2)
            plt.close(fig2)

            # Page 3: Матеріали + Електроди + Смуги
            fig3 = plt.figure(figsize=(11.69, 8.27))
            fig3.suptitle("Матеріали та електроди", fontsize=16, fontweight="bold", y=0.965)

            gs3 = gridspec.GridSpec(
                2, 2, figure=fig3,
                height_ratios=[1.0, 1.05],
                wspace=0.12, hspace=0.28,
                left=0.04, right=0.98, top=0.90, bottom=0.06
            )
            ax31 = fig3.add_subplot(gs3[0, 0])
            ax32 = fig3.add_subplot(gs3[0, 1])
            ax33 = fig3.add_subplot(gs3[1, :])

            areas_rows = [
                ["Площа 1 днища", fmt_m2(cum_area_bot)],
                ["Площа 2 днищ", fmt_m2(total_area_both_bottoms)],
                ["Площа ремонту циліндра", fmt_m2(area_cyl)],
                ["Реально ремонтована площа", fmt_m2(real_repair_area)],
                ["Фактична площа матеріалу (з нахлестом)", fmt_m2(used_material_area)],
                ["Листи 6×1.5 м (оцінка)", f"{total_sheets} шт"],
            ]

            electrodes_rows = [
                ["Діаметр електрода", f"{electrode_diam_mm:.1f} мм"],
                ["Питома витрата", f"{spec_consumption:.2f} кг/м"],
                ["Наплавлення (разом)", fmt_m(weld_deposition_total)],
                ["Потрібна маса електродів", fmt_kg(electrodes_mass)],
                ["Маса 1 пачки", fmt_kg(pack_mass)],
                ["Кількість пачок", f"{packs_needed} шт"],
            ]

            strips_rows = [[k.replace("×", "x"), f"{v} шт"] for k, v in strips_items]

            make_table_ax(ax31, "3) Площі та матеріали", areas_rows, font_size=10, title_pad=16)
            make_table_ax(ax32, "4) Електроди (оцінка)", electrodes_rows, font_size=10, title_pad=16)
            make_table_ax(
                ax33, "5) Відомість смуг (розкрій/кількість)",
                strips_rows,
                col_labels=("Розмір смуги", "К-сть"),
                font_size=9,
                title_pad=14,
                scale_y=1.20
            )

            pdf.savefig(fig3)
            plt.close(fig3)

        st.markdown("---")
        st.markdown("### Збереження результатів")
        st.download_button(
            label="Зберегти звіт як PDF",
            data=pdf_buffer.getvalue(),
            file_name="резервуар_звіт.pdf",
            mime="application/pdf"
        )

        cpng1, cpng2 = st.columns(2)
        with cpng1:
            buf = BytesIO()
            fig_bot.savefig(buf, format="png", dpi=200, bbox_inches="tight")
            st.download_button("Зберегти днище (PNG)", buf.getvalue(), "днище.png", "image/png")
        with cpng2:
            buf = BytesIO()
            fig_cyl.savefig(buf, format="png", dpi=200, bbox_inches="tight")
            st.download_button("Зберегти циліндр (PNG)", buf.getvalue(), "циліндр.png", "image/png")

        plt.close(fig_bot)
        plt.close(fig_cyl)

# =============================
# 2) ЛАТОЧНИЙ РЕМОНТ
# =============================
else:
    st.subheader("Латочний ремонт (прямокутні латки)")

    # ✅ фіксована ширина смуги як у повному ремонті: 500 мм
    STRIP_W_M = 0.5

    a1, a2, a3, a4 = st.columns(4)
    with a1:
        cut_dir = st.selectbox(
            "Різати латку смугами 0,5 м по",
            ["Висоті H", "Ширині W"],  # обираєш, як тобі потрібно технологічно
            index=0
        )
    with a2:
        overlap_cm = st.number_input("Нахлест між смугами, см", value=5.0, min_value=0.0, step=0.5, format="%.1f")
    with a3:
        seam_lines_on_lap = st.number_input("Ліній шва на нахлесті (1 або 2)", value=2, min_value=1, max_value=2, step=1)
    with a4:
        seam_lines_perimeter = st.number_input("Ліній шва по периметру латки (1 або 2)", value=1, min_value=1, max_value=2, step=1)

    b1, b2 = st.columns([1, 2])
    with b1:
        patches_count = st.number_input("Кількість латок", value=1, min_value=1, step=1)
    with b2:
        allowance_mm = st.number_input("Припуск по периметру (додається з кожного боку), мм", value=0.0, min_value=0.0, step=5.0, format="%.0f")

    st.caption("Розмір латки вводиться у мм. Латку розкроює на смуги шириною 0,5 м (500 мм) по вибраній стороні.")

    patches = []
    cols = st.columns(2)
    for i in range(int(patches_count)):
        with cols[i % 2]:
            st.markdown(f"**Латка #{i+1}**")
            w_mm = st.number_input(f"Ширина W, мм (латка #{i+1})", value=1000.0, min_value=10.0, step=10.0, format="%.0f", key=f"w_{i}")
            h_mm = st.number_input(f"Висота H, мм (латка #{i+1})", value=500.0, min_value=10.0, step=10.0, format="%.0f", key=f"h_{i}")
            patches.append((w_mm, h_mm))

    if st.button("Розрахувати"):
        overlap = overlap_cm / 100.0

        strips_dict = {}
        patch_rows = []
        strips_rows = []

        total_area = 0.0
        total_perimeter_geom = 0.0
        total_join_geom = 0.0  # сумарна довжина стиків смуг (геометрія)

        for idx, (w_mm, h_mm) in enumerate(patches, start=1):
            W = (w_mm + 2.0 * allowance_mm) / 1000.0
            H = (h_mm + 2.0 * allowance_mm) / 1000.0

            area_i = W * H
            total_area += area_i

            per_i = 2.0 * (W + H)
            total_perimeter_geom += per_i

            # ✅ визначаємо “розмір по різу” та “довжину смуги”
            if cut_dir == "Висоті H":
                cut_len = H          # по чому “крок 0,5”
                strip_len = W        # довжина смуги
                strip_label_len = W
            else:  # "Ширині W"
                cut_len = W
                strip_len = H
                strip_label_len = H

            n_strips = max(1, math.ceil(cut_len / STRIP_W_M))

            # ✅ стики між смугами: (n-1) * довжина смуги
            join_i = (n_strips - 1) * strip_len
            total_join_geom += join_i

            # ✅ запис смуг у відомість: (довжина смуги) × (фактична ширина по кроку)
            for s in range(n_strips):
                if s < n_strips - 1:
                    seg_w = STRIP_W_M
                else:
                    seg_w = cut_len - STRIP_W_M * (n_strips - 1)
                    if seg_w < 1e-9:
                        seg_w = STRIP_W_M

                # смуга: "довжина × ширина(по 0,5)"
                key = f"{strip_label_len:.2f}м × {seg_w:.2f}м"
                strips_dict[key] = strips_dict.get(key, 0) + 1

            patch_rows.append([
                f"Латка #{idx}",
                f"W={W:.2f} м; H={H:.2f} м",
                fmt_m2(area_i),
                fmt_m(per_i),
                f"{n_strips} смуг (0,50 м)",
                fmt_m(join_i)
            ])

        weld_geom_perimeter = total_perimeter_geom
        weld_geom_joins = total_join_geom
        weld_geom_total = weld_geom_perimeter + weld_geom_joins

        dep_perimeter = weld_geom_perimeter * seam_lines_perimeter * passes
        dep_joins = weld_geom_joins * seam_lines_on_lap * passes
        weld_deposition_total = dep_perimeter + dep_joins

        # ✅ додаткова площа перекриття від нахлесту між смугами:
        overlap_area = total_join_geom * overlap
        used_material_area = total_area + overlap_area

        sheet_area = 6.0 * 1.5
        total_sheets = math.ceil(used_material_area / sheet_area)

        electrodes_mass = weld_deposition_total * spec_consumption
        packs_needed = math.ceil(electrodes_mass / pack_mass)

        strips_items = sorted(strips_dict.items(), key=lambda x: x[0])
        for k, v in strips_items:
            strips_rows.append([k.replace("×", "x"), f"{v} шт"])

        # далі твій WEB/PDF вивід можна лишати як є (тільки поля "Висота смуги" прибрати або замінити на 0,5 м)


        # =============================
        # WEB ВИВІД
        # =============================
        st.markdown("## Результати розрахунків (латочний ремонт)")

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("### Підсумок")
            st.text(f"Тип: {tank_type}")
            st.text(f"ID: {tank_id}")
            st.text(f"Зав. №: {serial_no}")
            st.text(f"Дата: {repair_date.strftime('%d.%m.%Y')}")
            st.text("")
            st.text(f"Площа латок (геометрія): {total_area:.3f} м²")
            st.text(f"Додаткова площа перекриття (нахлест {overlap_cm:.1f} см): {overlap_area:.3f} м²")
            st.markdown(f"<p style='color:#b00020; font-weight:700;'>Фактична площа матеріалу: {used_material_area:.3f} м²</p>", unsafe_allow_html=True)
            st.text(f"Оцінка листів 6×1.5 м: {total_sheets} шт")
            st.text("")

            if note_weld:
                st.text(f"Шви (геометрія): {weld_geom_total:.2f} м")
                st.markdown(f"<p style='color:#006400; font-weight:700;'>Наплавлення: {weld_deposition_total:.2f} м</p>", unsafe_allow_html=True)
                st.text(f"  • периметр: {weld_geom_perimeter:.2f} м; стики смуг: {weld_geom_joins:.2f} м")
                st.text(f"  • ліній по периметру: {seam_lines_perimeter}; ліній на нахлесті: {seam_lines_on_lap}; проходів: {passes}")
            else:
                st.text(f"Загальна довжина швів: {weld_geom_total:.2f} м")

            st.markdown(f"<p style='color:#0033aa; font-weight:800;'>Електроди Ø{electrode_diam_mm:.1f} мм: {electrodes_mass:.2f} кг (~ {packs_needed} пач. по {pack_mass:.2f} кг)</p>", unsafe_allow_html=True)

        with c2:
            st.markdown("### Відомість латок")
            st.dataframe(
                {
                    "Латка": [r[0] for r in patch_rows],
                    "Розмір": [r[1] for r in patch_rows],
                    "Площа": [r[2] for r in patch_rows],
                    "Периметр": [r[3] for r in patch_rows],
                    "Смуги": [r[4] for r in patch_rows],
                    "Стики смуг": [r[5] for r in patch_rows],
                },
                use_container_width=True
            )

            st.markdown("### Відомість смуг (матеріал)")
            st.dataframe(
                {"Розмір смуги": [r[0] for r in strips_rows], "К-сть": [r[1] for r in strips_rows]},
                use_container_width=True
            )

        # =============================
        # PDF (2 сторінки, landscape)
        # =============================
        pdf_buffer = BytesIO()
        with PdfPages(pdf_buffer) as pdf:
            # Page 1: ремонт + шви
            fig1 = plt.figure(figsize=(11.69, 8.27))
            fig1.suptitle("ЗВІТ РОЗРАХУНКУ РЕМОНТУ (ЛАТОЧНИЙ)", fontsize=16, fontweight="bold", y=0.965)

            gs1 = gridspec.GridSpec(1, 2, figure=fig1, wspace=0.12, left=0.04, right=0.98, top=0.90, bottom=0.08)
            ax11 = fig1.add_subplot(gs1[0, 0])
            ax12 = fig1.add_subplot(gs1[0, 1])

            repair_rows = [
                ["Тип резервуара", tank_type],
                ["ID резервуара", tank_id],
                ["Заводський №", serial_no],
                ["Дата ремонту", repair_date.strftime("%d.%m.%Y")],
                ["Висота смуги", fmt_m(strip_h_m)],
                ["Нахлест між смугами", f"{overlap_cm:.1f} см"],
                ["Ліній шва по периметру", f"{seam_lines_perimeter}"],
                ["Ліній шва на нахлесті", f"{seam_lines_on_lap}"],
                ["Проходи", f"{passes}"],
                ["Площа латок (геометрія)", fmt_m2(total_area)],
                ["Дод. площа перекриття", fmt_m2(overlap_area)],
                ["Фактична площа матеріалу", fmt_m2(used_material_area)],
                ["Листи 6×1.5 м (оцінка)", f"{total_sheets} шт"],
            ]

            weld_rows = [
                ["Периметр латок (геометрія)", fmt_m(weld_geom_perimeter)],
                ["Стики смуг (геометрія)", fmt_m(weld_geom_joins)],
                ["Шви (геометрія) — разом", fmt_m(weld_geom_total)],
                ["Наплавлення по периметру", fmt_m(dep_perimeter)],
                ["Наплавлення на нахлесті", fmt_m(dep_joins)],
                ["Наплавлення — разом", fmt_m(weld_deposition_total)],
                ["Електроди Ø", f"{electrode_diam_mm:.1f} мм"],
                ["Питома витрата", f"{spec_consumption:.2f} кг/м"],
                ["Потрібна маса електродів", fmt_kg(electrodes_mass)],
                ["Маса 1 пачки", fmt_kg(pack_mass)],
                ["Кількість пачок", f"{packs_needed} шт"],
            ]

            make_table_ax(ax11, "1) Дані по ремонту", repair_rows, font_size=10, title_pad=16)
            make_table_ax(ax12, "2) Шви / наплавлення / електроди", weld_rows, font_size=10, title_pad=16)

            pdf.savefig(fig1)
            plt.close(fig1)

            # Page 2: латки + смуги
            fig2 = plt.figure(figsize=(11.69, 8.27))
            fig2.suptitle("Відомості (латки та смуги матеріалу)", fontsize=16, fontweight="bold", y=0.965)

            gs2 = gridspec.GridSpec(2, 2, figure=fig2,
                                    height_ratios=[1.0, 1.1],
                                    wspace=0.12, hspace=0.30,
                                    left=0.04, right=0.98, top=0.90, bottom=0.06)

            ax21 = fig2.add_subplot(gs2[0, :])
            ax22 = fig2.add_subplot(gs2[1, :])

            patch_table_rows = [[r[0], r[1], r[2], r[3], r[4], r[5]] for r in patch_rows]
            make_table_ax(
                ax21,
                "3) Відомість латок",
                patch_table_rows,
                col_labels=("Латка", "Розмір", "Площа", "Периметр", "Смуги", "Стики смуг"),
                font_size=9,
                title_pad=14,
                scale_y=1.15
            )

            make_table_ax(
                ax22,
                "4) Відомість смуг (розкрій/кількість)",
                strips_rows,
                col_labels=("Розмір смуги", "К-сть"),
                font_size=9,
                title_pad=14,
                scale_y=1.20
            )

            pdf.savefig(fig2)
            plt.close(fig2)

        st.markdown("---")
        st.markdown("### Збереження результатів")
        st.download_button(
            label="Зберегти звіт як PDF",
            data=pdf_buffer.getvalue(),
            file_name="резервуар_звіт_латки.pdf",
            mime="application/pdf"
        )
