import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
from datetime import date
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg

VERSION = "1.7.1"

# -----------------------------
# Streamlit
# -----------------------------
st.set_page_config(layout="wide")
st.title(f"Розрахунок ремонту резервуара та зварних швів — v{VERSION}")

# -----------------------------
# Helpers
# -----------------------------
def clamp(x, a, b):
    return max(a, min(b, x))

def primitive_F(R: float, y: float) -> float:
    y = clamp(y, -R, R)
    inside = max(0.0, R * R - y * y)
    return y * math.sqrt(inside) + R * R * math.asin(clamp(y / R, -1.0, 1.0))

def chord_len(R: float, y: float) -> float:
    if y <= -R or y >= R:
        return 0.0
    return 2.0 * math.sqrt(max(0.0, R*R - y*y))

def build_patterns(Wrem: float):
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

def make_table_ax(ax, title, rows, col_labels=("Параметр", "Значення"),
                  font_size=10, title_pad=14, scale_y=1.30, bbox=(0.0, 0.02, 1.0, 0.86)):
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

# -----------------------------
# Реквізити звіту (додано ID резервуара)
# -----------------------------
meta_c1, meta_c2, meta_c3, meta_c4 = st.columns([2, 1, 1, 1])
with meta_c1:
    tank_type = st.text_input("Тип резервуара", value="Паливний резервуар")
with meta_c2:
    tank_id = st.text_input("ID резервуара", value="—")  # ✅ НОВЕ
with meta_c3:
    serial_no = st.text_input("Заводський №", value="—")
with meta_c4:
    repair_date = st.date_input("Дата ремонту", value=date.today())

st.markdown("---")

# -----------------------------
# Режим ремонту
# -----------------------------
mode = st.radio(
    "Режим роботи",
    ["Ремонт резервуара (повний)", "Ремонт резервуара (латочний)"],
    horizontal=True
)
is_full = (mode == "Ремонт резервуара (повний)")
is_patch = not is_full

# -----------------------------
# Електроди / витрата (спільні)
# -----------------------------
st.markdown("### Електроди / витрата")
t1, t2, t3 = st.columns(3)
with t1:
    electrode_diam_mm = st.number_input("Діаметр електрода, мм", value=3.2, min_value=1.6, step=0.2, format="%.1f")
with t2:
    spec_consumption = st.number_input("Питома витрата електродів, кг/м наплавлення", value=0.25, min_value=0.01, step=0.01, format="%.2f")
with t3:
    pack_mass = st.number_input("Маса 1 пачки електродів, кг", value=2.5, min_value=0.5, step=0.5, format="%.1f")

st.markdown("---")

# -----------------------------
# Вхідні дані (повний ремонт)
# -----------------------------
st.markdown("### Вхідні дані")
col1, col2, col3, col4 = st.columns(4)
with col1:
    D = st.number_input("Діаметр резервуара, м", value=2.5, min_value=0.1, step=0.1, format="%.2f", disabled=is_patch)
with col2:
    L = st.number_input("Довжина резервуара, м", value=4.3, min_value=0.1, step=0.1, format="%.2f", disabled=is_patch)
with col3:
    Wrem = st.number_input("Ширина ремонтованої ділянки (по розгортці), м", value=1.0, min_value=0.1, step=0.1, format="%.2f", disabled=is_patch)
with col4:
    overlap_cm = st.number_input("Нахлест, см", value=5.0, min_value=0.0, step=0.5, format="%.1f", disabled=is_patch)

# -----------------------------
# Налаштування швів (повний ремонт)
# -----------------------------
st.markdown("### Налаштування під технологію (шви)")
w1, w2, w3, w4 = st.columns(4)
with w1:
    overlap_weld_lines = st.selectbox("Ліній шва на нахлесті (смуга-смуга)", [1, 2], index=1, disabled=is_patch)
with w2:
    edge_weld_lines = st.selectbox("Ліній шва по контуру ремонту (краї)", [1, 2], index=0, disabled=is_patch)
with w3:
    include_shell_to_bottom = st.checkbox("Враховувати шви приєднання циліндра до днищ", value=True, disabled=is_patch)
with w4:
    shell_to_bottom_lines = st.selectbox("Ліній шва на приєднанні (коло стику)", [1, 2], index=1, disabled=(is_patch or (not include_shell_to_bottom)))

# -----------------------------
# Латочний ремонт
# -----------------------------
st.markdown("---")
st.markdown("### Латочний ремонт (прямокутні латки)")

p1, p2, p3, p4 = st.columns([1.2, 1.2, 1.2, 1.4])
with p1:
    n_patches = st.number_input("Кількість латок, шт", value=2, min_value=1, step=1, disabled=is_full)
with p2:
    patch_weld_lines = st.selectbox("Ліній шва на 1 латці (по периметру)", [1, 2], index=1, disabled=is_full)
with p3:
    patch_allowance_pct = st.number_input("Запас матеріалу на підгін/обрізки, %", value=5.0, min_value=0.0, step=1.0, format="%.1f", disabled=is_full)
with p4:
    patch_note = st.text_input("Примітка (за бажанням)", value="Прямокутні латки по місцях дефектів", disabled=is_full)

patch_dims = []
if is_patch:
    st.markdown("#### Розміри латок (мм)")
    for i in range(int(n_patches)):
        cA, cB = st.columns(2)
        with cA:
            w_mm = st.number_input(f"Латка {i+1}: ширина, мм", value=300.0, min_value=10.0, step=10.0, format="%.0f", key=f"pw_{i}")
        with cB:
            h_mm = st.number_input(f"Латка {i+1}: висота, мм", value=200.0, min_value=10.0, step=10.0, format="%.0f", key=f"ph_{i}")
        patch_dims.append((w_mm, h_mm))

# -----------------------------
# Розрахунок
# -----------------------------
if st.button("Розрахувати"):
    # ============ ЛАТОЧНИЙ ============
    if is_patch:
        patches_m = []
        for (w_mm, h_mm) in patch_dims:
            patches_m.append((w_mm / 1000.0, h_mm / 1000.0))

        area_patches = sum(w*h for w, h in patches_m)
        area_used = area_patches * (1.0 + patch_allowance_pct/100.0)

        weld_per_patch = [2.0*(w+h) * patch_weld_lines for (w, h) in patches_m]
        weld_total = sum(weld_per_patch)

        electrodes_mass = weld_total * spec_consumption
        packs_needed = math.ceil(electrodes_mass / pack_mass)

        st.markdown("## Результати (латочний ремонт)")
        left_lines = [
            f"Тип резервуара: {tank_type}",
            f"ID резервуара: {tank_id}",  # ✅
            f"Заводський №: {serial_no}",
            f"Дата ремонту: {repair_date.strftime('%d.%m.%Y')}",
            f"Примітка: {patch_note}",
            "",
            f"Кількість латок: {int(n_patches)} шт",
            f"Площа латок (чиста): {fmt_m2(area_patches)}",
            f"Запас матеріалу: {patch_allowance_pct:.1f} %",
            f"Фактична площа матеріалу: {fmt_m2(area_used)}",
            "",
            f"Ліній шва на латці: {patch_weld_lines}",
            f"Загальна довжина наплавлення: {fmt_m(weld_total)}",
            "",
            f"Електроди Ø{electrode_diam_mm:.1f} мм: {fmt_kg(electrodes_mass)}  (~ {packs_needed} пач. по {fmt_kg(pack_mass)})",
        ]

        patch_rows = []
        for i, (w, h) in enumerate(patches_m, start=1):
            patch_rows.append([f"Латка {i}", f"{w*1000:.0f}×{h*1000:.0f} мм", fmt_m2(w*h), fmt_m(weld_per_patch[i-1])])

        c1, c2 = st.columns([1.2, 1.0])
        with c1:
            st.markdown("### Підсумки")
            for line in left_lines:
                if line.startswith("Фактична площа"):
                    st.markdown(f"<p style='color:#b00020; font-weight:700;'>{line}</p>", unsafe_allow_html=True)
                elif line.startswith("Загальна довжина"):
                    st.markdown(f"<p style='color:#006400; font-weight:800;'>{line}</p>", unsafe_allow_html=True)
                elif line.startswith("Електроди"):
                    st.markdown(f"<p style='color:#0033aa; font-weight:800;'>{line}</p>", unsafe_allow_html=True)
                else:
                    st.text(line)

        with c2:
            st.markdown("### Відомість латок")
            st.table([{"Латка": r[0], "Розмір": r[1], "Площа": r[2], "Наплавлення": r[3]} for r in patch_rows])

        pdf_buffer = BytesIO()
        with PdfPages(pdf_buffer) as pdf:
            fig1 = plt.figure(figsize=(11.69, 8.27))
            fig1.suptitle("ЗВІТ: ЛАТОЧНИЙ РЕМОНТ (ПРЯМОКУТНІ ЛАТКИ)", fontsize=16, fontweight="bold", y=0.965)

            gs = gridspec.GridSpec(2, 2, figure=fig1,
                                   height_ratios=[0.30, 0.70],
                                   wspace=0.12, hspace=0.25,
                                   left=0.04, right=0.98, top=0.90, bottom=0.06)

            ax_meta = fig1.add_subplot(gs[0, :])
            ax_meta.axis("off")
            meta_lines = [
                f"Тип резервуара: {tank_type}    |    ID: {tank_id}    |    Заводський №: {serial_no}    |    Дата: {repair_date.strftime('%d.%m.%Y')}",
                f"Кількість латок: {int(n_patches)} шт   |   Запас матеріалу: {patch_allowance_pct:.1f} %",
                f"Наплавлення: {weld_total:.2f} м   |   Електроди Ø{electrode_diam_mm:.1f} мм: {electrodes_mass:.2f} кг (~ {packs_needed} пач.)",
            ]
            ax_meta.text(0.01, 0.92, meta_lines[0], fontsize=11, fontweight="bold", va="top")
            ax_meta.text(0.01, 0.62, meta_lines[1], fontsize=10, va="top")
            ax_meta.text(0.01, 0.35, meta_lines[2], fontsize=10, va="top")
            if patch_note.strip():
                ax_meta.text(0.01, 0.10, f"Примітка: {patch_note}", fontsize=10, va="top")

            ax_left = fig1.add_subplot(gs[1, 0])
            ax_right = fig1.add_subplot(gs[1, 1])

            repair_rows = [
                ["Тип резервуара", tank_type],
                ["ID резервуара", tank_id],
                ["Заводський №", serial_no],
                ["Дата ремонту", repair_date.strftime("%d.%m.%Y")],
                ["Кількість латок", f"{int(n_patches)} шт"],
                ["Площа латок (чиста)", fmt_m2(area_patches)],
                ["Запас матеріалу", f"{patch_allowance_pct:.1f} %"],
                ["Фактична площа матеріалу", fmt_m2(area_used)],
            ]
            weld_rows = [
                ["Ліній шва на латці", f"{patch_weld_lines}"],
                ["Загальна довжина наплавлення", fmt_m(weld_total)],
                ["Питома витрата", f"{spec_consumption:.2f} кг/м"],
                ["Маса електродів", fmt_kg(electrodes_mass)],
                ["Маса пачки", fmt_kg(pack_mass)],
                ["Кількість пачок", f"{packs_needed} шт"],
            ]
            make_table_ax(ax_left, "1) Дані по ремонту", repair_rows, font_size=10, title_pad=16)
            make_table_ax(ax_right, "2) Шви та електроди", weld_rows, font_size=10, title_pad=16)
            pdf.savefig(fig1)
            plt.close(fig1)

            fig2 = plt.figure(figsize=(11.69, 8.27))
            fig2.suptitle("Відомість латок (деталізація)", fontsize=16, fontweight="bold", y=0.965)
            ax = fig2.add_subplot(111)
            rows = [[r[0], r[1], r[2], r[3]] for r in patch_rows]
            make_table_ax(ax, "3) Перелік латок", rows,
                          col_labels=("Латка", "Розмір", "Площа", "Наплавлення"),
                          font_size=10, title_pad=16, scale_y=1.20, bbox=(0.0, 0.02, 1.0, 0.90))
            pdf.savefig(fig2)
            plt.close(fig2)

        st.markdown("---")
        st.markdown("### Збереження результатів")
        st.download_button("Зберегти звіт як PDF", pdf_buffer.getvalue(), "резервуар_звіт_латки.pdf", "application/pdf")


    # ==========================================================
    # B) ПОВНИЙ РЕМОНТ (як ми робили)
    # ==========================================================
    else:
        # валідація
        if D <= 0 or L <= 0 or Wrem <= 0 or overlap_cm < 0:
            st.error("Введіть додатні D, L, Wrem та невід'ємний нахлест.")
            st.stop()

        overlap = overlap_cm / 100.0
        R = D / 2.0
        h_smuha = 0.5
        circumference = 2.0 * math.pi * R

        if Wrem > circumference + 1e-9:
            st.warning(f"Wrem = {Wrem:.2f} м більша за довжину кола {circumference:.2f} м. Розгортка буде умовною.")

        # ====== Днище: висота ремонту ======
        alpha = (Wrem / 2.0) / R
        Hcrit = R * (1.0 - math.cos(alpha))
        n_bot = max(1, math.ceil(Hcrit / h_smuha))
        y_cut = -R + Hcrit

        # ============================================
        # 1) ДНИЩЕ (графік + площі)
        # ============================================
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

            width = 0.0 if abs(y_ref) >= R else 2.0 * math.sqrt(max(0.0, R*R - y_ref*y_ref))
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

            ov_h = h_smuha + (overlap/2.0 if j < n_bot - 1 else 0.0)
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

            ax_bot.text(0, y_bot + h_smuha/2.0, f"S{j+1}\n{area_strip:.3f} м²",
                        ha="center", va="center", fontsize=8, zorder=7)

        cum_area_bot = sum(areas_bot)
        total_area_both_bottoms = 2.0 * cum_area_bot

        m = R * 0.1
        ax_bot.set_xlim(-R - m, R + m)
        ax_bot.set_ylim(-R - m, R + m)
        ax_bot.set_xlabel("x (м)")
        ax_bot.set_ylabel("y (м)")
        ax_bot.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

        # ============================================
        # 2) ЦИЛІНДР (графік 1:1)
        # ============================================
        full_rows = 1
        covered_height = h_smuha
        while covered_height + 1e-9 < L:
            full_rows += 1
            covered_height += (h_smuha - overlap)

        patterns = build_patterns(Wrem)

        fig_cyl, ax_cyl = plt.subplots(figsize=(6.2, 6.2))
        ax_cyl.set_aspect("equal", adjustable="box")
        ax_cyl.set_title("Розгорнута поверхня (масштаб 1:1, зона ремонту)", fontsize=12, fontweight="bold")

        ax_cyl.add_patch(plt.Rectangle((-circumference/2.0, 0), circumference, L,
                                       edgecolor="black", facecolor="#d0d0d0", linewidth=1, zorder=1))
        ax_cyl.axvline(0, color="red", linestyle="--", linewidth=1, alpha=0.7)
        ax_cyl.axhline(L, color="red", linestyle="--", linewidth=1, alpha=0.7)

        x_start = -Wrem / 2.0
        area_cyl_with_overlap = 0.0
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
                        ov_x = x_off - overlap/2.0
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

                ax_cyl.text(x_off + seg/2.0, y_off + base_h/2.0, f"{seg:.2f} м",
                            ha="center", va="center", fontsize=9, zorder=5)

                key_cyl = f"{round(seg, 2):>5.2f}м × {h_smuha:>3.2f}м"
                smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1

                x_off += seg

        ax_cyl.set_xlim(-Wrem/2.0 - 0.6, Wrem/2.0 + 0.6)
        ax_cyl.set_ylim(0, max_top_for_ylim)
        ax_cyl.set_xlabel("довжина поверхні (м)")
        ax_cyl.set_ylabel("висота (м)")
        ax_cyl.grid(True, linestyle="--", linewidth=0.5, alpha=0.25)

        # ============================================
        # 3) ПЛОЩІ
        # ============================================
        area_cyl = Wrem * L
        real_repair_area = total_area_both_bottoms + area_cyl

        area_bottoms_with_overlap_fact = 2.0 * sum(
            w * (h + (overlap/2.0 if idx < n_bot - 1 else 0.0))
            for idx, (w, h) in enumerate(zip(widths_bot, used_heights_bot))
        )
        used_material_area = area_bottoms_with_overlap_fact + area_cyl_with_overlap

        sheet_area = 6.0 * 1.5
        total_sheets = math.ceil(used_material_area / sheet_area)

        # ============================================
        # 4) ШВИ (прозора логіка)
        # ============================================
        # Горизонтальні: нижній край + верхній край + міжрядові стики (з лініями на нахлесті)
        weld_cyl_h = (2.0 * Wrem * edge_weld_lines) + ((full_rows - 1) * Wrem * overlap_weld_lines)

        # Вертикальні: лівий+правий край ремонту (контур)
        weld_cyl_v_edges = 2.0 * L * edge_weld_lines

        # Внутрішні вертикальні (між сегментами в рядах): вважаємо як стики смуг по висоті рядка.
        # Для технології беремо ліній = overlap_weld_lines (можеш змінити при потребі пізніше).
        weld_cyl_v_internal = 0.0
        for i in range(full_rows):
            y = i * (h_smuha - overlap)
            if y >= L:
                break
            h_eff = min(h_smuha, L - y)
            pat = patterns[i % len(patterns)]
            internal_joints = max(0, len(pat) - 1)
            weld_cyl_v_internal += internal_joints * h_eff * overlap_weld_lines

        weld_cyl_v = weld_cyl_v_edges + weld_cyl_v_internal
        weld_cyl_total = weld_cyl_h + weld_cyl_v

        # днище (геометричні шви смуг) — як раніше (по хордах + боковий контур)
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

        weld_bot_h_one = sum(chord_len(R, y) for y in levels) * overlap_weld_lines
        theta_top = math.asin(clamp(y_top / R, -1.0, 1.0))
        arc_len_one_side = R * (theta_top - (-math.pi / 2.0))
        weld_bot_side_one = (2.0 * arc_len_one_side) * edge_weld_lines

        weld_bot_one = weld_bot_h_one + weld_bot_side_one
        weld_bottoms_total = 2.0 * weld_bot_one

        weld_shell_to_bottom = 0.0
        if include_shell_to_bottom:
            weld_shell_to_bottom = 2.0 * circumference * shell_to_bottom_lines

        weld_total = weld_cyl_total + weld_bottoms_total + weld_shell_to_bottom

        # ============================================
        # 5) Електроди
        # ============================================
        electrodes_mass = weld_total * spec_consumption
        packs_needed = math.ceil(electrodes_mass / pack_mass)

        # ============================================
        # ВЕБ-ВИВІД
        # ============================================
        st.markdown("## Результати (повний ремонт)")

        left_lines = [
            f"Тип резервуара: {tank_type}",
            f"Заводський №: {serial_no}",
            f"Дата ремонту: {repair_date.strftime('%d.%m.%Y')}",
            f"Діаметр D: {D:.2f} м  |  Радіус R: {R:.2f} м",
            f"Довжина L: {L:.2f} м",
            f"Ширина ремонту Wrem: {Wrem:.2f} м",
            f"Нахлест: {overlap_cm:.1f} см",
            f"Висота ремонту днища H: {Hcrit:.3f} м",
            f"Кількість рядів на циліндрі: {full_rows} шт",
            "",
            f"Площа 1 днища: {cum_area_bot:.3f} м²",
            f"Площа 2 днищ: {total_area_both_bottoms:.3f} м²",
            f"Площа циліндричної частини: {area_cyl:.3f} м²",
            f"Реально ремонтована площа: {real_repair_area:.3f} м²",
            f"Фактична площа матеріалу (з нахлестом): {used_material_area:.3f} м²",
            f"Оцінка листів 6×1.5 м: {total_sheets} шт",
            "",
            f"Довжина швів (циліндр): {weld_cyl_total:.2f} м",
            f"Довжина швів (днища): {weld_bottoms_total:.2f} м",
            f"Приєднання днищ до циліндра: {weld_shell_to_bottom:.2f} м",
            f"Загальна довжина наплавлення: {weld_total:.2f} м",
            "",
            f"Електроди Ø{electrode_diam_mm:.1f} мм: {electrodes_mass:.2f} кг  (~ {packs_needed} пач. по {pack_mass:.2f} кг)",
        ]

        strips_items = sorted(smuhaDict.items(), key=lambda x: x[0])
        right_lines = ["Розмір смуги           К-сть"]
        for k, v in strips_items:
            right_lines.append(f"{k:>16s}   {v:>3d} шт")

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("### Графіки")
            st.pyplot(fig_bot, clear_figure=False)
            st.pyplot(fig_cyl, clear_figure=False)

        with c2:
            st.markdown("### Підсумкові дані")
            for line in left_lines:
                if line.startswith("Фактична площа"):
                    st.markdown(f"<p style='color:#b00020; font-weight:700;'>{line}</p>", unsafe_allow_html=True)
                elif line.startswith("Загальна довжина"):
                    st.markdown(f"<p style='color:#006400; font-weight:800;'>{line}</p>", unsafe_allow_html=True)
                elif line.startswith("Електроди"):
                    st.markdown(f"<p style='color:#0033aa; font-weight:800;'>{line}</p>", unsafe_allow_html=True)
                else:
                    st.text(line)

            st.markdown("#### Відомість смуг")
            for line in right_lines:
                st.text(line)

        # ============================================
        # PDF (3 сторінки)
        # ============================================
        pdf_buffer = BytesIO()
        with PdfPages(pdf_buffer) as pdf:
            # Page 1 (landscape): шапка + 2 графіки
            fig1 = plt.figure(figsize=(11.69, 8.27))
            fig1.suptitle("ЗВІТ РОЗРАХУНКУ РЕМОНТНИХ СМУГ", fontsize=16, fontweight="bold", y=0.965)

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
                f"Тип резервуара: {tank_type}    |    Заводський №: {serial_no}    |    Дата ремонту: {repair_date.strftime('%d.%m.%Y')}",
                f"D = {D:.2f} м;  L = {L:.2f} м;  Wrem = {Wrem:.2f} м;  Нахлест = {overlap_cm:.1f} см",
                f"Загальна довжина наплавлення: {weld_total:.2f} м    |    Матеріал (факт): {used_material_area:.3f} м²    |    Листи 6×1.5 м: {total_sheets} шт",
                f"Електроди Ø{electrode_diam_mm:.1f} мм: {electrodes_mass:.2f} кг  (~ {packs_needed} пач. по {pack_mass:.2f} кг; питома витрата {spec_consumption:.2f} кг/м)",
            ]
            ax_meta.text(0.01, 0.92, meta_lines[0], fontsize=11, fontweight="bold", va="top")
            ax_meta.text(0.01, 0.62, meta_lines[1], fontsize=10, va="top")
            ax_meta.text(0.01, 0.35, meta_lines[2], fontsize=10, va="top")
            ax_meta.text(0.01, 0.10, meta_lines[3], fontsize=10, va="top")

            img_bot = fig_to_img_array(fig_bot, dpi=160)
            img_cyl = fig_to_img_array(fig_cyl, dpi=160)
            ax_a.imshow(img_bot); ax_b.imshow(img_cyl)
            ax_a.set_title("Рис. 1 — Схема смуг на днищі", fontsize=11, fontweight="bold", pad=8)
            ax_b.set_title("Рис. 2 — Розгортка циліндра (масштаб 1:1)", fontsize=11, fontweight="bold", pad=8)
            ax_a.axis("off"); ax_b.axis("off")

            pdf.savefig(fig1); plt.close(fig1)

            # Page 2: ремонт + шви
            fig2 = plt.figure(figsize=(11.69, 8.27))
            fig2.suptitle("Розрахунок ремонту та зварних швів", fontsize=16, fontweight="bold", y=0.965)
            gs2 = gridspec.GridSpec(1, 2, figure=fig2, wspace=0.12, left=0.04, right=0.98, top=0.90, bottom=0.08)
            ax21 = fig2.add_subplot(gs2[0, 0]); ax22 = fig2.add_subplot(gs2[0, 1])

            repair_rows = [
                ["Тип резервуара", tank_type],
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
                ["Загальна довжина наплавлення", fmt_m(weld_total)],
            ]
            make_table_ax(ax21, "1) Дані по ремонту", repair_rows, font_size=10, title_pad=16)
            make_table_ax(ax22, "2) Довжина зварних швів (деталізація)", weld_rows, font_size=10, title_pad=16)
            pdf.savefig(fig2); plt.close(fig2)

            # Page 3: матеріали + електроди + смуги
            fig3 = plt.figure(figsize=(11.69, 8.27))
            fig3.suptitle("Матеріали та електроди", fontsize=16, fontweight="bold", y=0.965)

            gs3 = gridspec.GridSpec(2, 2, figure=fig3,
                                    height_ratios=[1.0, 1.05],
                                    wspace=0.12, hspace=0.28,
                                    left=0.04, right=0.98, top=0.90, bottom=0.06)
            ax31 = fig3.add_subplot(gs3[0, 0])
            ax32 = fig3.add_subplot(gs3[0, 1])
            ax33 = fig3.add_subplot(gs3[1, :])

            areas_rows = [
                ["Площа 1 днища", fmt_m2(cum_area_bot)],
                ["Площа 2 днищ", fmt_m2(total_area_both_bottoms)],
                ["Площа циліндричної частини (ремонт)", fmt_m2(area_cyl)],
                ["Реально ремонтована площа", fmt_m2(real_repair_area)],
                ["Фактична площа матеріалу (з нахлестом)", fmt_m2(used_material_area)],
                ["Листи 6×1.5 м (оцінка)", f"{total_sheets} шт"],
            ]
            electrodes_rows = [
                ["Діаметр електрода", f"{electrode_diam_mm:.1f} мм"],
                ["Питома витрата", f"{spec_consumption:.2f} кг/м"],
                ["Загальна довжина наплавлення", fmt_m(weld_total)],
                ["Потрібна маса електродів", fmt_kg(electrodes_mass)],
                ["Маса 1 пачки", fmt_kg(pack_mass)],
                ["Кількість пачок", f"{packs_needed} шт"],
            ]
            strips_rows = [[k.replace("×", "x"), f"{v} шт"] for k, v in strips_items]
            make_table_ax(ax31, "3) Площі та матеріали", areas_rows, font_size=10, title_pad=16)
            make_table_ax(ax32, "4) Електроди (оцінка)", electrodes_rows, font_size=10, title_pad=16)
            make_table_ax(ax33, "5) Відомість смуг (розкрій/кількість)", strips_rows,
                          col_labels=("Розмір смуги", "К-сть"), font_size=9, title_pad=14, scale_y=1.22)
            pdf.savefig(fig3); plt.close(fig3)

        # -----------------------------
        # Завантаження
        # -----------------------------
        st.markdown("---")
        st.markdown("### Збереження результатів")
        st.download_button("Зберегти звіт як PDF", pdf_buffer.getvalue(), "резервуар_звіт.pdf", "application/pdf")

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
