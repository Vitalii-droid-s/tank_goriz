import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
from datetime import date
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages

VERSION = "1.5"

# -----------------------------
# Streamlit
# -----------------------------
st.set_page_config(layout="wide")
st.title(f"Розрахунок смуг для резервуара — v{VERSION}")

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
    """
    - Wrem 4..6: "шахматка" як було
    - інакше: ріжемо на сегменти 3.0 м + залишок
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

def fmt_m(x):  # meters
    return f"{x:.2f} м"

def fmt_m2(x):  # square meters
    return f"{x:.3f} м²"

def fmt_kg(x):
    return f"{x:.2f} кг"

def make_table_ax(ax, title, rows, col_labels=("Параметр", "Значення"),
                  font_size=10, title_pad=14, scale_y=1.35):
    """
    Проф. таблиця без налізання тексту на заголовок:
    - заголовок (ax.set_title) + відступ
    - таблиця через bbox нижче заголовка
    """
    ax.axis("off")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=title_pad)

    tbl = ax.table(
        cellText=rows,
        colLabels=list(col_labels),
        colLoc="center",
        cellLoc="left",
        loc="center",
        bbox=[0.0, 0.02, 1.0, 0.86]  # важливо: опускаємо таблицю нижче заголовка
    )

    tbl.auto_set_font_size(False)
    tbl.set_fontsize(font_size)
    tbl.scale(1.0, scale_y)

    # оформлення header
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.8)
        if r == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#f2f2f2")
            cell.set_height(cell.get_height() * 1.10)

    return tbl

# -----------------------------
# Верхній блок реквізитів (для звіту)
# -----------------------------
meta_c1, meta_c2, meta_c3 = st.columns([2, 1, 1])
with meta_c1:
    tank_type = st.text_input("Тип резервуара", value="Паливний резервуар РГС-20")
with meta_c2:
    serial_no = st.text_input("Заводський №", value="215")
with meta_c3:
    repair_date = st.date_input("Дата ремонту", value=date.today())

# -----------------------------
# Вхідні параметри геометрії/технології
# -----------------------------
col1, col2, col3, col4 = st.columns(4)
with col1:
    D = st.number_input("Діаметр резервуара, м", value=2.5, min_value=0.1, step=0.1, format="%.2f")
with col2:
    L = st.number_input("Довжина резервуара, м", value=4.3, min_value=0.1, step=0.1, format="%.2f")
with col3:
    Wrem = st.number_input("Ширина ремонтованої ділянки (по розгортці), м", value=1.0, min_value=0.1, step=0.1, format="%.2f")
with col4:
    overlap_cm = st.number_input("Нахлест, см", value=5.0, min_value=0.0, step=0.5, format="%.1f")

t1, t2, t3 = st.columns(3)
with t1:
    electrode_diam_mm = st.number_input("Діаметр електрода, мм", value=3.2, min_value=1.6, step=0.2, format="%.1f")
with t2:
    spec_consumption = st.number_input("Питома витрата електродів, кг/м шва", value=0.25, min_value=0.01, step=0.01, format="%.2f")
with t3:
    pack_mass = st.number_input("Маса 1 пачки електродів, кг", value=2.5, min_value=0.5, step=0.5, format="%.1f")

include_shell_to_bottom = st.checkbox("Враховувати шви приєднання циліндричної частини до днищ", value=True)

# -----------------------------
# Run
# -----------------------------
if st.button("Розрахувати"):
    if D <= 0 or L <= 0 or Wrem <= 0 or overlap_cm < 0:
        st.error("Введіть додатні D, L, Wrem та невід'ємний нахлест.")
        st.stop()

    overlap = overlap_cm / 100.0
    R = D / 2.0
    h_smuha = 0.5

    circumference = 2.0 * math.pi * R
    if Wrem > circumference + 1e-9:
        st.warning(f"Wrem = {Wrem:.2f} м більша за довжину кола {circumference:.2f} м. Розгортка буде умовною.")

    # =========================
    # Днище: висота ремонту
    # =========================
    alpha = (Wrem / 2.0) / R
    Hcrit = R * (1.0 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))
    y_cut = -R + Hcrit

    # =========================
    # 1) Днище (графік + площа)
    # =========================
    widths_bot, used_heights_bot, areas_bot = [], [], []
    smuhaDict = {}

    fig_bot, ax_bot = plt.subplots(figsize=(6, 6))
    ax_bot.set_aspect("equal")
    ax_bot.set_title("Днище резервуара (смуги та зона нахлесту)", fontsize=12, fontweight="bold")

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

        # width at y_ref (твоя логіка)
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

        # облік смуг (дві симетричні)
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

    # =========================
    # 2) Циліндр (графік)
    # =========================
    full_rows = 1
    covered_height = h_smuha
    while covered_height + 1e-9 < L:
        full_rows += 1
        covered_height += (h_smuha - overlap)

    patterns = build_patterns(Wrem)

    fig_cyl, ax_cyl = plt.subplots(figsize=(8, 4.5))
    ax_cyl.set_aspect("equal", adjustable="box")
    ax_cyl.set_title("Розгорнута поверхня (масштаб 1:1, зона ремонту)", fontsize=12, fontweight="bold")

    # корпус (для контексту)
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

    ax_cyl.set_xlim(-Wrem/2.0 - 0.6, Wrem/2.0 + 0.6)  # щоб клієнт бачив, що зона ремонту = Wrem
    ax_cyl.set_ylim(0, max_top_for_ylim)
    ax_cyl.set_xlabel("довжина поверхні (м)")
    ax_cyl.set_ylabel("висота (м)")
    ax_cyl.grid(True, linestyle="--", linewidth=0.5, alpha=0.25)

    # =========================
    # 3) Площі
    # =========================
    area_cyl = Wrem * L
    real_repair_area = total_area_both_bottoms + area_cyl

    area_bottoms_with_overlap_fact = 2.0 * sum(
        w * (h + (overlap/2.0 if idx < n_bot - 1 else 0.0))
        for idx, (w, h) in enumerate(zip(widths_bot, used_heights_bot))
    )
    used_material_area = area_bottoms_with_overlap_fact + area_cyl_with_overlap

    sheet_area = 6.0 * 1.5
    total_sheets = math.ceil(used_material_area / sheet_area)

    # =========================
    # 4) Довжина зварних швів (ПРОФ. логіка)
    #    Рахуємо "довжину наплавлення" як суму швів по периметру латок:
    #    - циліндр: горизонтальні (між рядами) + вертикальні (між сегментами) + крайові
    #    - днища: хорди по рівнях + бічні дуги
    #    - приєднання днищ до циліндра: довжина кола * 2 (якщо увімкнено)
    # =========================
    # Циліндр — горизонтальні шви:
    # (full_rows+1)*Wrem — верх+низ + шви між рядами (уздовж ширини ремонту)
    weld_cyl_h = (full_rows + 1) * Wrem

    # Циліндр — вертикальні:
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
        # вертикальні: (внутрішні стики) + 2 крайові
        n_vertical = (len(pat) - 1) + 2
        weld_cyl_v += n_vertical * h_eff

    weld_cyl_total = weld_cyl_h + weld_cyl_v

    # Днище (1 шт.) — горизонтальні по рівнях (хорди):
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

    # Днище (1 шт.) — бічні дуги (по межі ремонту зліва+справа):
    theta_top = math.asin(clamp(y_top / R, -1.0, 1.0))
    arc_len_one_side = R * (theta_top - (-math.pi / 2.0))
    weld_bot_side_one = 2.0 * arc_len_one_side

    weld_bot_one = weld_bot_h_one + weld_bot_side_one
    weld_bottoms_total = 2.0 * weld_bot_one

    # Приєднання днищ до циліндра (2 кола)
    weld_shell_to_bottom = 0.0
    if include_shell_to_bottom:
        weld_shell_to_bottom = 2.0 * circumference

    weld_total = weld_cyl_total + weld_bottoms_total + weld_shell_to_bottom

    # =========================
    # 5) Електроди
    # =========================
    electrodes_mass = weld_total * spec_consumption
    packs_needed = math.ceil(electrodes_mass / pack_mass)

    # =========================
    # 6) Смуги (таблиця зі стор.5 -> переносимо на стор.4)
    # =========================
    strips_items = sorted(smuhaDict.items(), key=lambda x: x[0])

    # -----------------------------
    # UI output
    # -----------------------------
    st.markdown("## Результати розрахунків")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("### Графіки")
        st.pyplot(fig_bot, clear_figure=False)
        st.pyplot(fig_cyl, clear_figure=False)

    with c2:
        st.markdown("### Ключові результати")
        st.write(f"**Загальна довжина зварних швів:** {weld_total:.2f} м")
        st.write(f"**Фактична площа матеріалу (з нахлестом):** {used_material_area:.3f} м²")
        st.write(f"**Оцінка листів 6×1.5 м:** {total_sheets} шт")
        st.write(f"**Електроди Ø{electrode_diam_mm:.1f} мм:** {electrodes_mass:.2f} кг  (~ {packs_needed} пач. по {pack_mass:.2f} кг)")

    # -----------------------------
    # PDF report (професійна верстка)
    # -----------------------------
    pdf_buffer = BytesIO()

    with PdfPages(pdf_buffer) as pdf:
        # --- PAGE 1: титул/вхідні/ключові ---
        fig1 = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
        fig1.clf()
        ax = fig1.add_subplot(111)
        ax.axis("off")

        ax.text(0.5, 0.95, "ЗВІТ РОЗРАХУНКУ РЕМОНТНИХ СМУГ", ha="center", va="top",
                fontsize=18, fontweight="bold")
        ax.text(0.02, 0.90, f"Версія розрахунку: v{VERSION}", fontsize=11)

        ax.text(0.02, 0.86, f"Тип резервуара: {tank_type}", fontsize=12)
        ax.text(0.02, 0.83, f"Заводський №: {serial_no}", fontsize=12)
        ax.text(0.02, 0.80, f"Дата ремонту: {repair_date.strftime('%d.%m.%Y')}", fontsize=12)

        ax.text(0.02, 0.74, "Вхідні дані:", fontsize=13, fontweight="bold")
        ax.text(0.04, 0.71, f"D = {D:.2f} м;  L = {L:.2f} м;  Wrem = {Wrem:.2f} м;  нахлест = {overlap_cm:.1f} см", fontsize=11)
        ax.text(0.04, 0.68, f"Діаметр електрода = {electrode_diam_mm:.1f} мм;  питома витрата = {spec_consumption:.2f} кг/м;  пачка = {pack_mass:.2f} кг", fontsize=11)
        ax.text(0.04, 0.65, f"Врахування швів приєднання днищ до циліндра: {'Так' if include_shell_to_bottom else 'Ні'}", fontsize=11)

        ax.text(0.02, 0.58, "Ключові результати:", fontsize=13, fontweight="bold")
        ax.text(0.04, 0.55, f"Загальна довжина швів: {weld_total:.2f} м", fontsize=12)
        ax.text(0.04, 0.52, f"Фактична площа матеріалу (з нахлестом): {used_material_area:.3f} м²", fontsize=12)
        ax.text(0.04, 0.49, f"Оцінка листів 6×1.5 м: {total_sheets} шт", fontsize=12)
        ax.text(0.04, 0.46, f"Електроди: {electrodes_mass:.2f} кг  ≈  {packs_needed} пач. по {pack_mass:.2f} кг", fontsize=12)

        ax.text(0.02, 0.08,
                "Примітка: розрахунок витрати електродів є оціночним і базується на заданій питомій витраті (кг/м).\n"
                "Для технологічної карти (WPS) рекомендується узгодження з проектом ремонту/технологом.",
                fontsize=9)

        pdf.savefig(fig1)
        plt.close(fig1)

        # --- PAGE 2: графіки ---
        fig2 = plt.figure(figsize=(11.69, 8.27))  # A4 landscape
        gs = gridspec.GridSpec(1, 2, figure=fig2, wspace=0.20, left=0.05, right=0.98, top=0.88, bottom=0.08)
        fig2.suptitle("Графічні схеми ремонту", fontsize=16, fontweight="bold", y=0.95)

        ax21 = fig2.add_subplot(gs[0, 0])
        ax22 = fig2.add_subplot(gs[0, 1])

        # переносимо зображення з готових фігур
        fig_bot.canvas.draw()
        fig_cyl.canvas.draw()

        ax21.imshow(fig_bot.canvas.buffer_rgba())
        ax22.imshow(fig_cyl.canvas.buffer_rgba())
        ax21.axis("off")
        ax22.axis("off")

        pdf.savefig(fig2)
        plt.close(fig2)

        # --- PAGE 3: 2 таблиці (landscape) ---
        fig3 = plt.figure(figsize=(11.69, 8.27))
        fig3.suptitle("Розрахунок ремонту та швів", fontsize=16, fontweight="bold", y=0.96)

        gs3 = gridspec.GridSpec(1, 2, figure=fig3, wspace=0.12, left=0.04, right=0.98, top=0.88, bottom=0.08)
        ax31 = fig3.add_subplot(gs3[0, 0])
        ax32 = fig3.add_subplot(gs3[0, 1])

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
            ["Циліндр (латки): горизонтальні", fmt_m(weld_cyl_h)],
            ["Циліндр (латки): вертикальні", fmt_m(weld_cyl_v)],
            ["Циліндр (латки): разом", fmt_m(weld_cyl_total)],
            ["Днища (латки): разом", fmt_m(weld_bottoms_total)],
            ["Приєднання днищ до циліндра", fmt_m(weld_shell_to_bottom)],
            ["Загальна довжина швів", fmt_m(weld_total)],
        ]

        make_table_ax(ax31, "1) Дані по ремонту", repair_rows, font_size=10, title_pad=16)
        make_table_ax(ax32, "2) Довжина зварних швів (деталізація)", weld_rows, font_size=10, title_pad=16)

        pdf.savefig(fig3)
        plt.close(fig3)

        # --- PAGE 4: 2 таблиці зверху + СМУГИ (з колишньої стор.5) внизу ---
        fig4 = plt.figure(figsize=(11.69, 8.27))
        fig4.suptitle("Матеріали та електроди", fontsize=16, fontweight="bold", y=0.965)

        gs4 = gridspec.GridSpec(2, 2, figure=fig4,
                                height_ratios=[1.0, 1.05],
                                wspace=0.12, hspace=0.28,
                                left=0.04, right=0.98, top=0.89, bottom=0.06)

        ax41 = fig4.add_subplot(gs4[0, 0])
        ax42 = fig4.add_subplot(gs4[0, 1])
        ax43 = fig4.add_subplot(gs4[1, :])  # на всю ширину

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
            ["Загальна довжина швів", fmt_m(weld_total)],
            ["Потрібна маса електродів", fmt_kg(electrodes_mass)],
            ["Маса 1 пачки", fmt_kg(pack_mass)],
            ["Кількість пачок", f"{packs_needed} шт"],
        ]

        make_table_ax(ax41, "3) Площі та матеріали", areas_rows, font_size=10, title_pad=16)
        make_table_ax(ax42, "4) Електроди (оцінка)", electrodes_rows, font_size=10, title_pad=16)

        # Таблиця "Смуги" — переносимо сюди
        strips_rows = [[k.replace("×", "x"), f"{v} шт"] for k, v in strips_items]
        # щоб не була надто довга, збільшуємо читабельність
        make_table_ax(ax43, "5) Відомість смуг (розкрій/кількість)", strips_rows,
                      col_labels=("Розмір смуги", "К-сть"), font_size=9, title_pad=14, scale_y=1.25)

        pdf.savefig(fig4)
        plt.close(fig4)

    # -----------------------------
    # Download
    # -----------------------------
    st.markdown("---")
    st.markdown("### Збереження результатів")

    st.download_button(
        label="Зберегти звіт як PDF",
        data=pdf_buffer.getvalue(),
        file_name="резервуар_звіт.pdf",
        mime="application/pdf"
    )

    # PNG окремо
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
