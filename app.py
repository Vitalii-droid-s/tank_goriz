import streamlit as st
import matplotlib.pyplot as plt
import math
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages
from datetime import date

VERSION = "1.4"

# =========================
# Matplotlib: кирилиця + стабільний PDF
# =========================
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False

# =========================
# Streamlit UI
# =========================
st.set_page_config(layout="wide")
st.title(f"Розрахунок смуг для резервуара — v{VERSION}")

# ---- Шапка звіту (для клієнта)
h1, h2, h3, h4 = st.columns([1.3, 1.0, 0.8, 0.9])
with h1:
    tank_type = st.text_input("Тип резервуара", value="Паливний резервуар")
with h2:
    serial_no = st.text_input("Заводський №", value="")
with h3:
    repair_date = st.date_input("Дата виконання ремонту", value=date.today())
with h4:
    include_shell_to_bottom = st.checkbox("Враховувати шви приєднання днищ до циліндра", value=True)

st.markdown("---")

# ---- Вхідні геометричні дані
c1, c2, c3, c4 = st.columns(4)
with c1:
    D = st.number_input("Діаметр резервуара, м", value=4.0, min_value=0.1, step=0.1, format="%.1f")
with c2:
    L = st.number_input("Довжина резервуара, м", value=1.0, min_value=0.1, step=0.1, format="%.1f")
with c3:
    Wrem = st.number_input("Ширина ремонтованої ділянки, м", value=1.0, min_value=0.1, step=0.1, format="%.1f")
with c4:
    overlap_cm = st.number_input("Нахлест, см", value=5.0, min_value=0.0, step=0.5, format="%.1f")

# ---- Електроди (прозора модель)
st.markdown("### Електроди (розрахунок витрати)")
e1, e2, e3 = st.columns(3)
with e1:
    electrode_d_mm = st.number_input("Діаметр електрода, мм", value=3.2, min_value=1.6, max_value=6.0, step=0.1, format="%.1f")
with e2:
    pack_mass_kg = st.number_input("Маса 1 пачки, кг", value=2.5, min_value=0.5, step=0.1, format="%.1f")
with e3:
    # За твоїм прикладом: 20 м шва -> 2 пачки по 2.5 кг => 5/20=0.25 кг/м
    specific_kg_per_m = st.number_input("Питома витрата, кг/м (налаштовується)", value=0.25, min_value=0.01, step=0.01, format="%.2f")

st.markdown("---")


# =========================
# Допоміжні функції
# =========================
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
    eps = 1e-6
    if abs(Wrem - 1.0) < eps:
        return [[1.0]]
    if abs(Wrem - 2.0) < eps:
        return [[2.0]]
    if abs(Wrem - 3.0) < eps:
        return [[3.0]]
    if abs(Wrem - 4.0) < eps:
        return [[3.0, 1.0], [1.0, 3.0]]
    if abs(Wrem - 5.0) < eps:
        return [[3.0, 2.0], [2.0, 3.0]]
    if abs(Wrem - 6.0) < eps:
        return [[1.0, 3.0, 2.0], [2.0, 3.0, 1.0]]

    # Довільний Wrem: ріжемо на сегменти по 3 м + залишок
    segs = []
    remain = Wrem
    while remain > 3.0 + 1e-9:
        segs.append(3.0)
        remain -= 3.0
    if remain > 1e-9:
        segs.append(round(remain, 2))
    return [segs]

def fig_table(ax, title, rows, col_widths=(0.62, 0.38), font_size=10):
    """
    Малює табличку на переданому ax (без кириличних проблем).
    rows: list of [key, value]
    """
    ax.axis("off")
    ax.text(0.0, 1.02, title, fontsize=13, fontweight="bold", va="bottom")
    tbl = ax.table(
        cellText=rows,
        colLabels=["Параметр", "Значення"],
        colLoc="left",
        cellLoc="left",
        colWidths=list(col_widths),
        loc="upper left",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(font_size)
    tbl.scale(1, 1.45)

    # стиль заголовка
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.6)
        cell.set_edgecolor("#9aa0a6")
        if r == 0:
            cell.set_facecolor("#eef2f6")
            cell.set_text_props(weight="bold")

def fig_table_simple(ax, title, header, rows, col_widths, font_size=10):
    """
    Таблиця з довільними колонками (для смуг/електродів).
    header: list[str]
    rows: list[list[str]]
    """
    ax.axis("off")
    ax.text(0.0, 1.02, title, fontsize=13, fontweight="bold", va="bottom")
    tbl = ax.table(
        cellText=rows,
        colLabels=header,
        colLoc="left",
        cellLoc="left",
        colWidths=col_widths,
        loc="upper left",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(font_size)
    tbl.scale(1, 1.45)

    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.6)
        cell.set_edgecolor("#9aa0a6")
        if r == 0:
            cell.set_facecolor("#eef2f6")
            cell.set_text_props(weight="bold")


# =========================
# РОЗРАХУНОК
# =========================
if st.button("Розрахувати"):
    if D <= 0 or L <= 0 or Wrem <= 0 or overlap_cm < 0:
        st.error("Введіть додатні D, L, Wrem та невід'ємний нахлест.")
        st.stop()

    overlap = overlap_cm / 100.0
    R = D / 2.0
    h_smuha = 0.5
    circumference = 2.0 * math.pi * R

    # ---- Днище: висота ремонту
    alpha = (Wrem / 2.0) / R
    Hcrit = R * (1.0 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))
    y_cut = -R + Hcrit

    smuhaDict = {}

    # =========================
    # 1) ДНИЩЕ (графік + площа)
    # =========================
    widths_bot, used_heights_bot, areas_bot = [], [], []

    fig_bot, ax_bot = plt.subplots(figsize=(6, 6))
    ax_bot.set_aspect("equal")
    ax_bot.set_title("Днище резервуара (смуги та зона нахлесту)", fontsize=12, fontweight="bold")

    circle = plt.Circle((0, 0), R, edgecolor="black", facecolor="lightyellow", alpha=0.25, zorder=0)
    ax_bot.add_patch(circle)

    ax_bot.axvline(0, color="red", linestyle="--", linewidth=1, alpha=0.7)
    ax_bot.text(
        R * 0.55, -R + Hcrit + 0.03,
        f"H = {Hcrit:.3f} м",
        fontsize=10,
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

        # ширина смуги
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

        width_rnd = round(width, 2)
        key_bot = f"{width_rnd:>5.2f} м × {h_smuha:>3.2f} м"
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2  # симетрично

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

    cum_area_bot = sum(areas_bot)
    total_area_both_bottoms = 2.0 * cum_area_bot

    margin = R * 0.1
    ax_bot.set_xlim(-R - margin, R + margin)
    ax_bot.set_ylim(-R - margin, R + margin)
    ax_bot.set_xlabel("x, м")
    ax_bot.set_ylabel("y, м")
    ax_bot.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # =========================
    # 2) ЦИЛІНДР (графік + матеріал з нахлестом)
    # =========================
    full_rows = 1
    covered_height = h_smuha
    while covered_height + 1e-9 < L:
        full_rows += 1
        covered_height += (h_smuha - overlap)

    patterns = build_patterns(Wrem)

    fig_cyl, ax_cyl = plt.subplots(figsize=(10, 4))
    ax_cyl.set_aspect("equal", adjustable="box")
    ax_cyl.set_title("Розгорнута поверхня (масштаб 1:1, зона ремонту)", fontsize=12, fontweight="bold")

    # показуємо НЕ всю окружність, а зону ремонту з полем — щоб клієнт не думав що є «обман»
    pad_x = max(0.6, 0.15 * Wrem)
    x_view_left = -Wrem / 2.0 - pad_x
    x_view_right = Wrem / 2.0 + pad_x

    ax_cyl.add_patch(plt.Rectangle((x_view_left, 0), (x_view_right - x_view_left), L,
                                   edgecolor="black", facecolor="#d0d0d0", linewidth=1, zorder=1))

    ax_cyl.axvline(0, color="red", linestyle="--", linewidth=1, alpha=0.7)
    ax_cyl.axhline(L, color="red", linestyle="--", linewidth=1, alpha=0.7)

    x_start = -Wrem / 2.0
    area_cyl_with_overlap = 0.0

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

        for seg_idx, seg in enumerate(pat):
            base_h = h_smuha if is_last_row else visible_height

            ax_cyl.add_patch(plt.Rectangle((x_off, y_off), seg, base_h,
                                           edgecolor="black",
                                           facecolor=("orange" if rowNum % 2 == 0 else "lightgreen"),
                                           alpha=0.75, zorder=2))

            # боковий нахлест тільки для “шахматки” 4–6
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

            key_cyl = f"{round(seg, 2):>5.2f} м × {h_smuha:>3.2f} м"
            smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1

            x_off += seg

    ax_cyl.set_xlim(x_view_left, x_view_right)
    ax_cyl.set_ylim(0, max(L, full_rows * (h_smuha - overlap) + overlap))
    ax_cyl.set_xlabel("довжина поверхні (м)")
    ax_cyl.set_ylabel("висота (м)")
    ax_cyl.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # =========================
    # 3) ПЛОЩІ / ШВИ / МАТЕРІАЛИ
    # =========================
    area_cyl = Wrem * L
    real_repair_area = total_area_both_bottoms + area_cyl

    # ---- Площа матеріалу днищ з нахлестом (як було у твоїй логіці)
    bottoms_mat_area = 2.0 * sum(
        w * (h + (overlap / 2.0 if idx < n_bot - 1 else 0.0))
        for idx, (w, h) in enumerate(zip(widths_bot, used_heights_bot))
    )
    total_mat_area = bottoms_mat_area + area_cyl_with_overlap

    # ---- Листи 6×1.5
    sheet_area = 6.0 * 1.5
    total_sheets = math.ceil(total_mat_area / sheet_area)

    # ---- Шви циліндра (облік “ліній”)
    # Горизонтальні: верх+низ кожного ряду = (full_rows+1) * Wrem
    weld_h_cyl = (full_rows + 1) * Wrem

    def row_visible_height(i: int) -> float:
        y = i * (h_smuha - overlap)
        if y >= L:
            return 0.0
        return min(h_smuha, L - y)

    # Вертикальні: 2 крайові + внутрішні стики між сегментами
    weld_v_cyl = 0.0
    for i in range(full_rows):
        h_eff = row_visible_height(i)
        if h_eff <= 0:
            break
        pat = patterns[i % len(patterns)]
        n_vertical = (len(pat) - 1) + 2
        weld_v_cyl += n_vertical * h_eff

    weld_cyl_patch = weld_h_cyl + weld_v_cyl

    # ---- Днища (по твоїй моделі)
    levels = []
    k = 0
    while True:
        y_level = (-R) + k * h_smuha
        if y_level >= y_cut - 1e-9:
            break
        levels.append(y_level)
        k += 1
    levels.append(y_cut)

    weld_h_bot_one = sum(chord_len(R, y) for y in levels)
    theta_top = math.asin(clamp(y_cut / R, -1.0, 1.0))
    arc_len_one_side = R * (theta_top - (-math.pi / 2.0))
    weld_side_bot_one = 2.0 * arc_len_one_side
    weld_bot_one = weld_h_bot_one + weld_side_bot_one
    weld_bottoms_patch = 2.0 * weld_bot_one

    # ---- Додатково: шви приєднання днищ до циліндра (2 кільцеві шви)
    weld_shell_to_bottom = 0.0
    if include_shell_to_bottom:
        weld_shell_to_bottom = 2.0 * circumference

    total_weld = weld_cyl_patch + weld_bottoms_patch + weld_shell_to_bottom

    # =========================
    # 4) Електроди
    # =========================
    electrode_mass_total = total_weld * specific_kg_per_m
    packs_needed = math.ceil(electrode_mass_total / pack_mass_kg)

    # =========================
    # ВИВІД В APP
    # =========================
    st.markdown("## Результати")
    left, right = st.columns([1.2, 1.0])
    with left:
        st.markdown("### Графіки")
        st.pyplot(fig_bot, clear_figure=False)
        st.pyplot(fig_cyl, clear_figure=False)

    # Табличні дані в UI
    repair_table = [
        ["Тип резервуара", tank_type],
        ["Заводський №", serial_no if serial_no else "—"],
        ["Дата ремонту", repair_date.strftime("%d.%m.%Y")],
        ["Діаметр D", f"{D:.2f} м"],
        ["Радіус R", f"{R:.2f} м"],
        ["Довжина L", f"{L:.2f} м"],
        ["Ширина ремонту Wrem", f"{Wrem:.2f} м"],
        ["Нахлест", f"{overlap_cm:.1f} см"],
        ["Висота ремонту днища H", f"{Hcrit:.3f} м"],
        ["Кількість рядів на циліндрі", f"{full_rows} шт"],
    ]

    weld_table = [
        ["Циліндр (латки): горизонтальні", f"{weld_h_cyl:.2f} м"],
        ["Циліндр (латки): вертикальні", f"{weld_v_cyl:.2f} м"],
        ["Циліндр (латки): разом", f"{weld_cyl_patch:.2f} м"],
        ["Днища (латки): разом", f"{weld_bottoms_patch:.2f} м"],
        ["Приєднання днищ до циліндра", f"{weld_shell_to_bottom:.2f} м"],
        ["Загальна довжина швів", f"{total_weld:.2f} м"],
    ]

    area_table = [
        ["Площа 1 днища", f"{cum_area_bot:.3f} м²"],
        ["Площа 2 днищ", f"{total_area_both_bottoms:.3f} м²"],
        ["Площа циліндричної частини (ремонт)", f"{area_cyl:.3f} м²"],
        ["Реально ремонтована площа", f"{real_repair_area:.3f} м²"],
        ["Фактична площа матеріалу (з нахлестом)", f"{total_mat_area:.3f} м²"],
        ["Листи 6×1.5 м (оцінка)", f"{total_sheets} шт"],
    ]

    # смуги
    items = sorted(smuhaDict.items())
    strips_rows = [[k, str(v)] for k, v in items]

    electrode_table = [
        ["Діаметр електрода", f"{electrode_d_mm:.1f} мм"],
        ["Питома витрата", f"{specific_kg_per_m:.2f} кг/м"],
        ["Загальна довжина швів", f"{total_weld:.2f} м"],
        ["Потрібна маса електродів", f"{electrode_mass_total:.2f} кг"],
        ["Маса 1 пачки", f"{pack_mass_kg:.2f} кг"],
        ["Кількість пачок", f"{packs_needed} шт"],
    ]

    with right:
        st.markdown("### Підсумки (таблично)")
        st.table(repair_table)
        st.table(weld_table)
        st.table(area_table)
        st.markdown("### Смуги матеріалу")
        st.table(strips_rows if strips_rows else [["—", "—"]])
        st.markdown("### Електроди")
        st.table(electrode_table)

    # =========================
    # ПРОФЕСІЙНИЙ PDF (A4 Landscape, 2 таблиці/сторінка)
    # =========================
    pdf_buffer = BytesIO()
    with PdfPages(pdf_buffer) as pdf:
        # --- Сторінка 1: Титул
        fig1 = plt.figure(figsize=(11.69, 8.27))  # A4 landscape
        ax = fig1.add_axes([0, 0, 1, 1])
        ax.axis("off")

        ax.text(0.05, 0.92, "ЗВІТ РОЗРАХУНКУ РЕМОНТНИХ СМУГ", fontsize=20, fontweight="bold")
        ax.text(0.05, 0.87, f"Версія розрахунку: v{VERSION}", fontsize=11)
        ax.text(0.05, 0.82, f"Тип резервуара: {tank_type}", fontsize=13)
        ax.text(0.05, 0.78, f"Заводський №: {serial_no if serial_no else '—'}", fontsize=13)
        ax.text(0.05, 0.74, f"Дата ремонту: {repair_date.strftime('%d.%m.%Y')}", fontsize=13)

        ax.text(0.05, 0.66, "Вхідні дані:", fontsize=14, fontweight="bold")
        inputs_text = (
            f"D = {D:.2f} м;  L = {L:.2f} м;  Wrem = {Wrem:.2f} м;  нахлест = {overlap_cm:.1f} см\n"
            f"Діаметр електрода = {electrode_d_mm:.1f} мм;  питома витрата = {specific_kg_per_m:.2f} кг/м;  пачка = {pack_mass_kg:.2f} кг\n"
            f"Врахування швів приєднання днищ до циліндра: {'Так' if include_shell_to_bottom else 'Ні'}"
        )
        ax.text(0.05, 0.60, inputs_text, fontsize=12)

        ax.text(0.05, 0.48, "Ключові результати:", fontsize=14, fontweight="bold")
        kpi = (
            f"Загальна довжина швів: {total_weld:.2f} м\n"
            f"Фактична площа матеріалу (з нахлестом): {total_mat_area:.3f} м²\n"
            f"Оцінка листів 6×1.5 м: {total_sheets} шт\n"
            f"Електроди: {electrode_mass_total:.2f} кг ≈ {packs_needed} пач. по {pack_mass_kg:.2f} кг"
        )
        ax.text(0.05, 0.40, kpi, fontsize=13)

        ax.text(
            0.05, 0.18,
            "Примітка: розрахунок витрати електродів є оціночним і базується на вказаній питомій витраті (кг/м).\n"
            "Для повної технологічної карти витрати рекомендується узгодження з WPS/проектом ремонту.",
            fontsize=10
        )

        pdf.savefig(fig1, bbox_inches="tight")
        plt.close(fig1)

        # --- Сторінка 2: графіки поруч
        fig2 = plt.figure(figsize=(11.69, 8.27))
        gs = fig2.add_gridspec(1, 2, left=0.04, right=0.98, top=0.92, bottom=0.06, wspace=0.10)

        axA = fig2.add_subplot(gs[0, 0])
        axB = fig2.add_subplot(gs[0, 1])
        axA.axis("off"); axB.axis("off")

        # вставляємо готові фігури як зображення
        bufA = BytesIO(); fig_bot.savefig(bufA, format="png", dpi=220, bbox_inches="tight"); bufA.seek(0)
        bufB = BytesIO(); fig_cyl.savefig(bufB, format="png", dpi=220, bbox_inches="tight"); bufB.seek(0)

        imgA = plt.imread(bufA)
        imgB = plt.imread(bufB)
        axA.imshow(imgA); axB.imshow(imgB)

        fig2.suptitle("Графічні схеми ремонту", fontsize=16, fontweight="bold", y=0.97)
        pdf.savefig(fig2, bbox_inches="tight")
        plt.close(fig2)

        # --- Сторінка 3: 2 таблиці (ремонт + шви)
        fig3 = plt.figure(figsize=(11.69, 8.27))
        gs3 = fig3.add_gridspec(1, 2, left=0.04, right=0.98, top=0.92, bottom=0.06, wspace=0.12)

        ax31 = fig3.add_subplot(gs3[0, 0])
        ax32 = fig3.add_subplot(gs3[0, 1])

        fig_table(ax31, "1) Дані по ремонту", repair_table, col_widths=(0.60, 0.40), font_size=10)
        fig_table(ax32, "2) Довжина зварних швів (деталізація)", weld_table, col_widths=(0.68, 0.32), font_size=10)

        fig3.suptitle("Розрахунок ремонту та швів", fontsize=16, fontweight="bold", y=0.97)
        pdf.savefig(fig3, bbox_inches="tight")
        plt.close(fig3)

        # --- Сторінка 4: 2 таблиці (площі/матеріали + електроди)
        fig4 = plt.figure(figsize=(11.69, 8.27))
        gs4 = fig4.add_gridspec(1, 2, left=0.04, right=0.98, top=0.92, bottom=0.06, wspace=0.12)

        ax41 = fig4.add_subplot(gs4[0, 0])
        ax42 = fig4.add_subplot(gs4[0, 1])

        fig_table(ax41, "3) Площі та матеріали", area_table, col_widths=(0.70, 0.30), font_size=10)
        fig_table(ax42, "4) Електроди (оцінка)", electrode_table, col_widths=(0.66, 0.34), font_size=10)

        fig4.suptitle("Матеріали та електроди", fontsize=16, fontweight="bold", y=0.97)
        pdf.savefig(fig4, bbox_inches="tight")
        plt.close(fig4)

        # --- Сторінка 5: смуги (може бути багато рядків)
        fig5 = plt.figure(figsize=(11.69, 8.27))
        ax5 = fig5.add_axes([0.04, 0.06, 0.94, 0.86])

        header = ["Розмір смуги", "К-сть, шт"]
        rows = strips_rows if strips_rows else [["—", "—"]]

        # підбираємо ширини
        fig_table_simple(ax5, "5) Специфікація смуг матеріалу", header, rows,
                         col_widths=[0.75, 0.25], font_size=10)

        fig5.suptitle("Специфікація матеріалу", fontsize=16, fontweight="bold", y=0.97)
        pdf.savefig(fig5, bbox_inches="tight")
        plt.close(fig5)

    # =========================
    # Завантаження PDF
    # =========================
    st.markdown("---")
    st.markdown("### Експорт")
    st.download_button(
        label="Зберегти звіт як PDF",
        data=pdf_buffer.getvalue(),
        file_name="резервуар_звіт.pdf",
        mime="application/pdf"
    )

    # clean
    plt.close(fig_bot)
    plt.close(fig_cyl)
