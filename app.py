import streamlit as st
import matplotlib.pyplot as plt
import math
from io import BytesIO
from datetime import date
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# PDF (професійні таблиці)
from reportlab.lib.pagesizes import A4, landscape
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import mm

VERSION = "1.4"

# -----------------------------
# Streamlit
# -----------------------------
st.set_page_config(layout="wide")
st.title(f"Розрахунок смуг для резервуара — v{VERSION}")

# -----------------------------
# Допоміжні
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
    return 2.0 * math.sqrt(max(0.0, R * R - y * y))

def build_patterns(Wrem: float):
    """
    - Для Wrem 4..6: залишаємо «шахматку» (як було)
    - Для будь-якого іншого Wrem: ріжемо на 3 м + залишок
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

def fig_to_png_bytes(fig, dpi=200):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()

def make_table(data, col_widths=None, title=None):
    """
    data: list[list[str]]
    """
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "tbl_title",
        parent=styles["Heading4"],
        spaceAfter=6,
        spaceBefore=0
    )

    elements = []
    if title:
        elements.append(Paragraph(title, title_style))

    tbl = Table(data, colWidths=col_widths, hAlign="LEFT")
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f1f3f5")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#111827")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 10),

        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 1), (-1, -1), 9),

        ("ALIGN", (0, 0), (-1, 0), "CENTER"),
        ("ALIGN", (0, 1), (0, -1), "LEFT"),
        ("ALIGN", (1, 1), (-1, -1), "LEFT"),

        ("GRID", (0, 0), (-1, -1), 0.6, colors.HexColor("#9ca3af")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#fafafa")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    elements.append(tbl)
    return elements

def build_pdf_report(
    tank_type, serial_no, repair_date,
    fig_bot, fig_cyl,
    tbl_job_1, tbl_job_2,   # 2 таблиці на сторінку 1
    tbl_geom_1, tbl_geom_2, # 2 таблиці на сторінку 2
    tbl_weld_1, tbl_weld_2, # 2 таблиці на сторінку 3
    tbl_mat_1, tbl_mat_2,   # 2 таблиці на сторінку 4
):
    styles = getSampleStyleSheet()
    h1 = ParagraphStyle("h1", parent=styles["Heading1"], fontSize=16, spaceAfter=6)
    meta = ParagraphStyle("meta", parent=styles["Normal"], fontSize=10, textColor=colors.HexColor("#374151"))
    note = ParagraphStyle("note", parent=styles["Normal"], fontSize=9, textColor=colors.HexColor("#6b7280"))

    buf = BytesIO()
    doc = SimpleDocTemplate(
        buf,
        pagesize=landscape(A4),
        leftMargin=14*mm, rightMargin=14*mm,
        topMargin=12*mm, bottomMargin=12*mm
    )

    elements = []

    # ШАПКА
    elements.append(Paragraph("ЗВІТ РОЗРАХУНКУ РЕМОНТУ РЕЗЕРВУАРА (СМУГИ, ШВИ, МАТЕРІАЛИ)", h1))
    elements.append(Paragraph(
        f"<b>Тип резервуара:</b> {tank_type or '—'} &nbsp;&nbsp;&nbsp; "
        f"<b>Зав. №:</b> {serial_no or '—'} &nbsp;&nbsp;&nbsp; "
        f"<b>Дата виконання ремонту:</b> {repair_date.strftime('%d.%m.%Y')}",
        meta
    ))
    elements.append(Spacer(1, 6))

    # Сторінка з графіками (2 в ряд)
    img_bot = Image(BytesIO(fig_to_png_bytes(fig_bot, dpi=220)))
    img_cyl = Image(BytesIO(fig_to_png_bytes(fig_cyl, dpi=220)))

    # Підігнати під landscape A4
    img_bot.drawWidth = 130*mm
    img_bot.drawHeight = 130*mm
    img_cyl.drawWidth = 210*mm
    img_cyl.drawHeight = 110*mm

    elements.append(Paragraph("Схеми викладки (візуалізація)", styles["Heading3"]))
    elements.append(Spacer(1, 6))

    gtbl = Table([[img_bot, img_cyl]], colWidths=[135*mm, 250*mm])
    gtbl.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    elements.append(gtbl)
    elements.append(Spacer(1, 6))
    elements.append(Paragraph(
        "Примітка: масштаб відображення 1:1 забезпечується однаковим співвідношенням одиниць по осях (equal aspect).",
        note
    ))
    elements.append(PageBreak())

    # Далі — сторінки таблиць по 2 на сторінку
    def two_tables_page(title, t1, t2):
        elements.append(Paragraph(title, styles["Heading3"]))
        elements.append(Spacer(1, 6))
        t = Table([[t1, t2]], colWidths=[185*mm, 185*mm])
        t.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        elements.append(t)
        elements.append(PageBreak())

    two_tables_page("1) Дані ремонту та налаштування розрахунку", tbl_job_1, tbl_job_2)
    two_tables_page("2) Геометрія та площі", tbl_geom_1, tbl_geom_2)
    two_tables_page("3) Шви та наплавлення", tbl_weld_1, tbl_weld_2)
    two_tables_page("4) Матеріали та облік смуг", tbl_mat_1, tbl_mat_2)

    doc.build(elements)
    buf.seek(0)
    return buf.getvalue()

# -----------------------------
# Вхідні дані (додали шапку звіту)
# -----------------------------
st.markdown("### Дані для звіту")
c0, c1, c2, c3 = st.columns([1.2, 1, 1, 1])
with c0:
    tank_type = st.text_input("Тип резервуара", value="")
with c1:
    serial_no = st.text_input("Заводський №", value="")
with c2:
    repair_date = st.date_input("Дата виконання ремонту", value=date.today())
with c3:
    electrode_d_mm = st.number_input("Діаметр електрода, мм", value=3.2, min_value=1.6, step=0.1, format="%.1f")

st.markdown("---")
st.markdown("### Геометрія та технологія")

col1, col2, col3, col4 = st.columns(4)
with col1:
    D = st.number_input("Діаметр резервуара, м", value=2.5, min_value=0.1, step=0.1, format="%.1f")
with col2:
    L = st.number_input("Довжина резервуара, м", value=4.3, min_value=0.1, step=0.1, format="%.1f")
with col3:
    Wrem = st.number_input("Ширина ремонтованої ділянки, м", value=1.0, min_value=0.1, step=0.1, format="%.1f")
with col4:
    overlap_cm = st.number_input("Нахлест, см", value=5.0, min_value=0.0, step=0.5, format="%.1f")

col5, col6, col7, col8 = st.columns(4)
with col5:
    a_mm = st.number_input("Катет шва a, мм", value=4.0, min_value=1.0, step=0.5, format="%.1f")
with col6:
    eta = st.number_input("ККД наплавлення η", value=0.65, min_value=0.10, max_value=0.95, step=0.01, format="%.2f")
with col7:
    K_loss = st.number_input("Коеф. втрат K", value=1.10, min_value=1.00, max_value=2.00, step=0.01, format="%.2f")
with col8:
    pack_mass_kg = st.number_input("Маса 1 пачки, кг", value=2.5, min_value=0.5, step=0.5, format="%.1f")

col9, col10 = st.columns(2)
with col9:
    double_lap_line = st.checkbox("Подвійна лінія нахльосту (рахувати 2 шви по нахльосту)", value=True)
with col10:
    add_cyl_to_bottoms = st.checkbox("Додати приварку циліндричної частини до днищ (2×коло)", value=True)

# -----------------------------
# Розрахунок
# -----------------------------
if st.button("Розрахувати"):
    if D <= 0 or L <= 0 or Wrem <= 0 or overlap_cm < 0:
        st.error("Введіть додатні D, L, Wrem та невід'ємний нахлест.")
        st.stop()

    overlap = overlap_cm / 100.0  # м
    R = D / 2.0
    h_smuha = 0.5
    rho = 7850.0  # кг/м³ (сталь, для оцінки)

    circumference = 2.0 * math.pi * R
    if Wrem > circumference + 1e-9:
        st.warning(f"Wrem = {Wrem:.2f} м більша за довжину кола {circumference:.2f} м. Розгортка умовна.")

    # Висота ремонту по днищу (критична)
    alpha = (Wrem / 2.0) / R
    Hcrit = R * (1.0 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))

    # ============================================
    # 1) ДНИЩЕ (одне)
    # ============================================
    widths_bot = []
    used_heights_bot = []
    areas_bot = []
    smuhaDict = {}

    fig_bot, ax_bot = plt.subplots(figsize=(6.3, 6.3))
    ax_bot.set_aspect("equal")
    ax_bot.set_title("Днище (смуги та зона нахлесту)", fontsize=13, fontweight="bold")

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

        real_h = h_smuha
        if j == n_bot - 1:
            real_h = max(0.0, min(y_cut - y_bot, h_smuha))

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
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2  # два симетричних

        x_left = -width / 2.0

        ax_bot.add_patch(plt.Rectangle(
            (x_left, y_bot), width, h_smuha,
            edgecolor="black", facecolor="skyblue", alpha=0.45, zorder=2
        ))

        ov_h = h_smuha + (overlap/2.0 if j < n_bot - 1 else 0.0)
        ax_bot.add_patch(plt.Rectangle(
            (x_left, y_bot), width, ov_h,
            edgecolor="red", facecolor="none",
            linestyle="--", linewidth=1, alpha=0.7, zorder=3
        ))

        if real_h > 0:
            ax_bot.add_patch(plt.Rectangle(
                (x_left, y_bot), width, real_h,
                edgecolor="black", facecolor="skyblue", alpha=0.55, zorder=4
            ))

        extra_h = h_smuha - real_h
        if extra_h > 1e-9:
            ax_bot.add_patch(plt.Rectangle(
                (x_left, y_bot + real_h), width, extra_h,
                edgecolor="red", facecolor="none", hatch="///", alpha=0.5, zorder=5
            ))

        ax_bot.text(0, y_bot + h_smuha/2.0, f"S{j+1}\n{areas_bot[-1]:.3f} м²",
                    ha="center", va="center", fontsize=8, zorder=7)

    cum_area_bot = sum(areas_bot)
    total_area_both_bottoms = 2.0 * cum_area_bot

    margin = R * 0.1
    ax_bot.set_xlim(-R - margin, R + margin)
    ax_bot.set_ylim(-R - margin, R + margin)
    ax_bot.set_xlabel("x (м)")
    ax_bot.set_ylabel("y (м)")
    ax_bot.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # ============================================
    # 2) ЦИЛІНДР
    # ============================================
    full_rows = 1
    covered_height = h_smuha
    while covered_height + 1e-9 < L:
        full_rows += 1
        covered_height += (h_smuha - overlap)

    patterns = build_patterns(Wrem)

    fig_cyl, ax_cyl = plt.subplots(figsize=(10.5, 4.2))
    ax_cyl.set_aspect("equal", adjustable="box")  # 1:1
    ax_cyl.set_title("Розгорнута поверхня циліндра (чергування смуг, червоний — зона нахлесту)",
                     fontsize=13, fontweight="bold")

    ax_cyl.add_patch(plt.Rectangle(
        (-circumference/2.0, 0), circumference, L,
        edgecolor="black", facecolor="#d0d0d0", linewidth=1, zorder=1
    ))

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

        max_top_for_ylim = max(max_top_for_ylim, y_off + h_smuha)

        for seg_idx, seg in enumerate(pat):
            base_h = h_smuha if is_last_row else visible_height

            ax_cyl.add_patch(plt.Rectangle(
                (x_off, y_off), seg, base_h,
                edgecolor="black",
                facecolor=("orange" if (rowNum % 2 == 0) else "lightgreen"),
                alpha=0.7, zorder=2
            ))

            # Боковий нахлест — тільки для шахматки 4..6
            if 4.0 - 1e-6 <= Wrem <= 6.0 + 1e-6 and len(pat) > 1:
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
                ax_cyl.add_patch(plt.Rectangle(
                    (ov_x, ov_y_base), ov_w, ov_h,
                    edgecolor="red", facecolor="none",
                    linestyle="--", linewidth=1, alpha=0.7, zorder=3
                ))

            площа_cyl_з_нахлестом += ov_w * ov_h

            ax_cyl.text(x_off + seg/2.0, y_off + max(0.02, base_h/2.0),
                        f"{seg:.2f} м", ha="center", va="center", fontsize=9, zorder=5)

            key_cyl = f"{round(seg, 2):>5.2f} м × {h_smuha:>3.2f} м"
            smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1
            x_off += seg

    ax_cyl.set_xlim(-circumference/2.0, circumference/2.0)
    ax_cyl.set_ylim(0, max_top_for_ylim)
    ax_cyl.set_xlabel("довжина поверхні (м)")
    ax_cyl.set_ylabel("висота (м)")
    ax_cyl.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # ============================================
    # 3) Площі
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
    # 4) ШВИ (логічно і прозоро)
    # ============================================
    # 4.1 Циліндр: вертикальні шви
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
        # 2 крайові + внутрішні межі сегментів
        n_vertical = 2 + max(0, len(pat) - 1)
        weld_v_cyl += n_vertical * h_eff

    # 4.2 Циліндр: горизонтальні (між рядами + крайові)
    if double_lap_line:
        # 2 лінії на кожному нахльості між рядами + 2 крайові (низ/верх) по 1 лінії
        weld_h_cyl = 2.0 * (full_rows - 1) * Wrem + 2.0 * Wrem
    else:
        # як раніше (одна лінія на межу)
        weld_h_cyl = (full_rows + 1) * Wrem

    weld_cyl = weld_h_cyl + weld_v_cyl

    # 4.3 Днища (2 шт): горизонталі по хордах + бокові дуги (ліва+права)
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
    weld_side_bot_one = 2.0 * arc_len_one_side  # 2 боки

    weld_bot_one = weld_h_bot_one + weld_side_bot_one
    weld_both_bottoms = 2.0 * weld_bot_one

    # 4.4 Приварка циліндра до днищ (2×коло)
    weld_cyl_to_bottoms = (2.0 * circumference) if add_cyl_to_bottoms else 0.0

    # Геометрична довжина швів
    total_weld_geom = weld_cyl + weld_both_bottoms + weld_cyl_to_bottoms

    # Технологічна «довжина наплавлення» (для електродів)
    # - якщо double_lap_line=True: додатково врахували 2-гу лінію по нахльосту вже у weld_h_cyl
    # - технологічно для електродів беремо саме total_weld_geom
    total_deposition_length = total_weld_geom

    # ============================================
    # 5) ЕЛЕКТРОДИ (масовий метод)
    # ============================================
    a = a_mm / 1000.0
    area_fillet = (a * a) / 2.0  # м² на 1 м шва (трикутник)
    m_deposit_per_m = rho * area_fillet  # кг/м (метал шва)
    m_electrode_per_m = (m_deposit_per_m / max(1e-9, eta)) * K_loss  # кг/м електродів
    m_total = m_electrode_per_m * total_deposition_length
    packs = math.ceil(m_total / pack_mass_kg)

    # ============================================
    # ВИВІД в Streamlit
    # ============================================
    st.markdown("## Результати")
    cA, cB = st.columns([1.2, 1])
    with cA:
        st.markdown("### Графіки")
        st.pyplot(fig_bot, clear_figure=False)
        st.pyplot(fig_cyl, clear_figure=False)

    with cB:
        st.markdown("### Підсумок (коротко)")
        st.write(f"**Площа ремонту (фактична):** {реально_ремонт_площа:.3f} м²")
        st.write(f"**Площа матеріалу з нахлестами:** {фактична_площа_матеріалів:.3f} м²")
        st.write(f"**Листи 6×1.5 м (орієнтовно):** {total_sheets} шт")
        st.write(f"**Довжина швів (геометрична):** {total_weld_geom:.2f} м")
        st.write(f"**Електроди всього (оцінка):** {m_total:.2f} кг ≈ **{packs} пач(ки)** × {pack_mass_kg:.1f} кг")

    # ============================================
    # PDF: готуємо 8 таблиць (2 на сторінку)
    # ============================================
    # Табл. 1 (дані ремонту)
    job_1 = [
        ["Реквізит", "Значення"],
        ["Тип резервуара", tank_type or "—"],
        ["Заводський №", serial_no or "—"],
        ["Дата виконання ремонту", repair_date.strftime("%d.%m.%Y")],
        ["Версія розрахунку", VERSION],
    ]

    # Табл. 2 (вхідні параметри)
    job_2 = [
        ["Параметр", "Значення"],
        ["Діаметр резервуара D", f"{D:.2f} м"],
        ["Довжина резервуара L", f"{L:.2f} м"],
        ["Ширина ремонту Wrem", f"{Wrem:.2f} м"],
        ["Нахлест", f"{overlap_cm:.1f} см"],
        ["Висота смуги", f"{h_smuha:.2f} м"],
        ["Довжина кола", f"{circumference:.3f} м"],
        ["Рядів на циліндрі", f"{full_rows} шт"],
        ["Смуг на днищі (по висоті)", f"{n_bot} шт"],
    ]

    # Табл. 3 (геометрія/площі)
    geom_1 = [
        ["Показник", "Значення"],
        ["Площа одного днища", f"{cum_area_bot:.3f} м²"],
        ["Площа обох днищ", f"{total_area_both_bottoms:.3f} м²"],
        ["Площа циліндричної частини", f"{площа_cyl:.3f} м²"],
        ["Площа ремонту (фактична)", f"{реально_ремонт_площа:.3f} м²"],
        ["Площа матеріалу (з нахлестами)", f"{фактична_площа_матеріалів:.3f} м²"],
        ["Листи 6×1.5 м (орієнтовно)", f"{total_sheets} шт"],
    ]

    # Табл. 4 (пояснення патерну)
    pat_txt = " / ".join([", ".join([f"{x:.2f}" for x in row]) for row in patterns])
    geom_2 = [
        ["Пояснення", "Значення"],
        ["Сегментація по ширині Wrem", pat_txt],
        ["Примітка", "Для Wrem=4–6 м використана «шахматка», для інших — 3.00 м + залишок."],
    ]

    # Табл. 5 (шви — деталізація)
    weld_1 = [
        ["Показник", "Значення"],
        ["Циліндр: вертикальні шви", f"{weld_v_cyl:.2f} м"],
        ["Циліндр: горизонтальні шви", f"{weld_h_cyl:.2f} м"],
        ["Циліндр: разом", f"{weld_cyl:.2f} м"],
        ["Днища (2 шт): разом", f"{weld_both_bottoms:.2f} м"],
        ["Приварка циліндра до днищ (2×коло)", f"{weld_cyl_to_bottoms:.2f} м"],
        ["Довжина швів (геометрична)", f"{total_weld_geom:.2f} м"],
        ["Довжина наплавлення (для електродів)", f"{total_deposition_length:.2f} м"],
    ]

    # Табл. 6 (електроди — без «прикладу»)
    weld_2 = [
        ["Показник", "Значення"],
        ["Електрод", f"Ø {electrode_d_mm:.1f} мм"],
        ["Катет шва a", f"{a_mm:.1f} мм"],
        ["ККД наплавлення η", f"{eta:.2f}"],
        ["Коеф. втрат K", f"{K_loss:.2f}"],
        ["Витрата електродів на 1 м (оцінка)", f"{m_electrode_per_m:.3f} кг/м"],
        ["Потрібно електродів всього (оцінка)", f"{m_total:.2f} кг"],
        ["Пачки електродів", f"{packs} пач(ки) × {pack_mass_kg:.1f} кг"],
    ]

    # Табл. 7 (матеріали — смуги)
    items = sorted(smuhaDict.items(), key=lambda x: x[0])
    mat_rows = [["Розмір смуги", "К-сть, шт"]]
    for kkey, cnt in items:
        mat_rows.append([kkey, str(cnt)])

    # щоб не було занадто довго — обмежимо в PDF, але в Streamlit покажемо повністю
    mat_1 = mat_rows[:22] if len(mat_rows) > 22 else mat_rows
    if len(mat_rows) > 22:
        mat_1.append(["…", f"ще {len(mat_rows)-22} позицій"])

    # Табл. 8 (підсумок матеріалів/пояснення довіри)
    mat_2 = [
        ["Пункт контролю", "Пояснення"],
        ["Що враховано у швах", "Вертикальні шви смуг, горизонтальні межі рядів, шви днищ та (опц.) 2×коло приварки."],
        ["Нахлест", f"{overlap_cm:.1f} см. Опція «подвійна лінія» = 2 шви по межі нахльосту." if double_lap_line else "Одна лінія шва по межі."],
        ["Електроди", "Оцінка через масу наплавлення: ρ·(a²/2) на 1 м шва, з урахуванням η та K."],
    ]

    # Streamlit: повний список смуг
    st.markdown("---")
    st.markdown("### Смуги матеріалу (повний список)")
    st.dataframe(
        [{"Розмір смуги": k, "К-сть, шт": v} for k, v in items],
        use_container_width=True,
        hide_index=True
    )

    # ============================================
    # PDF генерація (A4 landscape, 2 таблиці/сторінку)
    # ============================================
    # В ReportLab краще фіксовані ширини колонок
    job_1_el = make_table(job_1, col_widths=[70*mm, 110*mm], title="Таблиця 1. Реквізити")
    job_2_el = make_table(job_2, col_widths=[80*mm, 95*mm], title="Таблиця 2. Вхідні параметри")

    geom_1_el = make_table(geom_1, col_widths=[95*mm, 80*mm], title="Таблиця 3. Площі")
    geom_2_el = make_table(geom_2, col_widths=[70*mm, 110*mm], title="Таблиця 4. Логіка сегментації")

    weld_1_el = make_table(weld_1, col_widths=[105*mm, 70*mm], title="Таблиця 5. Дідомка швів")
    weld_2_el = make_table(weld_2, col_widths=[105*mm, 70*mm], title="Таблиця 6. Електроди (оцінка)")

    mat_1_el = make_table(mat_1, col_widths=[115*mm, 55*mm], title="Таблиця 7. Облік смуг (скорочено)")
    mat_2_el = make_table(mat_2, col_widths=[70*mm, 110*mm], title="Таблиця 8. Пояснення для клієнта")

    # Витягуємо об’єкти таблиць (останній елемент у кожному списку — Table)
    tbl_job_1 = job_1_el[-1]
    tbl_job_2 = job_2_el[-1]
    tbl_geom_1 = geom_1_el[-1]
    tbl_geom_2 = geom_2_el[-1]
    tbl_weld_1 = weld_1_el[-1]
    tbl_weld_2 = weld_2_el[-1]
    tbl_mat_1 = mat_1_el[-1]
    tbl_mat_2 = mat_2_el[-1]

    pdf_bytes = build_pdf_report(
        tank_type, serial_no, repair_date,
        fig_bot, fig_cyl,
        tbl_job_1, tbl_job_2,
        tbl_geom_1, tbl_geom_2,
        tbl_weld_1, tbl_weld_2,
        tbl_mat_1, tbl_mat_2
    )

    st.download_button(
        label="Зберегти звіт (PDF, професійний)",
        data=pdf_bytes,
        file_name="резервуар_звіт.pdf",
        mime="application/pdf"
    )

    # PNG
    cc1, cc2 = st.columns(2)
    with cc1:
        st.download_button(
            label="Зберегти днище (PNG)",
            data=fig_to_png_bytes(fig_bot, dpi=220),
            file_name="днище.png",
            mime="image/png"
        )
    with cc2:
        st.download_button(
            label="Зберегти циліндр (PNG)",
            data=fig_to_png_bytes(fig_cyl, dpi=220),
            file_name="циліндр.png",
            mime="image/png"
        )

    plt.close(fig_bot)
    plt.close(fig_cyl)
