import streamlit as st
import matplotlib.pyplot as plt
import math
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages

VERSION = "1.3"

# -----------------------------
# Налаштування Streamlit
# -----------------------------
st.set_page_config(layout="wide")
st.title(f"Розрахунок смуг для резервуара — v{VERSION}")

# -----------------------------
# Допоміжні функції
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
    if y <= -R or y >= R:
        return 0.0
    return 2.0 * math.sqrt(max(0.0, R*R - y*y))

def build_patterns(Wrem: float):
    """
    Повертає список патернів (рядків), де кожен рядок — список сегментів по ширині (м).
    - Для Wrem 1..3 — один сегмент (Wrem)
    - Для 4..6 — “шахматка” як у тебе
    - Для довільних значень — сегментація на 3м + залишок
    """
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

    # Довільний Wrem: ріжемо на 3.0 м + залишок
    segs = []
    remain = Wrem
    while remain > 3.0 + 1e-9:
        segs.append(3.0)
        remain -= 3.0
    if remain > 1e-9:
        segs.append(round(remain, 2))
    return [segs]  # без шахматки


# -----------------------------
# Вхідні дані
# -----------------------------
col1, col2, col3, col4 = st.columns(4)
with col1:
    D = st.number_input("Діаметр резервуара, м", value=2.5, min_value=0.1, step=0.1, format="%.1f")
with col2:
    L = st.number_input("Довжина резервуара, м", value=4.3, min_value=0.1, step=0.1, format="%.1f")
with col3:
    Wrem = st.number_input("Ширина ремонтованої ділянки, м", value=1.0, min_value=0.1, step=0.1, format="%.1f")
with col4:
    overlap_cm = st.number_input("Нахлест, см", value=5.0, min_value=0.0, step=0.5, format="%.1f")


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

    circumference = 2.0 * math.pi * R
    if Wrem > circumference + 1e-9:
        st.warning(f"Wrem = {Wrem:.2f} м більша за довжину кола {circumference:.2f} м. Розгортка буде умовною.")

    # Висота ремонту по днищу (критична)
    alpha = (Wrem / 2.0) / R
    Hcrit = R * (1.0 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))

    # ============================================
    # 1) ДНИЩЕ (одне) — точна площа, симетрія x2
    # ============================================
    widths_bot = []
    used_heights_bot = []
    areas_bot = []
    smuhaDict = {}

    fig_bot, ax_bot = plt.subplots(figsize=(6, 6))
    ax_bot.set_aspect("equal")
    ax_bot.set_title(
        "Днище резервуара (синій — смуга, червоний — зона нахлесту)",
        fontsize=13, fontweight="bold"
    )

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

        # Фактична висота останньої смуги
        real_h = h_smuha
        if j == n_bot - 1:
            real_h = max(0.0, min(y_cut - y_bot, h_smuha))

        # y_ref для ширини (твоя логіка)
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

        # Точна площа смуги в колі (для ОДНОГО днища)
        y_top_clip = min(y_bot + real_h, y_cut)
        area_strip_exact = primitive_F(R, y_top_clip) - primitive_F(R, y_bot)
        areas_bot.append(area_strip_exact)

        # Облік смуг (дві симетричні)
        width_rnd = round(width, 2)
        key_bot = f"{width_rnd:>5.2f}м x {h_smuha:>3.2f}м"
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2

        x_left = -width / 2.0

        # Візуал: базова смуга
        ax_bot.add_patch(plt.Rectangle(
            (x_left, y_bot), width, h_smuha,
            edgecolor="black", facecolor="skyblue", alpha=0.45, zorder=2
        ))

        # Візуал: зона нахлесту (по вертикалі між рядами)
        ov_h = h_smuha + (overlap/2.0 if j < n_bot - 1 else 0.0)
        ax_bot.add_patch(plt.Rectangle(
            (x_left, y_bot), width, ov_h,
            edgecolor="red", facecolor="none",
            linestyle="--", linewidth=1, alpha=0.7, zorder=3
        ))

        # Візуал: реально використана частина (якщо обрізана)
        if real_h > 0:
            ax_bot.add_patch(plt.Rectangle(
                (x_left, y_bot), width, real_h,
                edgecolor="black", facecolor="skyblue", alpha=0.55, zorder=4
            ))

        # Візуал: надлишок верхньої смуги
        extra_h = h_smuha - real_h
        if extra_h > 1e-9:
            ax_bot.add_patch(plt.Rectangle(
                (x_left, y_bot + real_h), width, extra_h,
                edgecolor="red", facecolor="none", hatch="///", alpha=0.5, zorder=5
            ))

        ax_bot.text(
            0, y_bot + h_smuha/2.0,
            f"S{j+1}\n{areas_bot[-1]:.3f} м²",
            ha="center", va="center", fontsize=8, zorder=7
        )

    cum_area_bot = sum(areas_bot)                # одне днище
    total_area_both_bottoms = 2.0 * cum_area_bot # два днища

    margin = R * 0.1
    ax_bot.set_xlim(-R - margin, R + margin)
    ax_bot.set_ylim(-R - margin, R + margin)
    ax_bot.set_xlabel("x (м)")
    ax_bot.set_ylabel("y (м)")
    ax_bot.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # ============================================
    # 2) ЦИЛІНДР — з нахлестом по вертикалі + (опц.) боковий
    # ============================================
    # К-сть рядів з урахуванням нахлесту
    full_rows = 1
    covered_height = h_smuha
    while covered_height + 1e-9 < L:
        full_rows += 1
        covered_height += (h_smuha - overlap)

    patterns = build_patterns(Wrem)

    fig_cyl, ax_cyl = plt.subplots(figsize=(8, 6))
    ax_cyl.set_aspect("auto")
    ax_cyl.set_title(
        "Розгорнута поверхня циліндра (чергування смуг, червоний — зона нахлесту)",
        fontsize=13, fontweight="bold"
    )

    # Корпус
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

        if is_last_row:
            max_top_for_ylim = max(max_top_for_ylim, y_off + h_smuha)
        else:
            max_top_for_ylim = max(max_top_for_ylim, top_ov_end)

        for seg_idx, seg in enumerate(pat):
            base_h = h_smuha if is_last_row else visible_height

            ax_cyl.add_patch(plt.Rectangle(
                (x_off, y_off), seg, base_h,
                edgecolor="black",
                facecolor=("orange" if (rowNum % 2 == 0) else "lightgreen"),
                alpha=0.7, zorder=2
            ))

            # Якщо останній ряд виходить вище L — заштрихувати надлишок
            if is_last_row and (y_off + h_smuha > L):
                ax_cyl.add_patch(plt.Rectangle(
                    (x_off, L), seg, (y_off + h_smuha - L),
                    edgecolor="red", facecolor="none", hatch="///", alpha=0.5, zorder=4
                ))

            # Боковий нахлест — тільки для “шахматних” 4..6 (як у тебе)
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
                ax_cyl.add_patch(plt.Rectangle(
                    (ov_x, ov_y_base), ov_w, ov_h,
                    edgecolor="red", facecolor="none",
                    linestyle="--", linewidth=1, alpha=0.7, zorder=3
                ))

            площа_cyl_з_нахлестом += ov_w * ov_h

            ax_cyl.text(
                x_off + seg/2.0,
                y_off + (base_h/2.0 if base_h > 0 else 0.02),
                f"{seg:>5.2f}м",
                ha="center", va="center", fontsize=8, zorder=5
            )

            key_cyl = f"{round(seg, 2):>5.2f}м x {h_smuha:>3.2f}м"
            smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1

            x_off += seg

    ax_cyl.set_xlim(-circumference/2.0, circumference/2.0)
    ax_cyl.set_ylim(0, max_top_for_ylim)
    ax_cyl.set_xlabel("довжина поверхні (м)")
    ax_cyl.set_ylabel("висота (м)")
    ax_cyl.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # ============================================
    # 3) ШВИ + ПЛОЩІ + ЛИСТИ
    # ============================================
    площа_cyl = Wrem * L

    # Циліндр — шви
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
        n_vertical = (len(pat) - 1) + 2  # внутрішні + 2 крайові
        weld_v_cyl += n_vertical * h_eff

    weld_cyl = weld_h_cyl + weld_v_cyl

    # Днище — шви (1 днище)
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

    weld_bot_one = weld_h_bot_one + weld_side_bot_one
    weld_both_bottoms = 2.0 * weld_bot_one

    total_weld = weld_cyl + weld_both_bottoms

    реально_ремонт_площа = total_area_both_bottoms + площа_cyl

    площа_днищ_з_нахлестом_факт = 2.0 * sum(
        w * (h + (overlap/2.0 if idx < n_bot - 1 else 0.0))
        for idx, (w, h) in enumerate(zip(widths_bot, used_heights_bot))
    )
    фактична_площа_матеріалів = площа_днищ_з_нахлестом_факт + площа_cyl_з_нахлестом

    sheet_area = 6.0 * 1.5
    total_sheets = math.ceil(фактична_площа_матеріалів / sheet_area)

    left_lines = [
        f"Площа одного днища: {cum_area_bot:.3f} м²",
        f"Площа обох днищ:    {total_area_both_bottoms:.3f} м²",
        f"Площа циліндричної частини: {площа_cyl:.3f} м²",
        f"Реально ремонтована площа: {реально_ремонт_площа:.3f} м²",
        f"Фактична використана площа матеріалів (з нахлестом {overlap_cm} см): {фактична_площа_матеріалів:.3f} м²",
        f"Загальна довжина швів: {total_weld:.2f} м (днища: {weld_both_bottoms:.2f} м, циліндр: {weld_cyl:.2f} м)",
        f"Висота ремонтованої частини днища: {Hcrit:.3f} м",
        f"Довжина кола: {circumference:.3f} м",
        f"Кількість рядів на циліндрі: {full_rows} шт",
    ]

    items = sorted(smuhaDict.items())
    right_lines = ["Розмір смуги     К-сть"]
    for key, cnt in items:
        right_lines.append(f"{key:>12s}   {cnt:>3d} шт")
    right_lines.append("")
    right_lines.append(f"К-сть листів (6×1.5 м): {total_sheets} шт")
    right_lines.append(f"Використаний нахлест: {overlap_cm} см")

    # ============================================
    # ВИВІД
    # ============================================
    st.markdown("## Результати розрахунків")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("### Графіки")
        st.pyplot(fig_bot, clear_figure=False)
        st.pyplot(fig_cyl, clear_figure=False)

    with c2:
        st.markdown("### Підсумкові дані")
        st.markdown("#### Основні параметри")
        for line in left_lines:
            if line.startswith("Реально ремонтована площа"):
                st.markdown(f"<p style='color:blue;'>{line}</p>", unsafe_allow_html=True)
            elif line.startswith("Фактична використана площа матеріалів"):
                st.markdown(f"<p style='color:red;'>{line}</p>", unsafe_allow_html=True)
            elif line.startswith("Загальна довжина швів"):
                st.markdown(f"<p style='color:green;'>{line}</p>", unsafe_allow_html=True)
            else:
                st.text(line)

        st.markdown("#### Смуги матеріалу")
        for line in right_lines:
            st.text(line)

    st.markdown("---")
    st.markdown("### Збереження результатів")

    # PDF
    pdf_buffer = BytesIO()
    with PdfPages(pdf_buffer) as pdf:
        pdf.savefig(fig_bot)
        pdf.savefig(fig_cyl)

        fig_text = plt.figure(figsize=(8.27, 11.69))
        fig_text.clf()
        ax = fig_text.add_subplot(111)
        ax.axis("off")

        combined_text = "ПІДСУМКИ:\n\n" + "\n".join(left_lines) + "\n\n" + "\n".join(right_lines)
        ax.text(0.01, 0.99, combined_text, fontsize=10, va="top", ha="left", wrap=True)
        pdf.savefig(fig_text)

    st.download_button(
        label="Зберегти як PDF",
        data=pdf_buffer.getvalue(),
        file_name="резервуар_розрахунок.pdf",
        mime="application/pdf"
    )

    # PNG
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

    # Прибрати витоки пам’яті при багатьох перерахунках
    plt.close(fig_bot)
    plt.close(fig_cyl)
