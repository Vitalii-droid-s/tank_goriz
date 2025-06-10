import streamlit as st
import matplotlib.pyplot as plt
import math
from io import BytesIO
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

st.set_page_config(layout="wide")
st.title("Розрахунок смуг для резервуара")

# Вхідні дані
col1, col2, col3 = st.columns(3)
with col1:
    D = st.number_input("Діаметр резервуара, м", value=2.5, min_value=0.1, step=0.1, format="%.1f")
with col2:
    L = st.number_input("Довжина резервуара, м", value=4.3, min_value=0.1, step=0.1, format="%.1f")
with col3:
    Wrem = st.number_input("Ширина ремонтуємої ділянки, м", value=1.0, min_value=0.1, step=0.1, format="%.1f")

if st.button("Розрахувати"):
    R = D / 2.0
    h_smuha = 0.5
    alpha = (Wrem / 2) / R
    Hcrit = R * (1 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))
    
    # ============================================
    # 1) Днище
    # ============================================
    widths_bot = []
    areas_bot = []
    smuhaDict = {}
    
    fig_bot, ax_bot = plt.subplots(figsize=(6, 6))
    ax_bot.set_aspect('equal')
    
    # Малюємо контур кола фоном
    circle = plt.Circle(
        (0, 0), R,
        edgecolor='black',
        facecolor='lightyellow',
        alpha=0.3,
        zorder=0
    )
    ax_bot.add_patch(circle)
    
    # Додаємо вертикальну червону пунктирну лінію для x=0
    ax_bot.axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.7)
    # Підпис Hcrit на графіку
    ax_bot.text(
        R * 0.6,                  # x-координата праворуч
        -R + Hcrit + 0.03,        # y = верх зеленої смуги
        f"H = {Hcrit:.3f} м",
        fontsize=10,
        color='black',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'),
        zorder=10
    )
    
    # Початок від центра кола (y = 0)
    y_bot_global = -R  # починаємо з нижньої точки кола
    
    for j in range(n_bot):
        y_bot = y_bot_global + j * h_smuha
        y_top = y_bot + h_smuha
        
        # y_ref — для обчислення ширини
        if y_top <= 0:
            y_ref = y_top
        elif y_bot >= 0:
            y_ref = y_bot
        else:
            y_ref = 0.0

        if abs(y_ref) >= R:
            width = 0.0
        else:
            width = 2 * math.sqrt(R * R - y_ref * y_ref)

        width_rnd = round(width, 2)
        key_bot = f"{width_rnd:>5.2f}м x {h_smuha:>3.2f}м"
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2

        widths_bot.append(width)
        areas_bot.append(width * h_smuha)

        x_left = -width / 2
        x_right = width / 2

        # 1. Стандартна смуга
        rect = plt.Rectangle(
            (x_left, y_bot),
            width,
            h_smuha,
            edgecolor='black',
            facecolor='skyblue',
            alpha=0.5,
            zorder=2
        )
        ax_bot.add_patch(rect)

        # 2. Реальна смуга — лише в останній (верхній) смузі
        if j == n_bot - 1:
            y_cut = -R + Hcrit  # межа, до якої потрібно
            real_h = y_cut - y_bot
            real_h = min(real_h, h_smuha)

            if real_h > 0:
                rect_real = plt.Rectangle(
                    (x_left, y_bot),
                    width,
                    real_h,
                    edgecolor='black',
                    facecolor='skyblue',
                    alpha=0.5,
                    zorder=3
                )
                ax_bot.add_patch(rect_real)

            extra_h = h_smuha - real_h
            if extra_h > 0:
                rect_extra = plt.Rectangle(
                    (x_left, y_bot + real_h),
                    width,
                    extra_h,
                    edgecolor='red',
                    facecolor='none',
                    hatch='///',
                    alpha=0.5,
                    zorder=4
                )
                ax_bot.add_patch(rect_extra)

        # 4. Штрихування зон за межами кола
        dx = 0.01
        dy = 0.01
        x = x_left
        while x < x_right:
            y = y_bot
            while y < y_top:
                xc = x + dx / 2
                yc = y + dy / 2
                dist = math.hypot(xc, yc)
                if dist > R:
                    patch = plt.Rectangle(
                        (x, y),
                        dx,
                        dy,
                        facecolor='none',
                        edgecolor='red',
                        hatch='///',
                        linewidth=0.0,
                        alpha=0.5,
                        zorder=5
                    )
                    ax_bot.add_patch(patch)
                y += dy
            x += dx

        # 5. Підпис площі
        ax_bot.text(
            0,
            y_bot + h_smuha / 2,
            f"S{j + 1}\n{areas_bot[j]:.3f}м²",
            ha='center', va='center',
            fontsize=8,
            zorder=6
        )

    # Площа одного днища
    cum_area_bot = sum(areas_bot)
    # Площа обох днищ:
    total_area_both_bottoms = 2 * cum_area_bot

    # Налаштування осей для днища
    margin = R * 0.1
    ax_bot.set_xlim(-R - margin, R + margin)
    ax_bot.set_ylim(-R - margin, R + margin)
    ax_bot.set_xlabel('x (м)', fontsize=11, fontfamily='Arial')
    ax_bot.set_ylabel('y (м)', fontsize=11, fontfamily='Arial')
    ax_bot.set_title("Днище резервуара", fontsize=14, fontweight='bold', fontfamily='Arial')
    ax_bot.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)

    # ============================================
    # 2) Циліндр
    # ============================================
    circumference = 2 * math.pi * R
    full_rows = math.ceil(L / h_smuha)
    total_height = full_rows * h_smuha

    fig_cyl, ax_cyl = plt.subplots(figsize=(8, 6))
    ax_cyl.set_aspect('auto')

    # Білий фон поза резервуаром
    фон_біла_зона = plt.Rectangle(
        (-circumference / 2, 0),
        circumference,
        total_height,
        edgecolor='none',
        facecolor='white',
        zorder=0
    )
    ax_cyl.add_patch(фон_біла_зона)

    # Сіра заливка корпусу резервуара
    фон_резервуара = plt.Rectangle(
        (-circumference / 2, 0),
        circumference,
        L,
        edgecolor='black',
        facecolor='#d0d0d0',
        linewidth=1,
        zorder=1
    )
    ax_cyl.add_patch(фон_резервуара)

    # Додаємо вертикальну червону пунктирну лінію для x=0
    ax_cyl.axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.7)
    # Додаємо горизонтальну червону пунктирну лінію для y=L
    ax_cyl.axhline(L, color='red', linestyle='--', linewidth=1, alpha=0.7)

    # Ремонтна зона Wrem по центру
    x_start = -Wrem / 2

    if abs(Wrem - 1) < 1e-6:
        patterns = [[1]]
    elif abs(Wrem - 2) < 1e-6:
        patterns = [[2]]
    elif abs(Wrem - 3) < 1e-6:
        patterns = [[3]]
    elif abs(Wrem - 4) < 1e-6:
        patterns = [[3, 1], [1, 3]]
    elif abs(Wrem - 5) < 1e-6:
        patterns = [[3, 2], [2, 3]]
    elif abs(Wrem - 6) < 1e-6:
        patterns = [[1, 3, 2], [2, 3, 1]]
    else:
        patterns = [[Wrem]]

    for rowNum in range(full_rows):
        pattern = patterns[rowNum % len(patterns)]
        y_off = rowNum * h_smuha
        x_off = x_start

        visible_height = min(h_smuha, max(0.0, L - y_off))
        is_partial = (visible_height < h_smuha)

        for seg in pattern:
            if not is_partial:
                rect = plt.Rectangle(
                    (x_off, y_off),
                    seg,
                    h_smuha,
                    edgecolor='black',
                    facecolor='orange' if (rowNum % 2 == 0) else 'lightgreen',
                    alpha=0.7,
                    zorder=2
                )
                ax_cyl.add_patch(rect)
                y_label = y_off + h_smuha / 2
            else:
                if visible_height > 0:
                    rect_used = plt.Rectangle(
                        (x_off, y_off),
                        seg,
                        visible_height,
                        edgecolor='black',
                        facecolor='orange' if (rowNum % 2 == 0) else 'lightgreen',
                        alpha=0.7,
                        zorder=2
                    )
                    ax_cyl.add_patch(rect_used)

                extra_h = h_smuha - visible_height
                if extra_h > 0:
                    rect_extra = plt.Rectangle(
                        (x_off, y_off + visible_height),
                        seg,
                        extra_h,
                        edgecolor='red',
                        facecolor='none',
                        hatch='///',
                        alpha=0.5,
                        zorder=3
                    )
                    ax_cyl.add_patch(rect_extra)

                y_label = y_off + (visible_height / 2 if visible_height > 0 else 0)

            ax_cyl.text(
                x_off + seg / 2,
                y_label,
                f"{seg:>5.2f}м",
                ha='center', va='center',
                fontsize=8,
                zorder=4
            )
            key_cyl = f"{round(seg, 2):>5.2f}м x {h_smuha:>3.2f}м"
            smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1
            x_off += seg

    ax_cyl.set_xlim(-circumference / 2, circumference / 2)
    ax_cyl.set_ylim(0, total_height)
    ax_cyl.set_xlabel('довжина поверхні (м)', fontsize=11, fontfamily='Arial')
    ax_cyl.set_ylabel('висота (м)', fontsize=11, fontfamily='Arial')
    ax_cyl.set_title("Розгорнута поверхня циліндра", fontsize=14, fontweight='bold', fontfamily='Arial')
    ax_cyl.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)

    # ============================================
    # 3) Підрахунок підсумкових величин
    # ============================================
    # Площа циліндра враховує повні смуги:
    площа_cyl = full_rows * h_smuha * Wrem

    horz_bot = 2 * sum(widths_bot)
    vert_bot = n_bot * 2 * h_smuha
    weld_bot = horz_bot + vert_bot

    horz_cyl = full_rows * 2 * Wrem
    vert_cyl = 0.0
    for idx in range(full_rows):
        pattern = patterns[idx % len(patterns)]
        vert_cyl += 2
        vert_cyl += (len(pattern) - 1) * h_smuha
    weld_cyl = horz_cyl + vert_cyl

    total_weld = weld_bot + weld_cyl

    # Площа обох днищ + площа циліндра (за повними смугами)
    заг_площа = total_area_both_bottoms + площа_cyl

    висота_ділянки = D
    довжина_окружності = circumference
    sheet_area = 6 * 1.5
    total_sheets = math.ceil(заг_площа / sheet_area)

    # — Ліва частина: основні підсумки —
    left_lines = [
        f"Площа одного днища: {cum_area_bot:.3f} м²",
        f"Площа обох днищ:    {total_area_both_bottoms:.3f} м²",
        f"Площа циліндричної частини:     {площа_cyl:.3f} м²",
        f"Загальна площа:     {заг_площа:.3f} м²",
        f"Заг. довжина швів:  {total_weld:.2f} м",
        f"Висота ремонтуємої частини резервуара: {Hcrit:.3f} м", 
        f"Довжина кола: {довжина_окружності:.3f} м",
    ]

    # — Права частина: смуги + кількість листів —
    items = sorted(smuhaDict.items())
    right_lines = ["Розмір смуги     К-сть"]
    for key, cnt in items:
        right_lines.append(f"{key:>12s}   {cnt:>3d} шт")
    right_lines.append("")  # порожній рядок
    right_lines.append(f"К-сть листів (6×1.5м): {total_sheets} шт")

    # Відображення результатів
    st.markdown("## Результати розрахунків")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("### Графіки")
        st.pyplot(fig_bot)
        st.pyplot(fig_cyl)
    
    with col2:
        st.markdown("### Підсумкові дані")
        st.markdown("#### Основні параметри")
        for line in left_lines:
            st.text(line)
        
        st.markdown("#### Смуги матеріалу")
        for line in right_lines:
            st.text(line)
    
    # Кнопки для збереження
    st.markdown("---")
    st.markdown("### Збереження результатів")
    
    # Збереження графіків у PDF
    pdf_buffer = BytesIO()
    with PdfPages(pdf_buffer) as pdf:
        pdf.savefig(fig_bot)
        pdf.savefig(fig_cyl)
        
        # Додаємо текстову сторінку
        fig_text = plt.figure(figsize=(8.27, 11.69))
        fig_text.clf()
        ax = fig_text.add_subplot(111)
        ax.axis('off')
        
        combined_text = "ПІДСУМКИ:\n\n" + "\n".join(left_lines) + "\n\n" + "\n".join(right_lines)
        ax.text(0.01, 0.99, combined_text, fontsize=10, va='top', ha='left', wrap=True)
        pdf.savefig(fig_text)
    
    st.download_button(
        label="Зберегти як PDF",
        data=pdf_buffer.getvalue(),
        file_name="резервуар_розрахунок.pdf",
        mime="application/pdf"
    )
    
    # Збереження графіків окремо
    col1, col2 = st.columns(2)
    with col1:
        png_buffer_bot = BytesIO()
        fig_bot.savefig(png_buffer_bot, format="png")
        st.download_button(
            label="Зберегти днище (PNG)",
            data=png_buffer_bot.getvalue(),
            file_name="днище.png",
            mime="image/png"
        )
    with col2:
        png_buffer_cyl = BytesIO()
        fig_cyl.savefig(png_buffer_cyl, format="png")
        st.download_button(
            label="Зберегти циліндр (PNG)",
            data=png_buffer_cyl.getvalue(),
            file_name="циліндр.png",
            mime="image/png"
        )