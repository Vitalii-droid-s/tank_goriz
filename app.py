import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages

st.set_page_config(layout="wide")
st.title("🔧 Розрахунок смуг для резервуара")

col1, col2, col3 = st.columns(3)
with col1:
    D = st.number_input("Діаметр резервуара (м)", min_value=0.1, step=0.1, value=5.0)
with col2:
    L = st.number_input("Довжина резервуара (м)", min_value=0.1, step=0.1, value=10.0)
with col3:
    Wrem = st.number_input("Ширина ремонту (м)", min_value=0.1, step=0.1, value=2.0)

h_smuha = 0.5
R = D / 2.0

if st.button("🔢 Розрахувати"):
    alpha = (Wrem / 2) / R
    Hcrit = R * (1 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))

    widths_bot = []
    areas_bot = []
    smuhaDict = {}

    fig1, ax1 = plt.subplots(figsize=(4, 4))
    ax1.set_aspect('equal')
    circle = plt.Circle((0, 0), R, edgecolor='black', facecolor='lightyellow', alpha=0.3)
    ax1.add_patch(circle)
    ax1.axvline(0, color='red', linestyle='--')

    for j in range(1, n_bot + 1):
        y_bot = -R + (j - 1) * h_smuha
        y_top = y_bot + h_smuha
        y_ref = y_top if y_top <= 0 else (y_bot if y_bot >= 0 else 0.0)

        width = 0.0 if abs(y_ref) >= R else 2 * math.sqrt(R ** 2 - y_ref ** 2)
        width_rnd = round(width, 2)
        key_bot = f"{width_rnd:.2f}м x {h_smuha:.2f}м"
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2

        widths_bot.append(width)
        areas_bot.append(width * h_smuha)

        x_left = -width / 2
        rect = plt.Rectangle((x_left, y_bot), width, h_smuha, edgecolor='black', facecolor='skyblue', alpha=0.5)
        ax1.add_patch(rect)

        dx = dy = 0.01
        x = x_left
        while x < -x_left:
            y = y_bot
            while y < y_top:
                xc, yc = x + dx / 2, y + dy / 2
                if math.hypot(xc, yc) > R:
                    patch = plt.Rectangle((x, y), dx, dy, facecolor='none', edgecolor='red',
                                          hatch='///', linewidth=0.0, alpha=0.5)
                    ax1.add_patch(patch)
                y += dy
            x += dx

        ax1.text(0, y_bot + h_smuha / 2, f"S{j}\n{areas_bot[-1]:.2f}м²", ha='center', va='center', fontsize=7)

    ax1.set_title("Днище")
    ax1.set_xlim(-R - 0.2, R + 0.2)
    ax1.set_ylim(-R - 0.2, R + 0.2)

    fig1, ax1 = plt.subplots(figsize=(4, 4))
    ax1.set_aspect('equal')
    circle = plt.Circle((0, 0), R, edgecolor='black', facecolor='lightyellow', alpha=0.3)
    ax1.add_patch(circle)
    ax1.axvline(0, color='red', linestyle='--')

    for j in range(1, n_bot + 1):
        y_bot = -R + (j - 1) * h_smuha
        y_top = y_bot + h_smuha
        y_ref = y_top if y_top <= 0 else (y_bot if y_bot >= 0 else 0.0)

        width = 0.0 if abs(y_ref) >= R else 2 * math.sqrt(R ** 2 - y_ref ** 2)
        width_rnd = round(width, 2)
        key_bot = f"{width_rnd:.2f}м x {h_smuha:.2f}м"
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2

        widths_bot.append(width)
        areas_bot.append(width * h_smuha)

        x_left = -width / 2
        rect = plt.Rectangle((x_left, y_bot), width, h_smuha, edgecolor='black', facecolor='skyblue', alpha=0.5)
        ax1.add_patch(rect)

        dx = dy = 0.01
        x = x_left
        while x < -x_left:
            y = y_bot
            while y < y_top:
                xc, yc = x + dx / 2, y + dy / 2
                if math.hypot(xc, yc) > R:
                    patch = plt.Rectangle((x, y), dx, dy, facecolor='none', edgecolor='red',
                                          hatch='///', linewidth=0.0, alpha=0.5)
                    ax1.add_patch(patch)
                y += dy
            x += dx

        ax1.text(0, y_bot + h_smuha / 2, f"S{j}\n{areas_bot[-1]:.2f}м²", ha='center', va='center', fontsize=7)

    ax1.set_title("Днище")
    ax1.set_xlim(-R - 0.2, R + 0.2)
    ax1.set_ylim(-R - 0.2, R + 0.2)

    fig1, ax1 = plt.subplots(figsize=(4, 4))
    ax1.set_aspect('equal')
    circle = plt.Circle((0, 0), R, edgecolor='black', facecolor='lightyellow', alpha=0.3)
    ax1.add_patch(circle)
    ax1.axvline(0, color='red', linestyle='--')

    for j in range(1, n_bot + 1):
        y_bot = -R + (j - 1) * h_smuha
        y_top = y_bot + h_smuha
        y_ref = y_top if y_top <= 0 else (y_bot if y_bot >= 0 else 0.0)

        width = 0.0 if abs(y_ref) >= R else 2 * math.sqrt(R ** 2 - y_ref ** 2)
        width_rnd = round(width, 2)
        key_bot = f"{width_rnd:.2f}м x {h_smuha:.2f}м"
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2

        widths_bot.append(width)
        areas_bot.append(width * h_smuha)

        x_left = -width / 2
        rect = plt.Rectangle((x_left, y_bot), width, h_smuha, edgecolor='black', facecolor='skyblue', alpha=0.5)
        ax1.add_patch(rect)

        dx = dy = 0.01
        x = x_left
        while x < -x_left:
            y = y_bot
            while y < y_top:
                xc, yc = x + dx / 2, y + dy / 2
                if math.hypot(xc, yc) > R:
                    patch = plt.Rectangle((x, y), dx, dy, facecolor='none', edgecolor='red',
                                          hatch='///', linewidth=0.0, alpha=0.5)
                    ax1.add_patch(patch)
                y += dy
            x += dx

        ax1.text(0, y_bot + h_smuha / 2, f"S{j}\n{areas_bot[-1]:.2f}м²", ha='center', va='center', fontsize=7)

    ax1.set_title("Днище")
    ax1.set_xlim(-R - 0.2, R + 0.2)
    ax1.set_ylim(-R - 0.2, R + 0.2)

    circumference = 2 * math.pi * R
    full_rows = math.ceil(L / h_smuha)

    fig2, ax2 = plt.subplots(figsize=(5, 3))
    ax2.set_aspect('auto')
    ax2.set_xlim(-circumference / 2, circumference / 2)
    ax2.set_ylim(0, full_rows * h_smuha)
    ax2.axvline(0, color='red', linestyle='--')
    ax2.axhline(L, color='red', linestyle='--')
    ax2.set_title("Циліндр")

    patterns = {
        1: [[1]],
        2: [[2]],
        3: [[3]],
        4: [[3, 1], [1, 3]],
        5: [[3, 2], [2, 3]],
        6: [[1, 3, 2], [2, 3, 1]]
    }.get(int(Wrem), [[Wrem]])

    for rowNum in range(full_rows):
        pattern = patterns[rowNum % len(patterns)]
        y_off = rowNum * h_smuha
        for seg in pattern:
            x_start = -Wrem / 2 + sum(pattern[:pattern.index(seg)])
            rect = plt.Rectangle((x_start, y_off), seg, h_smuha, edgecolor='black', facecolor='orange' if rowNum % 2 == 0 else 'lightgreen'

    circumference = 2 * math.pi * R
    full_rows = math.ceil(L / h_smuha))
    fig2, ax2 = plt.subplots(figsize=(5, 3))
    ax2.set_aspect('auto')

    # --- ФОН: корпус резервуара сірий, за його межами білий ---
    y_top = full_rows * h_smuha

    фон_біла_зона = plt.Rectangle(
        (-circumference / 2, 0),
        circumference,
        y_top,
        edgecolor='none',
        facecolor='white',
        zorder=0
    )
    ax2.add_patch(фон_біла_зона)

    фон_резервуара = plt.Rectangle(
        (-circumference / 2, 0),
        circumference,
        L,
        edgecolor='black',
        facecolor='#d0d0d0',
        linewidth=1,
        zorder=1
    )
    ax2.add_patch(фон_резервуара)

    # Вертикальні і горизонтальні лінії
    ax2.axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax2.axhline(L, color='red', linestyle='--', linewidth=1, alpha=0.7)

    
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


    # Викладка смуг
    x_start = -Wrem / 2
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
                ax2.add_patch(rect)
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
                    ax2.add_patch(rect_used)
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
                    ax2.add_patch(rect_extra)
                y_label = y_off + (visible_height / 2 if visible_height > 0 else 0)

            ax2.text(
                x_off + seg / 2,
                y_label,
                f"{seg:.2f}м",
                ha='center', va='center',
                fontsize=7,
                zorder=4
            )
            key_cyl = f"{seg:.2f}м x {h_smuha:.2f}м"
            smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1
            x_off += seg

    ax2.set_xlim(-circumference / 2, circumference / 2)
    ax2.set_ylim(0, y_top)
    ax2.set_xlabel('довжина поверхні (м)', fontsize=11)
    ax2.set_ylabel('висота (м)', fontsize=11)
    ax2.set_title("Розгорнута поверхня циліндра", fontsize=13, weight='bold')
    ax2.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)
    st.pyplot(fig1)
    st.pyplot(fig2)


    cum_area_bot = sum(areas_bot)
    площа_cyl = full_rows * h_smuha * Wrem
    total_area = площа_cyl + 2 * cum_area_bot
    total_sheets = math.ceil(total_area / (6 * 1.5))

    left = f"""
**Площа одного днища:** {cum_area_bot:.2f} м²  
**Площа обох днищ:** {2 * cum_area_bot:.2f} м²  
**Площа циліндра:** {площа_cyl:.2f} м²  
**Загальна площа:** {total_area:.2f} м²  
**Кількість листів (6×1.5 м):** {total_sheets} шт
"""

    right = "**Розміри смуг:**\n"
    for k, v in sorted(smuhaDict.items()):
        right += f"- {k}: {v} шт\n"

    col_l, col_r = st.columns(2)
    with col_l:
        st.markdown(left)
    with col_r:
        st.markdown(right)

    buffer = BytesIO()
    with PdfPages(buffer) as pdf:
        pdf.savefig(fig1, bbox_inches='tight')
        pdf.savefig(fig2, bbox_inches='tight')

        fig_text, ax = plt.subplots(figsize=(8.27, 11.69))
        ax.axis('off')
        ax.text(0.01, 0.99, left + "\n" + right, fontsize=10, va='top', ha='left')
        pdf.savefig(fig_text)

    st.download_button(
        label="📄 Завантажити PDF-файл",
        data=buffer.getvalue(),
        file_name="резервуар_розрахунок.pdf",
        mime="application/pdf"
    )
