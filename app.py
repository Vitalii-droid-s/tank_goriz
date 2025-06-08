import streamlit as st
import matplotlib.pyplot as plt
import math
import numpy as np
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages

st.set_page_config(layout="wide")
st.title("Розрахунок смуг для резервуара")

# Ввід параметрів
st.sidebar.header("Вхідні параметри")
D = st.sidebar.number_input("Діаметр резервуара (м)", min_value=0.1, value=5.0, step=0.1)
L = st.sidebar.number_input("Довжина резервуара (м)", min_value=0.1, value=10.0, step=0.1)
Wrem = st.sidebar.number_input("Ширина ремонтної ділянки (м)", min_value=0.1, value=3.0, step=0.1)

if st.sidebar.button("Розрахувати"):
    R = D / 2
    h_smuha = 0.5
    alpha = (Wrem / 2) / R
    Hcrit = R * (1 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))

    widths_bot = []
    areas_bot = []

    fig1, ax1 = plt.subplots(figsize=(4, 4))
    ax1.set_aspect('equal')
    circle = plt.Circle((0, 0), R, edgecolor='black', facecolor='lightyellow', alpha=0.3)
    ax1.add_patch(circle)
    ax1.axvline(0, color='red', linestyle='--')

    for j in range(1, n_bot + 1):
        y_bot = -R + (j - 1) * h_smuha
        y_top = y_bot + h_smuha
        y_ref = y_top if y_top <= 0 else y_bot if y_bot >= 0 else 0
        width = 0 if abs(y_ref) >= R else 2 * math.sqrt(R**2 - y_ref**2)
        widths_bot.append(width)
        areas_bot.append(width * h_smuha)
        rect = plt.Rectangle((-width/2, y_bot), width, h_smuha, edgecolor='black', facecolor='skyblue', alpha=0.5)
        ax1.add_patch(rect)

        dx = dy = 0.01
        x = -width / 2
        while x < width / 2:
            y = y_bot
            while y < y_top:
                if math.hypot(x + dx/2, y + dy/2) > R:
                    patch = plt.Rectangle((x, y), dx, dy, facecolor='none', edgecolor='red', hatch='///', linewidth=0.0, alpha=0.5)
                    ax1.add_patch(patch)
                y += dy
            x += dx

        ax1.text(0, y_bot + h_smuha / 2, f"S{j}\n{areas_bot[j - 1]:.3f}м²", ha='center', va='center', fontsize=8)

    margin = R * 0.1
    ax1.set_xlim(-R - margin, R + margin)
    ax1.set_ylim(-R - margin, R + margin)
    ax1.set_title("Днище резервуара")
    st.pyplot(fig1)

    # === Побудова циліндра ===
    circumference = 2 * math.pi * R
    full_rows = math.ceil(L / h_smuha)
    total_height = full_rows * h_smuha

    fig2, ax2 = plt.subplots(figsize=(5, 3))
    ax2.set_aspect('auto')
    ax2.add_patch(plt.Rectangle((-circumference/2, 0), circumference, L, edgecolor='black', facecolor='#d0d0d0', zorder=1))
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

    for rowNum in range(full_rows):
        pattern = patterns[rowNum % len(patterns)]
        y_off = rowNum * h_smuha
        x_off = -Wrem / 2
        visible_height = min(h_smuha, max(0.0, L - y_off))
        is_partial = (visible_height < h_smuha)

        for seg in pattern:
            if not is_partial:
                ax2.add_patch(plt.Rectangle((x_off, y_off), seg, h_smuha,
                                            edgecolor='black', facecolor='orange' if (rowNum % 2 == 0) else 'lightgreen',
                                            alpha=0.7, zorder=2))
                y_label = y_off + h_smuha / 2
            else:
                if visible_height > 0:
                    ax2.add_patch(plt.Rectangle((x_off, y_off), seg, visible_height,
                                                edgecolor='black', facecolor='orange' if (rowNum % 2 == 0) else 'lightgreen',
                                                alpha=0.7, zorder=2))
                extra_h = h_smuha - visible_height
                if extra_h > 0:
                    ax2.add_patch(plt.Rectangle((x_off, y_off + visible_height), seg, extra_h,
                                                edgecolor='red', facecolor='none', hatch='///',
                                                alpha=0.5, zorder=3))
                y_label = y_off + (visible_height / 2 if visible_height > 0 else 0)

            ax2.text(x_off + seg / 2, y_label, f"{seg:.2f}м", ha='center', va='center', fontsize=8, zorder=4)
            x_off += seg

    ax2.set_xlim(-circumference / 2, circumference / 2)
    ax2.set_ylim(0, total_height)
    ax2.set_title("Розгорнута поверхня циліндра")
    st.pyplot(fig2)


    # Зберігаємо у session_state
    st.session_state['fig1'] = fig1
    st.session_state['fig2'] = fig2
    st.session_state['summary_text'] = f"""ПІДСУМКИ:

Площа одного днища: {cum_area_bot:.3f} м²
Площа обох днищ:    {2 * cum_area_bot:.3f} м²
Площа циліндра:     {площа_cyl:.3f} м²
Загальна площа:     {2 * cum_area_bot + площа_cyl:.3f} м²
Заг. довжина швів:  {total_weld:.2f} м
Висота ділянки:     {D:.2f} м
Довжина окружності: {2 * math.pi * (D/2):.3f} м"""



    # PDF export
    
# Збереження в PDF
st.subheader("⬇️ Збереження результатів у PDF")
if st.button("Завантажити PDF"):
    if 'fig1' in st.session_state and 'fig2' in st.session_state:
        buffer = BytesIO()
        with PdfPages(buffer) as pdf:
            pdf.savefig(st.session_state['fig1'], bbox_inches='tight')
            pdf.savefig(st.session_state['fig2'], bbox_inches='tight')

            # Сторінка з підсумками
            fig_text = plt.figure(figsize=(8.27, 11.69))
            fig_text.clf()
            ax_text = fig_text.add_subplot(111)
            ax_text.axis('off')
            ax_text.text(0.01, 0.99, st.session_state.get('summary_text', ''),
                         fontsize=10, va='top', ha='left', wrap=True)
            pdf.savefig(fig_text)

                    file_name="резервуар_розрахунок.pdf",
            mime="application/pdf"
        )
    else:
        st.warning("⚠️ Спочатку натисніть 'Розрахувати'")

                           file_name="резервуар_розрахунок.pdf", mime="application/pdf")

                    file_name="резервуар_розрахунок.pdf",
            mime="application/pdf"
        )


        st.download_button(
            label="📄 Завантажити PDF-файл",
            data=buffer.getvalue(),
            file_name="резервуар_розрахунок.pdf",
            mime="application/pdf"
        )
