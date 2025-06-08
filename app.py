import streamlit as st
import matplotlib.pyplot as plt
import math
import numpy as np

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

    fig1, ax1 = plt.subplots(figsize=(6, 6))
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
    # Побудова таблиці смуг
    circumference = 2 * math.pi * R
    full_rows = math.ceil(L / h_smuha)

    # Словник смуг
    smuhaDict = {}
    for j in range(1, n_bot + 1):
        width_rnd = round(widths_bot[j - 1], 2)
        key_bot = f"{width_rnd:>5.2f}м x {h_smuha:>3.2f}м"
        smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2

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
        for seg in pattern:
            key_cyl = f"{round(seg, 2):>5.2f}м x {h_smuha:>3.2f}м"
            smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1

    площа_cyl = full_rows * h_smuha * Wrem
    total_area_both_bottoms = 2 * sum(areas_bot)
    заг_площа = total_area_both_bottoms + площа_cyl
    total_sheets = math.ceil(заг_площа / (6 * 1.5))

    st.markdown("### 📋 Підсумки та кількість смуг")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"- **Площа циліндра:** `{площа_cyl:.3f} м²`")
        st.markdown(f"- **Загальна площа:** `{заг_площа:.3f} м²`")
        st.markdown(f"- **К-сть листів (6×1.5 м):** `{total_sheets}`")

    with col2:
        st.markdown("**Кількість смуг:**")
        for k, v in sorted(smuhaDict.items()):
            st.markdown(f"- `{k}` — `{v}` шт")
    

    # Підсумки
    cum_area = sum(areas_bot)
    total_area = cum_area * 2
    st.subheader("Підсумки")
    st.text(f"Площа одного днища: {cum_area:.3f} м²")
    st.text(f"Площа обох днищ: {total_area:.3f} м²")
