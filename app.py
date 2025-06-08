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

    # Підсумки
    cum_area = sum(areas_bot)
    total_area = cum_area * 2
    st.subheader("Підсумки")
    st.text(f"Площа одного днища: {cum_area:.3f} м²")
    st.text(f"Площа обох днищ: {total_area:.3f} м²")
