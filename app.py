import streamlit as st
import matplotlib.pyplot as plt
import math
import numpy as np

st.set_page_config(page_title="Резервуар", layout="wide")

st.title("Розрахунок смуг для резервуара")

# --- Ввідні параметри ---
col1, col2, col3 = st.columns(3)
D = col1.number_input("Діаметр резервуара, м", min_value=0.1, value=12.0, step=0.1)
L = col2.number_input("Довжина резервуара, м", min_value=0.1, value=6.0, step=0.1)
Wrem = col3.number_input("Ширина ремонтуємої ділянки, м", min_value=0.1, value=3.0, step=0.5)

# --- Кнопка ---
if st.button("Розрахувати"):
    R = D / 2
    h_smuha = 0.5

    alpha = (Wrem / 2) / R
    Hcrit = R * (1 - math.cos(alpha))
    n_bot = max(1, math.ceil(Hcrit / h_smuha))

    widths_bot = []
    areas_bot = []

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")
    ax.add_patch(plt.Circle((0, 0), R, color='lightyellow', alpha=0.3))
    ax.axvline(0, color='red', linestyle='--')

    for j in range(1, n_bot + 1):
        y_bot = -R + (j - 1) * h_smuha
        y_top = y_bot + h_smuha
        y_ref = y_top if y_top <= 0 else y_bot if y_bot >= 0 else 0
        width = 2 * math.sqrt(max(0, R**2 - y_ref**2))
        widths_bot.append(width)
        areas_bot.append(width * h_smuha)

        x_left = -width / 2
        rect = plt.Rectangle((x_left, y_bot), width, h_smuha, color='skyblue', alpha=0.5)
        ax.add_patch(rect)

        # штриховка за межами кола
        dx = dy = 0.01
        x = x_left
        while x < -R or x > R:
            x += dx
        x = x_left
        while x < -R or x + dx > R:
            y = y_bot
            while y < y_top:
                if math.hypot(x + dx/2, y + dy/2) > R:
                    patch = plt.Rectangle((x, y), dx, dy, facecolor='none', edgecolor='red', hatch='///', linewidth=0.0)
                    ax.add_patch(patch)
                y += dy
            x += dx

        ax.text(0, y_bot + h_smuha / 2, f"S{j}\n{areas_bot[-1]:.2f}м²", ha='center')

    ax.set_xlim(-R * 1.2, R * 1.2)
    ax.set_ylim(-R * 1.2, R * 1.2)
    ax.set_title("Днище резервуара")
    st.pyplot(fig)

    st.success(f"Площа одного днища: {sum(areas_bot):.2f} м²")
