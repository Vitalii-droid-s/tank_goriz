import streamlit as st
import matplotlib.pyplot as plt
import math
from io import BytesIO

st.set_page_config(layout="wide")
st.title("Розрахунок ремонту резервуара")

# Вхідні дані
D = st.number_input("Діаметр резервуара, м", value=2.5, step=0.1)
L = st.number_input("Довжина резервуара, м", value=4.3, step=0.1)
Wrem = st.number_input("Ширина ремонтуємої частини, м", value=1.0, step=0.1)
h_smuha = 0.5  # постійна висота смуги

R = D / 2
Hcrit = Wrem
n_bot = math.ceil(Hcrit / h_smuha)
n_cyl = math.ceil(L / 3)  # кількість 3-метрових смуг по довжині

# Графік днища
fig_bot, ax_bot = plt.subplots(figsize=(6, 6))
ax_bot.set_aspect('equal')
ax_bot.set_xlim(-R - 0.5, R + 0.5)
ax_bot.set_ylim(-R - 0.2, R + 0.2)

# Коло днища
circle = plt.Circle((0, 0), R, color='lightgrey', zorder=0)
ax_bot.add_patch(circle)

# Межа ремонту (червона лінія)
ax_bot.axhline(-R + Hcrit, color='red', linestyle='--', linewidth=1)
ax_bot.text(R * 0.6, -R + Hcrit + 0.03, f"H = {Hcrit:.3f} м", fontsize=10,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'), zorder=10)

# Побудова смуг на днищі
smuhaDict_bot = {}
cum_area_bot = 0.0
y_bot_global = -R
for j in range(n_bot):
    y_bot = y_bot_global + j * h_smuha
    y_top = y_bot + h_smuha
    y_ref = y_bot
    if abs(y_ref) >= R:
        width = 0.0
    else:
        width = 2 * math.sqrt(R * R - y_ref * y_ref)

    width_rnd = round(width, 2)
    key = f"{width_rnd:>5.2f}м x {h_smuha:>3.2f}м"
    smuhaDict_bot[key] = smuhaDict_bot.get(key, 0) + 2
    area = width * h_smuha
    cum_area_bot += area
    x_left = -width / 2

    rect = plt.Rectangle((x_left, y_bot), width, h_smuha,
                         edgecolor='black', facecolor='skyblue', alpha=0.5, zorder=2)
    ax_bot.add_patch(rect)

    if j == n_bot - 1:
        y_cut = -R + Hcrit
        real_h = y_cut - y_bot
        real_h = min(real_h, h_smuha)
        if real_h > 0:
            rect_real = plt.Rectangle((x_left, y_bot), width, real_h,
                                      edgecolor='black', facecolor='skyblue', alpha=0.5, zorder=3)
            ax_bot.add_patch(rect_real)
        extra_h = h_smuha - real_h
        if extra_h > 0:
            rect_extra = plt.Rectangle((x_left, y_bot + real_h), width, extra_h,
                                       edgecolor='red', facecolor='none', hatch='///', alpha=0.5, zorder=4)
            ax_bot.add_patch(rect_extra)

# Побудова циліндричної частини
fig_cyl, ax_cyl = plt.subplots(figsize=(8, 2))
ax_cyl.set_xlim(0, L)
ax_cyl.set_ylim(0, D)
ax_cyl.set_aspect('auto')
ax_cyl.set_title("Циліндрична частина")

smuhaDict_cyl = {}
cum_area_cyl = 0.0
for i in range(n_cyl):
    for j in range(int(D / h_smuha)):
        x0 = i * 3
        y0 = j * h_smuha
        rect = plt.Rectangle((x0, y0), 3, h_smuha,
                             edgecolor='black', facecolor='lightgreen', alpha=0.6)
        ax_cyl.add_patch(rect)
        key = f"3.00м x {h_smuha:>3.2f}м"
        smuhaDict_cyl[key] = smuhaDict_cyl.get(key, 0) + 1
        cum_area_cyl += 3 * h_smuha

# Вивід
col1, col2 = st.columns([1, 2])
with col1:
    st.markdown("### Висновки")
    st.text(f"Висота ремонтуємої частини резервуара: {Hcrit:.3f} м")
    st.text(f"Площа одного днища: {cum_area_bot:.3f} м²")
    for key, val in smuhaDict_bot.items():
        st.text(f"{key}: {val} шт")
    st.markdown("---")
    st.text(f"Площа циліндричної частини: {cum_area_cyl:.3f} м²")
    for key, val in smuhaDict_cyl.items():
        st.text(f"{key}: {val} шт")

    buf_png = BytesIO()
    fig_bot.savefig(buf_png, format="png")
    st.download_button("⬇️ Завантажити PNG (днище)", data=buf_png.getvalue(), file_name="dnysche.png", mime="image/png")

    buf_pdf = BytesIO()
    fig_bot.savefig(buf_pdf, format="pdf")
    st.download_button("⬇️ Завантажити PDF (днище)", data=buf_pdf.getvalue(), file_name="dnysche.pdf", mime="application/pdf")

with col2:
    st.pyplot(fig_bot)
    st.pyplot(fig_cyl)
