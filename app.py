import streamlit as st
import matplotlib.pyplot as plt
import math
from io import BytesIO
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

VERSION = "1.2"

st.set_page_config(layout="wide")
st.title(f"Розрахунок смуг для резервуара — v{VERSION}")

# Вхідні дані
col1, col2, col3, col4 = st.columns(4)
with col1:
    D = st.number_input("Діаметр резервуара, м", value=2.5, min_value=0.1, step=0.1, format="%.1f")
with col2:
    L = st.number_input("Довжина резервуара, м", value=4.3, min_value=0.1, step=0.1, format="%.1f")
with col3:
    Wrem = st.number_input("Ширина ремонтуємої ділянки, м", value=1.0, min_value=0.1, step=0.1, format="%.1f")
with col4:
    overlap_cm = st.number_input("Нахлест, см", value=5.0, min_value=0.0, step=0.5, format="%.1f")

if st.button("Розрахувати"):
    if D <= 0 or L <= 0 or Wrem <= 0 or overlap_cm < 0:
        st.error("Введіть, будь ласка, додатні числа для D, L та Wrem і невід'ємне значення для нахлеста.")
    else:
        overlap = overlap_cm / 100.0  # Переводимо см в метри
        R = D / 2.0
        h_smuha = 0.5
        alpha = (Wrem / 2) / R
        Hcrit = R * (1 - math.cos(alpha))
        n_bot = max(1, math.ceil(Hcrit / h_smuha))
        
        # ============================================
        # 1) Днище
        # ============================================
        widths_bot = []          # фактичні ширини смуг на обраному y
        used_heights_bot = []    # фактичні висоти смуг (останню «обрізаємо»)
        areas_bot = []           # площі смуг з урахуванням «обрізки»
        smuhaDict = {}
        
        fig_bot, ax_bot = plt.subplots(figsize=(6, 6))
        ax_bot.set_aspect('equal')
        ax_bot.set_title("Днище резервуара (синій - основна смуга, червоний - нахлест)", 
                        fontsize=14, fontweight='bold', fontfamily='Arial')
        
        # Контур кола фоном
        circle = plt.Circle((0, 0), R, edgecolor='black', facecolor='lightyellow', alpha=0.3, zorder=0)
        ax_bot.add_patch(circle)

        # Вертикальна лінія x=0 і підпис Hcrit
        ax_bot.axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.7)
        ax_bot.text(
            R * 0.60, -R + Hcrit + 0.03,
            f"H = {Hcrit:.3f} м",
            fontsize=10, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'),
            zorder=10
        )

        # Починаємо з нижньої точки кола
        y_bot_global = -R

        for j in range(n_bot):
            y_bot = y_bot_global + j * h_smuha
            y_top = y_bot + h_smuha

            # y_ref для обчислення ширини
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
            smuhaDict[key_bot] = smuhaDict.get(key_bot, 0) + 2  # по дві симетричні смуги

            x_left = -width / 2
            x_right = width / 2

            # --- ФАКТИЧНА ВИСОТА СМУГИ ---
            real_h = h_smuha
            if j == n_bot - 1:  # верхня смуга може бути «обрізаною»
                y_cut = -R + Hcrit
                real_h = max(0.0, min(y_cut - y_bot, h_smuha))

            widths_bot.append(width)
            used_heights_bot.append(real_h)

            # точна площа смуги між y_bot і y_top_clip у колі радіуса R
            def _F(y):  # первісна для 2*sqrt(R^2 - y^2)
                y = max(-R, min(R, y))
                return y * math.sqrt(max(0.0, R*R - y*y)) + R*R * math.asin(max(-1.0, min(1.0, y / R)))

            y_top_clip = min(y_bot + real_h, -R + Hcrit)  # обрізати верх до потрібної висоти
            area_strip_exact = _F(y_top_clip) - _F(y_bot)
            areas_bot.append(area_strip_exact)

            # 1) Стандартна зона смуги (фон)
            base_rect = plt.Rectangle(
                (x_left, y_bot), width, h_smuha,
                edgecolor='black', facecolor='skyblue', alpha=0.5, zorder=2
            )
            ax_bot.add_patch(base_rect)

            # 2) Візуалізація нахлеста (ЛИШЕ по вертикалі між рядами)
            ov_h = h_smuha + (overlap/2 if j < n_bot - 1 else 0.0)
            overlap_rect = plt.Rectangle(
                (x_left, y_bot),         # без бокового зсуву
                width,                   # без +overlap по ширині
                ov_h,                    # тільки вертикальний нахлест
                edgecolor='red', facecolor='none',
                linestyle='--', linewidth=1, alpha=0.7, zorder=3
            )
            ax_bot.add_patch(overlap_rect)

            # 3) Реально використана частина (поверх «фону»)
            if real_h > 0:
                rect_real = plt.Rectangle(
                    (x_left, y_bot), width, real_h,
                    edgecolor='black', facecolor='skyblue', alpha=0.5, zorder=4
                )
                ax_bot.add_patch(rect_real)

            # 4) Невикористана частина верхньої смуги — штрихуванням
            extra_h = h_smuha - real_h
            if extra_h > 0:
                rect_extra = plt.Rectangle(
                    (x_left, y_bot + real_h), width, extra_h,
                    edgecolor='red', facecolor='none', hatch='///', alpha=0.5, zorder=5
                )
                ax_bot.add_patch(rect_extra)

            # 5) Штрихування ділянок поза колом
            dx = 0.01
            dy = 0.01
            x = x_left
            while x < x_right:
                y = y_bot
                while y < y_top:
                    xc = x + dx / 2
                    yc = y + dy / 2
                    if math.hypot(xc, yc) > R:
                        ax_bot.add_patch(plt.Rectangle(
                            (x, y), dx, dy,
                            facecolor='none', edgecolor='red',
                            hatch='///', linewidth=0.0, alpha=0.5, zorder=6
                        ))
                    y += dy
                x += dx

            # 6) Підпис площі смуги (за фактичною висотою real_h)
            ax_bot.text(
                0, y_bot + h_smuha / 2,
                f"S{j + 1}\n{areas_bot[-1]:.3f}м²",
                ha='center', va='center', fontsize=8, zorder=7
            )

        # Площа одного днища (за фактичною висотою верхньої смуги)
        cum_area_bot = sum(areas_bot)
        # Площа обох днищ
        total_area_both_bottoms = 2 * cum_area_bot

        # Оформлення осей
        margin = R * 0.1
        ax_bot.set_xlim(-R - margin, R + margin)
        ax_bot.set_ylim(-R - margin, R + margin)
        ax_bot.set_xlabel('x (м)', fontsize=11, fontfamily='Arial')
        ax_bot.set_ylabel('y (м)', fontsize=11, fontfamily='Arial')
        ax_bot.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)

        # ============================================
        # 2) Циліндр (останній ряд повністю; верх – штрих; без верхнього нахльосту)
        # ============================================
        circumference = 2 * math.pi * R

        # К-сть рядів з урахуванням нахльосту
        full_rows = 1
        covered_height = h_smuha
        epsilon = 1e-6
        while covered_height + epsilon < L:
            full_rows += 1
            covered_height += h_smuha - overlap  # кожний новий ряд перекриває попередній на overlap

        fig_cyl, ax_cyl = plt.subplots(figsize=(8, 6))
        ax_cyl.set_aspect('auto')
        ax_cyl.set_title("Розгорнута поверхня циліндра (кольори - чергування смуг, червоний - нахлест)", 
                        fontsize=14, fontweight='bold', fontfamily='Arial')

        # Сіра зона корпуса (рівно L)
        ax_cyl.add_patch(plt.Rectangle(
            (-circumference / 2, 0), circumference, L,
            edgecolor='black', facecolor='#d0d0d0', linewidth=1, zorder=1
        ))

        # Осьові лінії
        ax_cyl.axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.7)
        ax_cyl.axhline(L, color='red', linestyle='--', linewidth=1, alpha=0.7)

        # Ремонтна зона по центру
        x_start = -Wrem / 2

        # Візерунки
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

        # Акумулятор площі циліндра з нахльостом
        площа_cyl_з_нахлестом = 0.0

        # Малювання рядів
        max_top_for_ylim = L
        for rowNum in range(full_rows):
            y_off = rowNum * (h_smuha - overlap)  # зсув рядів з урахуванням нахльосту
            if y_off >= L + overlap/2:
                break

            pattern = patterns[rowNum % len(patterns)]
            x_off = x_start

            # Видима частина ряду в межах L
            visible_height = max(0.0, min(h_smuha, L - y_off))
            # Нижній додаток від нахльосту (не нижче 0)
            extra_down = min(overlap/2, max(0.0, y_off))
            # Чи це останній ряд?
            is_last_row = (rowNum == full_rows - 1)

            # Верхня межа прямокутника нахльосту:
            # - для звичайних рядів: + overlap/2 зверху;
            # - для останнього ряду: БЕЗ верхнього нахльосту.
            top_ov_end = y_off + visible_height + (0.0 if is_last_row else overlap/2)
            ov_y_base = max(0.0, y_off - extra_down)
            ov_h = max(0.0, top_ov_end - ov_y_base)

            # Для ліміту осі Y показуємо повну висоту останнього ряду
            if is_last_row:
                max_top_for_ylim = max(max_top_for_ylim, y_off + h_smuha)
            else:
                max_top_for_ylim = max(max_top_for_ylim, top_ov_end)

            for seg_idx, seg in enumerate(pattern):
                # 1) БАЗОВА смуга:
                #    - для останнього ряду малюємо ПОВНУ висоту h_smuha,
                #    - для інших — лише видиму частину в межах L.
                base_h = h_smuha if is_last_row else visible_height
                ax_cyl.add_patch(plt.Rectangle(
                    (x_off, y_off), seg, base_h,
                    edgecolor='black',
                    facecolor=('orange' if (rowNum % 2 == 0) else 'lightgreen'),
                    alpha=0.7, zorder=2
                ))

                # 1.1) Якщо останній ряд виходить вище L — заштрихувати надлишок
                if is_last_row and (y_off + h_smuha > L):
                    ax_cyl.add_patch(plt.Rectangle(
                        (x_off, L), seg, (y_off + h_smuha - L),
                        edgecolor='red', facecolor='none', hatch='///', alpha=0.5, zorder=4
                    ))

                # 2) Нахльост (червоний пунктир)
                if Wrem in [4, 5, 6] and 0 < seg_idx < len(pattern) - 1:
                    ov_x = x_off - overlap/2
                    ov_w = seg + overlap
                else:
                    ov_x = x_off
                    ov_w = seg

                # НЕ малюємо верхній нахльост для останнього ряду (ov_h уже без +overlap/2)
                if ov_h > 0:
                    ax_cyl.add_patch(plt.Rectangle(
                        (ov_x, ov_y_base), ov_w, ov_h,
                        edgecolor='red', facecolor='none',
                        linestyle='--', linewidth=1, alpha=0.7, zorder=3
                    ))

                # 3) НАКОПИЧЕННЯ площі матеріалу з нахльостом
                площа_cyl_з_нахлестом += ov_w * ov_h

                # 4) Підпис
                ax_cyl.text(
                    x_off + seg / 2,
                    y_off + (base_h/2 if base_h > 0 else 0.02),
                    f"{seg:>5.2f}м", ha='center', va='center', fontsize=8, zorder=5
                )

                # 5) Облік смуг
                key_cyl = f"{round(seg, 2):>5.2f}м x {h_smuha:>3.2f}м"
                smuhaDict[key_cyl] = smuhaDict.get(key_cyl, 0) + 1

                x_off += seg

        # Межі осей: показати повний останній ряд
        ax_cyl.set_xlim(-circumference / 2, circumference / 2)
        ax_cyl.set_ylim(0, max_top_for_ylim)
        ax_cyl.set_xlabel('довжина поверхні (м)', fontsize=11, fontfamily='Arial')
        ax_cyl.set_ylabel('висота (м)', fontsize=11, fontfamily='Arial')
        ax_cyl.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)

        # ============================================
        # 3) Підрахунок підсумкових величин і вивід
        # ============================================

        # Площа циліндричної частини без нахлесту — рівно по довжині ремонту
        площа_cyl = Wrem * L

        # ============================================
        # ШВИ (циліндр + обидва днища) — ТОЧНИЙ РОЗРАХУНОК
        # ============================================

        # Допоміжна: довжина хорди кола на рівні y (для днищ)
        def _chord_len(y: float) -> float:
            if y <= -R or y >= R:
                return 0.0
            return 2.0 * math.sqrt(max(0.0, R*R - y*y))

        # ------------------------------
        # 1) ЦИЛІНДРИЧНА ЧАСТИНА
        # ------------------------------
        # Горизонтальні шви: нижній край + (міжрядові) + верхній край
        weld_h_cyl = (full_rows + 1) * Wrem

        # Вертикальні шви: у кожному ряду 2 крайових + (len(pattern) - 1) внутрішніх.
        def _row_visible_height(i: int) -> float:
            y_off = i * (h_smuha - overlap)
            if y_off >= L:
                return 0.0
            return min(h_smuha, L - y_off)

        weld_v_cyl = 0.0
        for i in range(full_rows):
            h_eff = _row_visible_height(i)
            if h_eff <= 0:
                break
            pat = patterns[i % len(patterns)]
            n_vertical = (len(pat) - 1) + 2   # внутрішні + 2 крайові
            weld_v_cyl += n_vertical * h_eff

        weld_cyl = weld_h_cyl + weld_v_cyl

        # ------------------------------
        # 2) ДНИЩЕ (ОДНЕ) — ЄДИНІ ШВИ
        # ------------------------------
        # Верхня межа ремонту днища (за Hcrit)
        y0 = -R
        y_top = -R + Hcrit  # тут ремонтуємо до цієї висоти

        # 2.1 Горизонтальні шви по хордах (без подвоєння міжсмугових):
        levels = []
        k = 0
        while True:
            y_level = y0 + k * h_smuha
            if y_level >= y_top - 1e-9:
                break
            levels.append(y_level)
            k += 1
        # додаємо верхню межу (може не збігатися з кратним кроком)
        levels.append(y_top)

        weld_h_bot_one = sum(_chord_len(y) for y in levels)

        # 2.2 Бічні шви (ліва+права кромки патча) — як дуги кола від y=-R до y=y_top
        theta_top = math.asin(max(-1.0, min(1.0, y_top / R)))   # ∈ [-π/2, π/2]
        arc_len_one_side = R * (theta_top - (-math.pi / 2))
        weld_side_bot_one = 2.0 * arc_len_one_side  # ліва+права

        # Шви одного днища:
        weld_bot_one = weld_h_bot_one + weld_side_bot_one

        # Обидва днища:
        weld_both_bottoms = 2.0 * weld_bot_one

        # ------------------------------
        # 3) ПІДСУМОК
        # ------------------------------
        total_weld = weld_cyl + weld_both_bottoms

        # Реально ремонтуєма площа (без нахлестів): обидва днища + циліндр
        реально_ремонт_площа = total_area_both_bottoms + площа_cyl

        # Фактична використана площа матеріалів (з нахлестом):
        # без бокового нахльосту; для верхньої смуги — без верхнього нахльосту
        площа_днищ_з_нахлестом_факт = 2 * sum(
            w * (h + (overlap/2 if idx < n_bot - 1 else 0.0))
            for idx, (w, h) in enumerate(zip(widths_bot, used_heights_bot))
        )

        фактична_площа_матеріалів = площа_днищ_з_нахлестом_факт + площа_cyl_з_нахлестом

        # Кількість листів 6×1.5 м
        sheet_area = 6 * 1.5
        total_sheets = math.ceil(фактична_площа_матеріалів / sheet_area)

        # — Ліва частина: основні підсумки —
        left_lines = [
            f"Площа одного днища: {cum_area_bot:.3f} м²",
            f"Площа обох днищ:    {total_area_both_bottoms:.3f} м²",
            f"Площа циліндричної частини: {площа_cyl:.3f} м²",
            f"Реально ремонтуєма площа: {реально_ремонт_площа:.3f} м²",
            f"Фактична використана площа матеріалів (з нахлестом {overlap_cm} см): {фактична_площа_матеріалів:.3f} м²",
            f"Заг. довжина швів:  {total_weld:.2f} м",
            f"Висота ремонтуємої частини резервуара: {Hcrit:.3f} м",
            f"Довжина кола: {circumference:.3f} м",
        ]

        # — Права частина: смуги та кількість листів —
        items = sorted(smuhaDict.items())
        right_lines = ["Розмір смуги     К-сть"]
        for key, cnt in items:
            right_lines.append(f"{key:>12s}   {cnt:>3d} шт")
        right_lines.append("")
        right_lines.append(f"К-сть листів (6×1.5м): {total_sheets} шт")
        right_lines.append(f"Використаний нахлест: {overlap_cm} см")

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
                if line.startswith("Реально ремонтуєма площа"):
                    st.markdown(f"<p style='color:blue;'>{line}</p>", unsafe_allow_html=True)
                elif line.startswith("Фактична використана площа матеріалів"):
                    st.markdown(f"<p style='color:red;'>{line}</p>", unsafe_allow_html=True)
                elif line.startswith("Заг. довжина швів"):
                    st.markdown(f"<p style='color:green;'>{line}</p>", unsafe_allow_html=True)
                else:
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