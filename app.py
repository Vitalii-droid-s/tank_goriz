import streamlit as st
import matplotlib.pyplot as plt
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages

st.set_page_config(layout="centered")
st.title("Тест кнопки PDF")

# Створимо простий графік
fig, ax = plt.subplots()
ax.plot([1, 2, 3], [1, 4, 9])
st.pyplot(fig)

# Кнопка PDF
st.subheader("⬇️ Збереження PDF")
if st.button("Створити PDF"):
    buffer = BytesIO()
    with PdfPages(buffer) as pdf:
        pdf.savefig(fig, bbox_inches='tight')

    st.download_button(
        label="📄 Завантажити PDF-файл",
        data=buffer.getvalue(),
        file_name="резервуар_розрахунок.pdf",
        mime="application/pdf"
    )
