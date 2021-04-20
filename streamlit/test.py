import streamlit as st
import numpy as np

# Use the non-interactive Agg backend, which is recommended as a
# thread-safe backend.
# See https://matplotlib.org/3.3.2/faq/howto_faq.html#working-with-threads.
import matplotlib as mpl
mpl.use("agg")

import matplotlib.pyplot as plt


m = st.sidebar.slider("slope", min_value=0, max_value=10, value=5, step=1)
c = st.sidebar.slider("intercept", min_value=0., max_value=10., value=5., step=0.1)

# These two elements will be filled in later, so we create a placeholder
# for them using st.empty()
frame_text = st.sidebar.empty()
image = st.empty()



fig, ax = plt.subplots()

x = np.array([-1,1])

ax.plot(x, m*x+c)

ax.set_ylim([-5, 5])

plt.title("y=mx+c")
plt.xlabel("x")
plt.ylabel("y")


st.pyplot(fig)
