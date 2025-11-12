import matplotlib.pyplot as plt

from shiny.express import  input, render, ui
import numpy as np
import seaborn as sns
import math


ui.page_opts(title="Delayed-normalization model demo")

def gammaPDF(t, tau, n):
    pdf = (t / tau) ** (n - 1) * np.exp(-t / tau) / (tau * math.factorial(n-1))
    pdf = pdf / np.sum(pdf)
    return pdf

with ui.layout_sidebar():
    with ui.sidebar():
        with ui.accordion():
            with ui.accordion_panel("Settings"):
                ui.input_slider("s", ui.tags.span("Stimulus duration (s)", style="font-size: smaller;"), 0, 1.2, 0.4, step=0.01)
            with ui.accordion_panel("Free parameters"):
                ui.input_slider("tau1", ui.tags.span("tau1", style="font-size: smaller;"), 0.001, 0.3, 0.05, step=0.01)
                ui.input_slider("tau2", ui.tags.span("tau2", style="font-size: smaller;"), 0.001, 0.5, 0.1, step=0.01)
                ui.input_slider("w", ui.tags.span("weight", style="font-size: smaller;"), 0.0, 0.8, 0.0, step=0.01)
                ui.input_slider("n", ui.tags.span("n", style="font-size: smaller;"), 0.5, 2.5, 1.0, step=0.1)
                ui.input_slider("sigma", ui.tags.span("sigma", style="font-size: smaller;"), 0.01, 0.5, 0.1, step=0.01)

    @render.plot(alt="Predicted neural response")  
    def plot():  
        tau1 = input.tau1()
        tau2 = input.tau2()
        w = input.w()
        n = input.n()
        sigma = input.sigma()
        
        # One-pulse stimulus
        dur = input.s() # in seconds
        samples  = 1000;
        finer_t = np.arange(0, 3*samples, 1) # Max out to 3 seconds
        t_axis = finer_t / samples
        stimSeq = np.zeros(len(finer_t))
        stimSeq[(finer_t > 0) & (finer_t <= dur * samples)] = 1

        # Impulse response function
        irf_pos = gammaPDF(t_axis, tau1, 2)
        irf_neg = gammaPDF(t_axis, tau1*1.5, 2)
        irf = irf_pos - w * irf_neg

        # Delayed response function
        lowpass = np.exp(-t_axis / tau2)
        lowpass = lowpass / np.sum(lowpass)

        # Numerator
        linrsp = np.convolve(stimSeq, irf, mode='full')
        linrsp = linrsp[0:len(t_axis)]
        numrsp = np.abs(linrsp) ** n

        # Denominator
        poolrsp = np.convolve(linrsp, lowpass, mode='full')
        poolrsp = poolrsp[0:len(t_axis)]
        demrsp = sigma ** n + np.abs(poolrsp) ** n
        
        # Normalization response
        normrsp = np.divide(numrsp, demrsp)

        # Plot the results

        sns.set_palette("colorblind") 
        colors = sns.color_palette()

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(t_axis, normrsp, color=colors[0], label='Predicted neural response')
        ax.plot(t_axis, stimSeq, color=colors[1], label='Stimulus', linestyle='--')
        ax.set_xlim([-0.1, 1.5])
        ax.set_ylim([0, 3])
        ax.legend(loc='upper right')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Response')
\
        return fig  
