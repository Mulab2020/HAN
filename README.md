# Code for *A Thalamus–Brainstem Attractor Network Drives History-Biased Decisions*

This repository contains code related to the paper:

> **A Thalamus–Brainstem Attractor Network Drives History-Biased Decisions**

The code provided here focuses on the **data preprocessing, figure plotting, and model simulation part** of the study.
The scripts for model simulation of the hierarchical attractor network are implemented in Python (using Jupyter notebooks) and are organized to reproduce specific analyses and visualizations of neural activity trajectories. All other scripts are in MATLAB.

## Contents
- **`Behavioral_model`**  
  Code for recapitulating the behavioral model (reactive agent and serial-dependent agent) in the paper.
  
- **`Ephys_data_preprocess`**  
  Ephys signal filtering and bout extraction codes.
  
- **`Imaging_data_preprocess`**  
  Preprocessing code for whole-brain calcium imaging.
  
- **`Plot_code`**  
  Code for drawing the figures in the manuscript.
  
- **`trajectory_plot.ipynb`**  
  Code for plotting average neural activity trajectories across different conditions.

- **`tranInh_trial_plot.ipynb`**  
  Code for analyzing and plotting trial-level inhibitory/transfer dynamics.

- **`example_trial_plot.ipynb`**  
  Example script showing single-trial activity and plotting pipeline.

## Notes

- The folder structure may evolve as additional scripts and analysis modules are included.  
- Dependencies: Python 3.11, Jupyter Notebook, NumPy, Matplotlib, BrainPy (and other standard scientific Python packages); MATLAB

