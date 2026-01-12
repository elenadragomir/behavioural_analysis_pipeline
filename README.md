# behavioural_analysis_pipeline: Closed-loop Phenotyping Framework
End-to-end Python pipeline for closed-loop behavioral tracking, feature engineering, and population-level statistical evaluation.
Currently optimized for OMR experiments in larval zebrafish, with modules for additional paradigms in development. 

## 🔬 System Components

### 1. Experimental Control & Stimulus Delivery (`/protocols`)
The Python scripts in this directory define the experimental logic and real-time closed-loop feedback systems.
* **Framework:** Powered by [Stytra](https://github.com/portugueslab/stytra), an open-source framework for visual stimulation.
* **Dynamic Protocols:** Includes contrast-sensitivity mapping (`free_omr_lr_cl_contrasts_v01.py`) and relative rotation stimulus delivery.
* **Closed-Loop Logic:** Implements real-time updates to stimulus orientation based on the animal's live heading to maintain consistent relative sensory input.

### 2. Behavioral Feature Extraction (`/analysis`)
These metrics are used to quantify how sensory input (visual gratings) is transformed into motor output (swimming bouts), allowing for the identification of subtle behavioral phenotypes across different groups.
The analytical pipeline is split into two stages to manage large-scale datasets efficiently.
* **Framework:** Utilizes [Bouter](https://github.com/portugueslab/bouter) for handling freely swimming experiment data structures.
* **Individual Processing (`01_feature_extraction.ipynb`):** Iterates through raw experiment folders to identify discrete behavioral "bouts" and extract features such as **Inter-Bout Interval (IBI)**, **Bout Duration**, and **Mean Velocity**.
* **Population Statistics (`02_population_statistics.ipynb`):** Aggregates data across cohorts to identify group-level trends. I leverage the [Bouter](https://github.com/portugueslab/bouter) package to extract higher-order metrics like **Heading Error** and **Angular Velocity** using circular statistics and noise reduction techniques to ensure signal clarity.
  
## 📊 Key Behavioral Features Engineered
* **Kinematics:** Displacement, instantaneous velocity, and turn-rate kinetics.
* **Decision Metrics:** Accuracy of response relative to stimulus direction (Heading Error).
* **Activity Patterns:** Frequency and temporal distribution of swimming bouts (IBIs) across different stimulus conditions.

  
## 📂 Shared Utility Library (`/shared_utils`)

> **Note:** This module is currently under active development to centralize logic across all behavioral paradigms.

I am migrating core analytical functions into a reusable library, `utilities_free_swim.py`, to ensure consistency as the framework expands. This library includes:

* **Bout Classification**: Logic to categorize movement into turns or forward swims based on angular displacement.
* **Spatial Metrics**: Functions to quantify environmental preference (**Time-in-Light**) and edge preference (**Thigmotaxis**).
* **Temporal Dynamics**: Calculations for average **Inter-Bout Intervals (IBI)** and detection of **Freezing Episodes**.
* **Data Quality**: Utilities like `time_out` to track and report instances where the animal is not detected (e.g., at the arena edge).


###  Data Infrastructure
* High-Frequency Processing: Designed to process high-frequency (100Hz+) tracking data 
* Storage formats: Managed via HDF5 for raw kinematic traces and JSON for metadata and population-level summaries.

