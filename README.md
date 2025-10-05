# 🌍 SpaceApps Challenge -- Asteroid Impact Simulator

## Overview

This project was developed for the **NASA Space Apps Challenge**.\
It is an **interactive, educational web app** built with
[Streamlit](https://streamlit.io) that allows K-12 students and the
public to explore what happens when an asteroid impacts Earth.

The app provides an approachable interface to: - Experiment with
asteroid parameters (size, speed, composition, angle, location). -
Visualize outcomes such as **energy released**, **airburst vs crater
formation**, **blast effects**, and **seismic impact**. - Compare
"before" and "after" outcomes of hypothetical **deflection
strategies**. - Learn about the **science of impacts** through simple
explanations and key equations.

The focus is on **education, clarity, and safety**:
- Sensitive outputs (like human impacts) are hidden by default, but can
be toggled for advanced discussion.\
- Confidence badges show which outputs are **high confidence
(formula-based)** vs **simplified approximations**.

------------------------------------------------------------------------

## ✨ Features

-   **Explore Tab** -- Adjust asteroid properties with sliders and
    instantly see results.\
-   **Defend Earth Tab** -- Simulate deflection strategies and compare
    side-by-side impacts.\
-   **Learn Tab** -- Age-appropriate explanations of impact science and
    important equations.\
-   **Toggle Confidence Levels** -- Teachers/students can show/hide
    accuracy badges.\
-   **K-12 Friendly Mode** -- Simplified language, no graphic details,
    clear disclaimers.\
-   **NASA Feeds Integration** -- Pulls in latest asteroid-related news
    and events.

------------------------------------------------------------------------

## 🧮 Scientific Basis

The simulator is inspired by real research including: - The [Earth
Impact Effects Program (Collins, Melosh, Marcus,
2005)](https://impact.ese.ic.ac.uk/ImpactEarth/). - Published
NASA/planetary defense white papers.\
- Simplified blast scaling laws and seismic energy conversion.

**What's accurate**\
- Mass and kinetic energy are direct physics calculations.\
- Airburst vs crater thresholds are modeled after established
atmospheric breakup studies.

**What's approximate**\
- Blast radii and seismic magnitudes are simplified scaling laws (good
for order-of-magnitude).\
- Population exposure is based on rough density overlays, not real-time
demographics.

------------------------------------------------------------------------

## 🚀 Tech Stack

-   [Python 3.11+](https://www.python.org/)\
-   [Streamlit](https://streamlit.io) -- interactive web app framework\
-   [Pandas](https://pandas.pydata.org/) -- data handling\
-   [Matplotlib](https://matplotlib.org/) -- plotting\
-   [Requests](https://docs.python-requests.org/) -- NASA data feeds

------------------------------------------------------------------------

## 📦 Installation & Usage

### Prerequisites

-   Python 3.11 or later
-   `pip` package manager

### Clone the repository

``` bash
git clone https://github.com/<your-team>/SpaceAppsChallenge.git
cd SpaceAppsChallenge
```

### Install dependencies

``` bash
pip install -r requirements.txt
```

### Run the app

``` bash
streamlit run main.py
```

The app will start on `http://localhost:8501`.

------------------------------------------------------------------------

## 📚 Educational Use

Teachers can use this app to: - Demonstrate the importance of asteroid
monitoring and planetary defense.\
- Discuss physics concepts: mass, velocity, energy, scaling laws.\
- Explore real-world NASA data and connect science to current events.

Built-in safeguards (like toggling sensitive outputs) make it suitable
for **classroom discussions** at all levels.

------------------------------------------------------------------------

## 🏆 Hackathon Context

-   **Challenge:** Planetary Defense & Asteroid Impacts\
-   **Event:** NASA Space Apps Challenge\
-   **Team Name:** `<Your Team Name>`\
-   **Year:** 2025

We wanted to make asteroid science accessible and engaging for
**students and the general public**, while staying grounded in real
physics.

------------------------------------------------------------------------

## 📸 Screenshots

(Add screenshots of Explore Tab, Defend Earth Tab, Learn Tab here)

------------------------------------------------------------------------

## 🔮 Future Work

-   Integrate with NASA's **real asteroid database (NeoWs API)** for
    live objects.\
-   Add **3D crater and impact visualizations**.\
-   Expand language support for global classrooms.\
-   Improve scaling laws with more precise models.

------------------------------------------------------------------------

## 🤝 Acknowledgments

-   NASA Space Apps Challenge for the inspiration.\
-   Earth Impact Effects Program for scientific reference.\
-   Open-source libraries: Streamlit, Pandas, Matplotlib.\
-   All mentors, volunteers, and teammates!

------------------------------------------------------------------------

## 📜 License

This project is licensed under the MIT License -- see the
[LICENSE](LICENSE) file for details.
