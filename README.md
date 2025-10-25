# Personalized Brain Decoding of Spontaneous Pain in Individuals With Chronic Pain

This repository is for sharing the codes that are used in the following paper.

> Jae-Joong Lee, Seongwoo Jo, Sungkun Cho, Choong-Wan Woo\*, Personalized Brain Decoding of Spontaneous Pain in Individuals With Chronic Pain, _in revision_

---

## System Requirements

### General Software

* **MATLAB**

### Experiment

* **Psychtoolbox 3.0**

### Preprocessing

* **FSL 6.0**
* **FreeSurfer 7.2**
* **AFNI 23.0**
* **Connectome Workbench 1.5.0**
* **ciftify**
  ➤ [https://github.com/edickie/ciftify](https://github.com/edickie/ciftify)

### Analysis

* **Connectome Workbench 1.5.0**
* **MRIcroGL 1.2**
* **MSCcodebase**
  ➤ [https://github.com/MidnightScanClub/MSCcodebase](https://github.com/MidnightScanClub/MSCcodebase)
* **cocoanCORE**
  ➤ [https://github.com/cocoanlab/cocoanCORE](https://github.com/cocoanlab/cocoanCORE)
* **Schaefer 500-region parcellation**
  ➤ [https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal)
* **Infomap 0.15.7**
  ➤ [https://github.com/mapequation/infomap](https://github.com/mapequation/infomap)

---

## Installation

Clone this repository:

```bash
git clone https://github.com/cocoanlab/DEIPP.git
```

---

## How to Run

| Stage         | Script                                                        | Notes                                       |
| ------------- | ------------------------------------------------------------- | ------------------------------------------- |
| Experiment    | `DEIPP.m`                                                     | Psychtoolbox-based fMRI task                |
| Preprocessing | `preprocessing_cocoan_anat.m` & `preprocessing_cocoan_func.m` | Anatomical and functional preprocessing     |
| Analysis      | `analysis_generate_command.m`                                 | Main analysis                               |


---

## Support / Questions

For inquiries or assistance, please contact:

* **Choong-Wan Woo** — [choongwan.woo@gmail.com](mailto:choongwan.woo@gmail.com)
* **Jae-Joong Lee** — [jaejoonglee92@gmail.com](mailto:jaejoonglee92@gmail.com)

