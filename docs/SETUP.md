# Setup without git

How to get the project code and run it on your computer, without using git. You only download files from GitHub. Tasks are in `docs/TEAM_TASKS.md`; how to hand in your results is under "Handoff rules" there.

**Bohyeon:** you don't need the EEG data or most of these steps. Your inputs are CSV files (see [Only need the CSVs?](#only-need-the-csvs) at the end).

---

## 1. Download the project

1. Open the repository page on GitHub.
2. Click the green **Code** button, then **Download ZIP**.
3. Unzip it somewhere convenient. You'll get a folder called `AI-based-Neural-Data-Science-main`. This is the "project folder" below.

**When Chuong announces an update,** download the ZIP again and unzip it to a new folder. Copy your `paths.data_root` line (step 3) into the new `config.yaml`, because the downloaded copy has Chuong's path.

To get a single file instead (for example one CSV), open it on GitHub and click the **Download raw file** button (the down-arrow icon at the top right of the file view).

## 2. Install Python and the packages

The project was tested with **Python 3.10**; 3.11 and 3.12 should also work.

### macOS

1. Install Python 3.10 from <https://www.python.org/downloads/macos/> (or `brew install python@3.10` if you use Homebrew).
2. Open **Terminal**, go to the project folder, and create an environment:

   ```bash
   cd ~/Downloads/AI-based-Neural-Data-Science-main      # your project folder
   python3.10 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   pip install -r requirements.txt
   ```

### Windows

1. Install Python 3.10 from <https://www.python.org/downloads/windows/>. In the installer, tick **"Add python.exe to PATH"**.
2. Open **PowerShell**, go to the project folder, and create an environment:

   ```powershell
   cd $HOME\Downloads\AI-based-Neural-Data-Science-main   # your project folder
   py -3.10 -m venv .venv
   .venv\Scripts\Activate.ps1
   python -m pip install --upgrade pip
   pip install -r requirements.txt
   ```

   If PowerShell refuses to run `Activate.ps1`, run this once and try again:
   `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`

Installing takes a few minutes. **Every time you open a new terminal,** go to the project folder and activate the environment again (`source .venv/bin/activate` on macOS, `.venv\Scripts\Activate.ps1` on Windows) before running anything.

## 3. Point the code at your dataset

Open `config.yaml` in a text editor and change `paths.data_root` to the folder on your computer that contains `subjects.csv`, `continuous_data/`, `epoched_data/` and `trial_info/`:

```yaml
paths:
  data_root: "/Users/yourname/Data/Dataset_ALL"        # macOS
  # data_root: "C:/Users/yourname/Data/Dataset_ALL"    # Windows: use forward slashes
```

On Windows, write the path with **forward slashes** (`C:/Users/...`). Backslashes inside double quotes are read as escape characters and break the path.

Don't change anything else in `config.yaml`: all analysis rules live there and are agreed by the team.

## 4. Check that it works

From the project folder, with the environment active:

```bash
python notebooks/example_backfit.py
```

This loads `config.yaml` and the K = 5 templates, keeps sub-P001's correct Simon trials, filters the cue-locked epochs at 2–20 Hz, resamples them to 250 Hz, back-fits the templates and prints each state's share of time 0–800 ms after the stimulus. Expected output:

```text
Templates K = 5: states ['A', 'B', 'C', 'D', 'E']
sub-P001: 322 of 384 trials kept
Labels: (322, 875)
  A: 11.4%
  B: 11.0%
  C: 20.3%
  D: 29.0%
  E: 28.4%
```

Add a subject ID to try another one, for example `python notebooks/example_backfit.py sub-C004`. The script is short and commented; copy it as the starting point for your own analysis.

**In your own code**, load things the same way:

```python
import sys
sys.path.insert(0, "path/to/project/folder")
from src import io as pio, microstates as ms

cfg = pio.load_config()
templates = ms.load_templates(cfg["microstates"]["primary_k"], cfg)
```

Notebooks in `notebooks/` are plain Python files with `# %%` cell markers. In VS Code they run cell by cell ("Run Cell"); from a terminal, run them with `python notebooks/<file>.py`.

## Only need the CSVs?

For analyses on the feature tables (Bohyeon's RQ1 and RQ3 tasks), no EEG data and no templates are needed:

1. Download the CSVs from `results/features/` on GitHub (or the whole ZIP): `microstates_task.csv` (primary), the sensitivity files (`_nosmooth`, `_qcincluded`, `_K4`, `microstates_rest_noblink.csv`) and `behavior.csv`.
2. Group, age, sex, MoCA and UPDRS are **not** in these files. They come from `subjects.csv` in your copy of the dataset (`Dataset_ALL/subjects.csv`; a small file, no EEG). It is not on GitHub because the dataset is never committed. Merge on `subject_id`.
3. Install the packages as in step 2 (pandas, statsmodels and scikit-learn are what you'll use), or use any Python environment that has them.
