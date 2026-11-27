#%%
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.power import NormalIndPower
from statsmodels.stats.proportion import confint_proportions_2indep, proportion_confint

pd.set_option("display.width", 160)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 200)


@dataclass(frozen=True)
class Config:
    data_file: str = "data parasetamol jurnal.xlsx"
    verified_file: str = "data_terverifikasi.xlsx"
    output_dir: str = "outputs"
    seed: int = 20261127
    n_time: int = 20
    run_length: int = 10
    simulations: int = 500
    blocks: tuple = (5, 10)
    sbp_corrections: tuple = ((1113, 113), (1114, 114))
    draws: int = 200_000
    first_data_row: int = 5


CFG = Config()
RNG = np.random.default_rng(CFG.seed)
OUT = Path(CFG.output_dir)
OUT.mkdir(exist_ok=True)
VITALS = {"sbp": "Systolic BP, mmHg", "dbp": "Diastolic BP, mmHg", "hr": "Heart rate, beats/min", "spo2": "SpO2, %"}
ARMS = {"paracetamol": "Paracetamol", "kontrol": "Control"}
COLORS = {"Paracetamol": "#0072B2", "Control": "#D55E00"}
TABLES = {}


def section(title):
    print(f"\n{title}\n{'=' * len(title)}")


def keep(name, frame):
    TABLES[name] = frame
    return frame


def cols(var):
    return [f"{var}{t}" for t in range(1, CFG.n_time + 1)]


#%%
def locate_data():
    candidates = [Path("..") / "verification" / CFG.verified_file, Path(CFG.verified_file), Path(CFG.data_file),
                  Path("..") / "source" / CFG.data_file, Path("data") / CFG.data_file]
    for path in candidates:
        if path.exists():
            return path
    try:
        from google.colab import files
    except ImportError as error:
        raise FileNotFoundError(f"Place {CFG.data_file} next to the notebook or in ../source") from error
    uploaded = files.upload()
    return Path(next(iter(uploaded)))


def load_raw(path):
    book = openpyxl.load_workbook(path, data_only=True)
    names = ["excel_row", "arm", "no", "name", "age", "weight", "height", "diagnosis", "anaesthesia"]
    names += [c for var in VITALS for c in cols(var)] + ["wrench"]
    records = []
    for sheet in book.worksheets:
        for r in range(CFG.first_data_row, sheet.max_row + 1):
            values = [sheet.cell(r, c).value for c in range(1, len(names) - 1)]
            if any(v is not None for v in values):
                records.append([r, ARMS[sheet.title.strip().lower()]] + values)
    frame = pd.DataFrame(records, columns=names).drop(columns="name")
    frame.insert(0, "study_id", [f"P{i:03d}" for i in range(1, len(frame) + 1)])
    return frame


source_path = locate_data()
raw = load_raw(source_path)
section("Data import")
print(f"source file {source_path.name}")
print(f"records {len(raw)}  variables {raw.shape[1]}  direct identifiers removed at import")
print(raw.groupby("arm").size().to_string())

#%%
def to_number(series):
    text = series.map(lambda v: v.replace(",", ".").strip() if isinstance(v, str) else v)
    return pd.to_numeric(text, errors="coerce")


def clean(frame):
    data = frame.copy()
    log = []

    def note(rule, variable, mask):
        for _, row in data.loc[mask].iterrows():
            log.append({"rule": rule, "variable": variable, "study_id": row.study_id, "arm": row.arm,
                        "excel_row": row.excel_row, "raw_value": frame.loc[row.name, variable]})

    numeric = ["age", "weight", "height", "wrench"] + [c for var in VITALS for c in cols(var)]
    for column in numeric:
        note("decimal comma converted", column, frame[column].map(lambda v: isinstance(v, str) and "," in v))
        data[column] = to_number(data[column])
        note("missing in source", column, frame[column].isna())
    for column in cols("sbp"):
        for wrong, right in CFG.sbp_corrections:
            mask = data[column] == wrong
            note(f"team-reported entry error corrected to {right}", column, mask)
            data.loc[mask, column] = right
    limits = {"age": (15, 55), "weight": (35, 180), "height": (130, 200), "wrench": (0, 4)}
    limits.update({c: (50, 250) for c in cols("sbp")})
    limits.update({c: (25, 150) for c in cols("dbp")})
    limits.update({c: (35, 200) for c in cols("hr")})
    limits.update({c: (70, 100) for c in cols("spo2")})
    for column, (low, high) in limits.items():
        bad = data[column].notna() & ~data[column].between(low, high)
        note("implausible value set to missing", column, bad)
        data.loc[bad, column] = np.nan
    data["diagnosis"] = data.diagnosis.astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
    data["bmi"] = data.weight / (data.height / 100) ** 2
    data["shiver"] = np.where(data.wrench.notna(), (data.wrench >= 1).astype(float), np.nan)
    data["severe"] = np.where(data.wrench.notna(), (data.wrench >= 3).astype(float), np.nan)
    data["treated"] = (data.arm == "Paracetamol").astype(int)
    return data, pd.DataFrame(log)


df, cleaning_log = clean(raw)
section("Cleaning audit trail")
print(keep("cleaning_summary", cleaning_log.groupby(["rule", "variable"]).size().rename("n").reset_index()).to_string(index=False))
keep("cleaning_log", cleaning_log)
section("Completeness")
completeness = pd.DataFrame({v: [int(df[cols(v)].size), int(df[cols(v)].isna().sum().sum())] for v in VITALS},
                            index=["measurements", "missing"]).T
print(keep("completeness", completeness).to_string())
print(df[["age", "weight", "height", "bmi", "wrench"]].isna().groupby(df.arm).sum().to_string())

#%%
ELIGIBILITY_RULES = {
    "labour, ruptured membranes or failed induction": r"kpsw|kpd|lepsw|partur|gagal (?:drip|induksi)|drip gagal|drip$|gagal jam|kala|konprem|kontraksi",
    "gestation below 37 weeks": r"(?<!\d)(?:29-30|34-35|35-36|36-37) ?(?:mg+|minggu)|(?<!\d)(?:34|36) ?(?:mg+|minggu)",
    "twin pregnancy": r"gemel",
    "severe pre-eclampsia or chronic hypertension (ASA II to be confirmed)": r"peb|pes$|ped$|hipertensi",
    "placenta praevia or fetal demise or anomaly (review)": r"previa|iuvd|iuf$|\+ ?iuf\b|edward",
}


def eligibility(frame):
    text = frame.diagnosis.str.lower()
    flags = pd.DataFrame({label: text.str.contains(rule, regex=True) for label, rule in ELIGIBILITY_RULES.items()})
    flags["age above 45 years"] = frame.age > 45
    return flags


flags = eligibility(df)
core = ["labour, ruptured membranes or failed induction", "gestation below 37 weeks", "twin pregnancy", "age above 45 years"]
df["protocol_ineligible"] = flags[core].any(axis=1)
section("Eligibility screen from free-text diagnosis (requires source verification)")
counts = flags.groupby(df.arm).sum().T
counts.loc["any core criterion, unique women"] = df.groupby("arm").protocol_ineligible.sum()
print(keep("eligibility_counts", counts).to_string())
listing = df.loc[flags.any(axis=1), ["study_id", "arm", "excel_row", "age", "diagnosis"]].copy()
listing["criteria"] = flags.loc[listing.index].apply(lambda r: "; ".join(r.index[r.values]), axis=1)
keep("eligibility_listing", listing)
print(f"\nwomen with at least one flag {len(listing)}; row-level listing is exported locally and not printed")

#%%
def sequence_matrix(var, start=1):
    values = raw[cols(var)[start - 1:]].apply(to_number)
    return values.fillna(-1).to_numpy().astype(int)


def exact_members(matrix):
    keys = list(map(tuple, matrix))
    tally = Counter(keys)
    return np.array([tally[k] > 1 for k in keys]), keys, tally


def run_members(matrix, length):
    members = set()
    for s in range(matrix.shape[1] - length + 1):
        flagged = exact_members(matrix[:, s:s + length])[0]
        members |= set(np.flatnonzero(flagged))
    return members


def simulate_ar1(matrix, n_sim):
    n, t = matrix.shape
    mean = matrix.mean(0)
    residual = matrix - mean
    person = residual.mean(1)
    within = residual - person[:, None]
    rho = np.sum(within[:, 1:] * within[:, :-1]) / np.sum(within[:, :-1] ** 2)
    for _ in range(n_sim):
        noise = RNG.normal(0, within.std() * np.sqrt(1 - rho ** 2), (n, t))
        series = np.empty((n, t))
        series[:, 0] = RNG.normal(0, within.std(), n)
        for j in range(1, t):
            series[:, j] = rho * series[:, j - 1] + noise[:, j]
        yield np.rint(mean + RNG.normal(0, person.std(), n)[:, None] + series).astype(int)


def simulate_blocks(matrix, n_sim, width):
    n, t = matrix.shape
    for _ in range(n_sim):
        yield np.column_stack([matrix[RNG.integers(0, n, n), s:s + width] for s in range(0, t, width)])


def null_summary(observed, matrices):
    exact, runs = [], []
    for sample in matrices:
        exact.append(exact_members(sample)[0].sum())
        runs.append(len(run_members(sample, CFG.run_length)))
    exact = np.array(exact)
    return {"duplicated mean": round(float(exact.mean()), 1), "duplicated max": int(exact.max()),
            "run mean": round(float(np.mean(runs)), 1), "run max": int(np.max(runs)),
            "Monte Carlo P (duplicated)": round(float((1 + np.sum(exact >= observed)) / (1 + len(exact))), 4)}


section("Integrity screen 1: repeated individual time series (raw values, SpO2 excluded)")
rows, groups, nulls = [], [], []
for var in ["sbp", "dbp", "hr"]:
    matrix = sequence_matrix(var)
    member, keys, tally = exact_members(matrix)
    across = {k for k in tally if df.arm[[key == k for key in keys]].nunique() == 2}
    rows.append({"variable": VITALS[var], "patients": len(df), "with an identical 20-point series": int(member.sum()),
                 "with an identical T2 to T20 series": int(exact_members(sequence_matrix(var, 2))[0].sum()),
                 "largest identical group": max(tally.values()), "series present in both arms": len(across),
                 f"sharing {CFG.run_length}+ values at the same time points": len(run_members(matrix, CFG.run_length))})
    complete = matrix[(matrix >= 0).all(1) & (matrix < 1000).all(1)]
    observed = int(exact_members(complete)[0].sum())
    models = {"AR(1) with patient effect": simulate_ar1(complete, CFG.simulations)}
    models.update({f"{w}-point block resampling of these data": simulate_blocks(complete, CFG.simulations, w) for w in CFG.blocks})
    for model, matrices in models.items():
        nulls.append({"variable": var, "model": model, "complete records": len(complete), "observed duplicated": observed,
                      "observed run": len(run_members(complete, CFG.run_length)), **null_summary(observed, matrices)})
    for g, (key, size) in enumerate([x for x in tally.most_common() if x[1] > 1], 1):
        hit = df.loc[[k == key for k in keys]]
        groups.append({"variable": var, "group": g, "size": size,
                       "paracetamol_excel_rows": ", ".join(map(str, hit.excel_row[hit.arm == "Paracetamol"])),
                       "control_excel_rows": ", ".join(map(str, hit.excel_row[hit.arm == "Control"])),
                       "series": " ".join(map(str, key))})
integrity = keep("integrity_series", pd.DataFrame(rows))
keep("integrity_null_models", pd.DataFrame(nulls))
keep("integrity_copy_groups", pd.DataFrame(groups))
print(integrity.T.to_string(header=False))
print("\nnull models fitted to records without missing or four-digit entries")
print(TABLES["integrity_null_models"].to_string(index=False))
exact_any = np.column_stack([exact_members(sequence_matrix(v))[0] for v in ["sbp", "dbp", "hr"]])
run_any = set().union(*[run_members(sequence_matrix(v), CFG.run_length) for v in ["sbp", "dbp", "hr"]])
print(f"\nexact duplicate in at least one pressure or heart-rate series {exact_any.any(1).sum()} of {len(df)}")
print(f"exact duplicate in all three series {exact_any.all(1).sum()} of {len(df)}")
print(f"patients with at least one shared run {len(run_any)} of {len(df)}")

#%%
section("Integrity screen 2: records, recurrent errors, digits and outcome coding")
profile = df[[c for v in VITALS for c in cols(v)]].fillna(-1)
twins = df.loc[profile.duplicated(keep=False), ["study_id", "arm", "excel_row", "age", "weight", "height", "diagnosis", "wrench"]]
keep("duplicate_records", twins)
print(f"records sharing all 80 haemodynamic values {len(twins)}")
print(twins.groupby("arm").size().to_string())
errors = cleaning_log[cleaning_log.raw_value.isin([1113, 1114])]
print("\nrecurrent four-digit systolic entries")
print(keep("four_digit_corrections", errors.groupby(["raw_value", "variable", "arm"]).size().rename("n").reset_index()).to_string(index=False))
digits = pd.DataFrame({v: np.bincount(df[cols(v)].stack().astype(int) % 10, minlength=10) for v in ["sbp", "dbp", "hr"]})
print("\nterminal digit frequencies")
print(keep("terminal_digits", digits).T.to_string())
grades = pd.crosstab(df.arm, df.wrench.fillna(-1)).reindex(columns=[-1, 0, 1, 2, 3, 4], fill_value=0)
grades.columns = ["missing", "grade 0", "grade 1", "grade 2", "grade 3", "grade 4"]
print("\nrecorded Wrench grades")
print(keep("wrench_grades", grades).to_string())

#%%
ar1 = TABLES["integrity_null_models"].query("model.str.startswith('AR')", engine="python")
gate = pd.DataFrame([
    ("integrity", "Repeated vital-sign series within the AR(1) null range", bool((ar1["observed duplicated"] <= ar1["duplicated max"]).all())),
    ("integrity", "No complete haemodynamic record is repeated", len(twins) == 0),
    ("integrity", "No recurrent identical impossible value", len(errors) < 2),
    ("plausibility", "All five Wrench grades observed", bool((grades.iloc[:, 1:].sum() > 0).all())),
    ("plausibility", "No eligibility flag from free-text diagnosis", int(df.protocol_ineligible.sum()) == 0),
], columns=["type", "check", "passed"])
DATA_VERIFIED = bool(gate.passed[gate.type == "integrity"].all())
STATUS = "verified" if DATA_VERIFIED else "UNVERIFIED"
section("Integrity gate")
print(keep("integrity_gate", gate).to_string(index=False))
print("\nSTATUS:", "PASSED" if DATA_VERIFIED else "FAILED. Every estimate below is unverified and must not be published before source-record verification.")

#%%
def heatmap(var, ax):
    matrix = sequence_matrix(var)
    member, keys, tally = exact_members(matrix)
    order = sorted(range(len(df)), key=lambda i: (-tally[keys[i]], keys[i]))
    image = ax.imshow(np.where(matrix[order] < 0, np.nan, matrix[order]), aspect="auto", cmap="viridis", interpolation="nearest")
    ax.axhline(member.sum() - 0.5, color="red", lw=1)
    ax.set(title=VITALS[var], xlabel="Time point")
    plt.colorbar(image, ax=ax, shrink=0.7)


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.show()


fig, axes = plt.subplots(1, 3, figsize=(13, 6), constrained_layout=True)
for var, ax in zip(["sbp", "dbp", "hr"], axes):
    heatmap(var, ax)
axes[0].set_ylabel("Patients, ordered by size of identical group (above red line: duplicated)")
save(fig, "fig_integrity_heatmap")

fig, axes = plt.subplots(2, 3, figsize=(12, 6.5), constrained_layout=True)
for ax, (var, label) in zip(axes.ravel(), [("age", "Age, years"), ("weight", "Weight, kg"), ("height", "Height, cm"), ("bmi", "BMI, kg/m2")]):
    for arm, colour in COLORS.items():
        ax.hist(df.loc[df.arm == arm, var].dropna(), bins=15, alpha=0.55, color=colour, label=arm)
    ax.set(xlabel=label, ylabel="Women")
    ax.legend(frameon=False)
digits.plot.bar(ax=axes[1, 1], rot=0)
axes[1, 1].set(xlabel="Terminal digit", ylabel="Measurements")
(grades.iloc[:, 1:].div(grades.iloc[:, 1:].sum(1), axis=0) * 100).T.plot.bar(ax=axes[1, 2], rot=30, color=[COLORS[a] for a in grades.index])
axes[1, 2].set(ylabel="Percent")
save(fig, "fig_eda")

#%%
def smd(a, b):
    pooled = np.sqrt((a.var() + b.var()) / 2)
    return (a.mean() - b.mean()) / pooled if pooled > 0 else 0.0


def describe_continuous(name, label):
    a, b = (df.loc[df.arm == arm, name].dropna() for arm in COLORS)
    fmt = lambda x: f"{x.mean():.1f} ({x.std():.1f}); {x.median():.1f} [{x.quantile(.25):.1f}, {x.quantile(.75):.1f}]; n={len(x)}"
    return {"characteristic": label, "Paracetamol": fmt(a), "Control": fmt(b), "SMD": round(smd(a, b), 2),
            "Shapiro P": f"{min(stats.shapiro(a).pvalue, stats.shapiro(b).pvalue):.3g}"}


def describe_binary(mask, label):
    k = mask.groupby(df.arm).sum()
    n = df.groupby("arm").size()
    p = k / n
    pooled = np.sqrt((p * (1 - p)).mean())
    return {"characteristic": label, **{a: f"{int(k[a])} ({100 * p[a]:.1f}%)" for a in COLORS},
            "SMD": round(float((p["Paracetamol"] - p["Control"]) / pooled), 2) if pooled > 0 else 0.0, "Shapiro P": ""}


text = df.diagnosis.str.lower()
df["previous_cs"] = text.str.contains(r"bekas\s*\+?\s*sc|bekas sc")
df["cpd"] = text.str.contains(r"cpd|cpp|panggul sempit")
df["malpresentation"] = text.str.contains(r"sunsang|sungsang|let su|lintang|letak muka|letak kaki")
gravidity = text.str.extract(r"^g\s*([0-9i]+)")[0].replace({"i": "1", "22": "2"})
df["gravidity"] = pd.to_numeric(gravidity, errors="coerce")
baseline = [describe_continuous(v, l) for v, l in [
    ("age", "Age, years"), ("weight", "Weight, kg"), ("height", "Height, cm"), ("bmi", "BMI, kg/m2"),
    ("sbp1", "Baseline systolic BP, mmHg"), ("dbp1", "Baseline diastolic BP, mmHg"), ("hr1", "Baseline heart rate"), ("spo21", "Baseline SpO2, %")]]
baseline += [describe_binary(m, l) for m, l in [
    (df.gravidity == 1, "Gravidity 1"), (df.gravidity == 2, "Gravidity 2"), (df.gravidity >= 3, "Gravidity 3 or more"),
    (df.gravidity.isna(), "Gravidity not recorded"), (df.previous_cs, "Previous caesarean delivery"),
    (df.cpd, "Cephalopelvic disproportion"), (df.malpresentation, "Malpresentation"),
    (df.protocol_ineligible, "Did not meet eligibility criteria")]]
section("Baseline characteristics")
print(keep("table1", pd.DataFrame(baseline)).to_string(index=False))
