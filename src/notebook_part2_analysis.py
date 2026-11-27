#%%
def two_by_two(frame, outcome):
    use = frame.dropna(subset=[outcome])
    g = use.groupby("treated")[outcome].agg(["sum", "size"]).reindex([1, 0])
    return int(g.loc[1, "sum"]), int(g.loc[1, "size"]), int(g.loc[0, "sum"]), int(g.loc[0, "size"])


def effects(a, n1, c, n0):
    table = np.array([[a, n1 - a], [c, n0 - c]])
    p1, p0 = a / n1, c / n0
    rr_lo, rr_hi = confint_proportions_2indep(a, n1, c, n0, compare="ratio", method="log")
    or_lo, or_hi = confint_proportions_2indep(a, n1, c, n0, compare="odds-ratio", method="logit")
    rd_lo, rd_hi = confint_proportions_2indep(c, n0, a, n1, compare="diff", method="newcomb")
    arr = p0 - p1
    ci1 = proportion_confint(a, n1, method="wilson")
    ci0 = proportion_confint(c, n0, method="wilson")
    chi = stats.chi2_contingency(table, correction=False)
    return {"paracetamol": f"{a}/{n1} ({100 * p1:.1f}%; {100 * ci1[0]:.1f} to {100 * ci1[1]:.1f})",
            "control": f"{c}/{n0} ({100 * p0:.1f}%; {100 * ci0[0]:.1f} to {100 * ci0[1]:.1f})",
            "RR (95% CI)": f"{p1 / p0:.2f} ({rr_lo:.2f} to {rr_hi:.2f})",
            "ARR % (95% CI)": f"{100 * arr:.1f} ({100 * rd_lo:.1f} to {100 * rd_hi:.1f})",
            "NNT (95% CI)": f"{1 / arr:.1f} ({1 / rd_hi:.1f} to {1 / rd_lo:.1f})" if rd_lo > 0 else "not defined",
            "OR (95% CI)": f"{(a * (n0 - c)) / ((n1 - a) * c):.2f} ({or_lo:.2f} to {or_hi:.2f})",
            "chi-square": round(chi.statistic, 1), "P chi-square": f"{chi.pvalue:.2e}",
            "P Fisher": f"{stats.fisher_exact(table).pvalue:.2e}"}


def scenario(frame, value_treated, value_control):
    filled = frame.copy()
    fill = pd.Series(np.where(filled.treated == 1, value_treated, value_control), index=filled.index)
    filled["shiver"] = filled.shiver.fillna(fill)
    return filled


analyses = {
    "Primary: all recorded women with an outcome": df,
    "Missing outcome counted as no shivering": scenario(df, 0, 0),
    "Missing outcome counted as shivering": scenario(df, 1, 1),
    "Excluding women not meeting eligibility criteria": df[~df.protocol_ineligible],
    "Also excluding fetal distress or PTM notes (post hoc)": df[~df.protocol_ineligible & ~df.diagnosis.str.lower().str.contains(r"gawat|vital|\bptm\b")],
}
primary = pd.DataFrame({k: effects(*two_by_two(v, "shiver")) for k, v in analyses.items()})
section("Primary outcome: any shivering (unverified)" if not DATA_VERIFIED else "Primary outcome: any shivering")
print(keep("primary_outcome", primary.T.reset_index(names="analysis").assign(status=STATUS)).set_index("analysis").T.to_string())
print("\nA per-protocol analysis cannot be run because adherence to dose, timing and protocol steps was not recorded.")

#%%
def fragility(a, n1, c, n0):
    flips = 0
    while a < n1 and stats.fisher_exact([[a, n1 - a], [c, n0 - c]]).pvalue < 0.05:
        a, flips = a + 1, flips + 1
    return flips


def poisson_rr(formula, frame):
    terms = [t.strip() for t in formula.split("~")[1].split("+")]
    use = frame.dropna(subset=["shiver"] + terms)
    fit = smf.glm(formula, use, family=sm.families.Poisson()).fit(cov_type="HC1")
    lo, hi = np.exp(fit.conf_int().loc["treated"])
    return {"model": formula, "n": len(use), "RR": round(np.exp(fit.params["treated"]), 3), "low": round(lo, 3),
            "high": round(hi, 3), "P": f"{fit.pvalues['treated']:.2e}"}


a, n1, c, n0 = two_by_two(df, "shiver")
section("Robustness")
adjusted = pd.DataFrame([poisson_rr(f, df) for f in [
    "shiver ~ treated", "shiver ~ treated + age + bmi", "shiver ~ treated + age + bmi + previous_cs + sbp1 + hr1"]])
print(keep("adjusted_models", adjusted.assign(status=STATUS)).to_string(index=False))
post1 = RNG.beta(1 + a, 1 + n1 - a, CFG.draws)
post0 = RNG.beta(1 + c, 1 + n0 - c, CFG.draws)
bayes = pd.Series({"posterior median RR": np.median(post1 / post0),
                   "95% credible low": np.quantile(post1 / post0, .025), "95% credible high": np.quantile(post1 / post0, .975),
                   "P(RR < 1)": np.mean(post1 < post0), "P(ARR > 15 points)": np.mean(post0 - post1 > .15),
                   "fragility index": fragility(a, n1, c, n0)})
print(keep("bayes_fragility", bayes.rename("value").to_frame().assign(status=STATUS)).to_string())

subgroups = []
for label, flag in [("Age 35 or older", df.age >= 35), ("BMI 30 or higher", (df.bmi >= 30).astype(float).where(df.bmi.notna())), ("Previous caesarean", df.previous_cs)]:
    work = df.assign(sub=flag.astype(float)).dropna(subset=["shiver", "sub"])
    fit = smf.glm("shiver ~ treated * sub", work, family=sm.families.Poisson()).fit(cov_type="HC1")
    for level in (0, 1):
        res = effects(*two_by_two(work[work["sub"] == level], "shiver"))
        subgroups.append({"subgroup": label, "present": bool(level), "RR (95% CI)": res["RR (95% CI)"],
                          "interaction P": round(fit.pvalues["treated:sub"], 3)})
print(keep("subgroups_exploratory", pd.DataFrame(subgroups).assign(status=STATUS)).to_string(index=False))

#%%
section("Secondary outcomes: severity")
x, y = (df.loc[df.arm == arm, "wrench"].dropna() for arm in COLORS)
mw = stats.mannwhitneyu(x, y, alternative="two-sided")
delta = 2 * mw.statistic / (len(x) * len(y)) - 1
table3x2 = pd.crosstab(df.arm, df.wrench)
print(f"Mann-Whitney P {mw.pvalue:.2e}  Cliff delta {delta:.3f}  probabilistic index {mw.statistic / (len(x) * len(y)):.3f}")
print(f"chi-square for observed grades P {stats.chi2_contingency(table3x2).pvalue:.2e}")
severity = pd.DataFrame({
    "grade 3 or 4, all women": effects(*two_by_two(df, "severe")),
    "grade 3 or 4 among shiverers (post-randomisation subset)": effects(*two_by_two(df[df.shiver == 1], "severe"))})
print(keep("severity", severity.T.reset_index(names="analysis").assign(status=STATUS)).set_index("analysis").T.to_string())

#%%
def long_form(var):
    wide = df[["study_id", "arm", "treated"] + cols(var)]
    out = wide.melt(["study_id", "arm", "treated"], var_name="time", value_name="value").dropna()
    out["time"] = out.time.str[len(var):].astype(int)
    out["baseline"] = out.study_id.map(df.set_index("study_id")[f"{var}1"])
    return out


section("Haemodynamics (series failed integrity screening, shown for audit only)")
mixed = []
for var in ["sbp", "dbp", "hr"]:
    data = long_form(var)
    for label, formula, frame in [("group + time", "value ~ treated + C(time)", data),
                                  ("baseline adjusted, T2 to T20", "value ~ treated + baseline + C(time)", data[data.time > 1])]:
        fit = smf.mixedlm(formula, frame, groups=frame.study_id).fit(reml=True)
        lo, hi = fit.conf_int().loc["treated"]
        re_var = float(fit.cov_re.iloc[0, 0])
        mixed.append({"variable": VITALS[var], "model": label, "difference": round(fit.params["treated"], 2), "low": round(lo, 2),
                      "high": round(hi, 2), "P": f"{fit.pvalues['treated']:.3g}",
                      "ICC": round(re_var / (re_var + fit.scale), 2), "status": STATUS})
    ordered = data[data.time > 1].sort_values(["study_id", "time"])
    structure = sm.cov_struct.Autoregressive(grid=True)
    gee = smf.gee("value ~ treated + baseline + C(time)", "study_id", ordered, time=ordered.time, cov_struct=structure).fit()
    lo, hi = gee.conf_int().loc["treated"]
    mixed.append({"variable": VITALS[var], "model": "GEE AR(1), baseline adjusted, T2 to T20", "difference": round(gee.params["treated"], 2),
                  "low": round(lo, 2), "high": round(hi, 2), "P": f"{gee.pvalues['treated']:.3g}",
                  "ICC": np.nan, "status": STATUS})
print(keep("repeated_measures", pd.DataFrame(mixed)).to_string(index=False))

sbp = df[cols("sbp")]
events = {
    "SBP fall of more than 20% from baseline": sbp.iloc[:, 1:].lt(0.8 * df.sbp1, axis=0).any(axis=1),
    "SBP below 90 mmHg": sbp.lt(90).any(axis=1),
    "Heart rate below 60": df[cols("hr")].lt(60).any(axis=1),
    "SpO2 below 95%": df[cols("spo2")].lt(95).any(axis=1),
}
rows = []
for label, mask in events.items():
    k = mask.groupby(df.arm).sum()
    n = df.groupby("arm").size()
    p = stats.fisher_exact([[k[a], n[a] - k[a]] for a in COLORS]).pvalue
    rows.append({"event": label, **{a: f"{int(k[a])}/{n[a]}" for a in COLORS}, "P Fisher": round(p, 3)})
print(keep("haemodynamic_events", pd.DataFrame(rows).assign(status=STATUS)).to_string(index=False))

removed = long_form("sbp").merge(errors.assign(time=errors.variable.str[3:].astype(int))[["study_id", "time"]], how="left", indicator=True)
sensitivity = []
for label, frame in [("four-digit entries corrected", removed), ("four-digit entries set to missing", removed[removed["_merge"] == "left_only"])]:
    fit = smf.mixedlm("value ~ treated + C(time)", frame, groups=frame.study_id).fit(reml=True)
    low, high = fit.conf_int().loc["treated"]
    sensitivity.append({"systolic data": label, "difference": round(fit.params["treated"], 2), "low": round(low, 2), "high": round(high, 2), "status": STATUS})
print(keep("correction_sensitivity", pd.DataFrame(sensitivity)).to_string(index=False))

pointwise = pd.DataFrame({f"{v} difference": df[df.treated == 1][cols(v)].mean().to_numpy()
                          - df[df.treated == 0][cols(v)].mean().to_numpy() for v in ["sbp", "dbp", "hr"]},
                         index=[f"T{t}" for t in range(1, CFG.n_time + 1)])
print("\nmean difference, paracetamol minus control, at each time point (descriptive)")
print(keep("pointwise_differences", pointwise.round(1).assign(status=STATUS)).to_string())

#%%
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
for ax, var in zip(axes, ["sbp", "dbp", "hr"]):
    for arm, colour in COLORS.items():
        values = df.loc[df.arm == arm, cols(var)]
        mean, half = values.mean(), 1.96 * values.sem()
        t = np.arange(1, CFG.n_time + 1)
        ax.plot(t, mean, color=colour, marker="o", ms=3, label=arm)
        ax.fill_between(t, mean - half, mean + half, color=colour, alpha=0.2)
    ax.set(xlabel="Time point", ylabel=VITALS[var])
    ax.legend(frameon=False)
fig.suptitle("UNVERIFIED haemodynamic series. Not for publication.", color="red")
save(fig, "fig_haemodynamics_audit")

#%%
section("Design evaluation")
p_control, p_treated = 0.35, 0.20
pbar = (p_control + p_treated) / 2
za, zb = stats.norm.ppf(0.975), stats.norm.ppf(0.80)
n_group = (za * np.sqrt(2 * pbar * (1 - pbar)) + zb * np.sqrt(p_control * (1 - p_control) + p_treated * (1 - p_treated))) ** 2 / (p_control - p_treated) ** 2
h = 2 * np.arcsin(np.sqrt(p_control)) - 2 * np.arcsin(np.sqrt(p_treated))
arcsine = NormalIndPower().solve_power(effect_size=h, alpha=0.05, power=0.80)
design = pd.Series({"n per group, pooled-variance formula": np.ceil(n_group), "n per group, arcsine method": np.ceil(arcsine),
                    "with 8% attrition": np.ceil(np.ceil(n_group) / 0.92), "analysed with outcome": int(df.shiver.notna().sum())})
print(keep("sample_size", design.rename("value").to_frame()).to_string())

with pd.ExcelWriter(OUT / "analysis_results.xlsx") as writer:
    for name, table in TABLES.items():
        table.to_excel(writer, sheet_name=name[:31], index=name in {"completeness", "eligibility_counts", "terminal_digits", "wrench_grades", "bayes_fragility", "sample_size", "pointwise_differences"})
df.drop(columns=["diagnosis", "excel_row", "no"]).to_csv(OUT / "analysis_dataset_pseudonymised.csv", index=False)
section("Run complete")
print(f"integrity gate passed: {DATA_VERIFIED}")
print(sorted(p.name for p in OUT.iterdir()))
