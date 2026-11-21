# STROBE / Record Checklist for the German Varicella Implementation

Mapped to the active illustrative implementation (`experiment/code/application.py`, 400 German counties in $G=16$ federal states, 2015--2019, $T=261$ estimation weeks; replacing the earlier archived US mask-order prototype):

| Item | Location in Manuscript / Supporting Information |
|---|---|
| Title & abstract | `sections/abstract.tex` (states count data, sparse surveillance panel, illustrative scope) |
| Background & rationale | `sections/introduction.tex` ($\S 1$), `sections/application.tex` ($\S 4$ opening) |
| Study design & objectives | `sections/application.tex` ($\S 4.1$: computational and diagnostic sparse-count testbed) |
| Setting, dates & eligibility | `sections/application.tex` ($\S 4.1$), `experiment/data/data_contract.md`, Web Appendix~D |
| Data sources & SHA-256 provenance | `sections/backmatter.tex`, `experiment/data/fetch_germany.py`, `experiment/data/checksums_germany.sha256` |
| Variables & exposure coding | `sections/application.tex` ($\S 4.1$: one-week-lagged state summer holiday $\ge 4$/7 days, population $\omega_i$, adjacency $W$; sensitivity grids Table~7b) |
| Bias, reporting & diagnostic checks | `sections/application.tex` ($\S 4.3$--$\S 4.4$), Table~3 (`tab:envelope`), Web Appendices~D--F (Tables~11--18) |
| Quantitative variables & zero handling | `sections/application.tex` ($\S 4.1$--$\S 4.2$, 58.8\% county-week zeros retained without deletion; NUTS-2 aggregation) |
| Statistical methods & few-cluster inference | `sections/methods.tex` ($\S 2$), `sections/application.tex` ($\S 4.3$: $t_{15}$, WCR-score Webb6/Rademacher, jackknife, pairs bootstrap, Stock--Wright $S_1$/$S_3$, A4/A4') |
| Descriptive & main working results | `sections/application.tex` ($\S 4.1$--$\S 4.2$), `tables/table_7_application.tex`, `figures/fig_03_application.pdf` |
| Limitations & scope statement | `sections/application.tex` ($\S 4.4$ Scope statement: notification-week reporting, Assumption~A7 failure, non-absorbing seasonal holidays, non-causal framing) |
| Funding & reproducibility | `sections/backmatter.tex`, `scripts/rebuild_manifest.py`, `tables/manifest.json` |
