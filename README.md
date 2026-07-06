# codex-test

A MATLAB test repository for qualitatively reproducing a Fig. 4-style combined-resonance frequency response.

## Recommended workflow

Run these commands from the repository root in MATLAB:

```matlab
No4_MCK_w_fig4
No6_RK_fig4
```

`No4_MCK_w_fig4` keeps the Galerkin derivation framework and writes `No4_coeff_fig4.mat`. For the current Scheme A workflow, the saved `coef` uses a dimensionless tuning model, while `raw` keeps the symbolic Galerkin coefficients for reference.

`No6_RK_fig4` loads `coef`, draws one paper-style smooth combined-resonance response curve, prints curve diagnostics, and saves the MATLAB figure file `Fig4_qualitative_response.fig`.
