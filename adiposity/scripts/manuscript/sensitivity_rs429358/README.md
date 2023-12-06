Scripts in this folder:

Following reviewer request, check the rs429358 associations for survival bias more closely.

1. **association_with_followup.R** - Linear regression for association of rs429358 allele dosage with number of follow-up measures of adiposity and length of follow-up.
2. **association_with_slope_no_age_adjustment.R** - Compare the association of rs429358 allele with b1 BLUP from LMEs, unadjusted for age and age-squared.
3. **plot_effect_sizes.R** - Plot results from models above
4. **young_vs_old.R** - Double-check that early mortality of 'C' allele carriers cannot be associated with weight-lowering by 'C' allele. To do this, plot weight-change slopes separately in young ages (30-50yrs) and old ages (60-80yrs) -- slopes are lower in older people so early mortality of 'C' carriers cannot be responsible for this association. 
