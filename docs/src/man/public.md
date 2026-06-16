# Public elements

The functions and types below are exported by `MassJ`. Full signatures and
docstrings are on the [References](@ref) page; the manual pages linked from each
group give worked examples.

**I/O and selection**
- [`info`](@ref)
- [`load`](@ref)
- [`load_usi`](@ref)
- [`save`](@ref)
- [`average`](@ref)
- [`extract`](@ref)
- [`chromatogram`](@ref)
- [`mobilogram`](@ref)
- [`ionogram`](@ref)

**Processing**
- [`smooth`](@ref)
- [`baseline_correction`](@ref)
- [`centroid`](@ref)
- [`chrom_peaks`](@ref) / [`chrom_peak`](@ref)
- [`detect_features`](@ref) / [`feature_table`](@ref)

**Scientific calculations**
- [`formula`](@ref) / [`masses`](@ref)
- [`isotopic_distribution`](@ref) / [`isotope_table`](@ref) / [`simulate`](@ref)
- [`deconv`](@ref)
- [`adduct_mz`](@ref) / [`neutral_mass`](@ref)
- [`calibrate`](@ref)
- [`assign_formula`](@ref) / [`score_isotope_pattern`](@ref)

**Energy-resolved yields**
- [`yields`](@ref) / [`integrate_window`](@ref) / [`read_peaklist`](@ref)
- [`normalize_tic`](@ref) / [`normalize_external`](@ref)
- [`combine_yields`](@ref) / [`shift_x`](@ref) / [`scale_yields`](@ref) / [`recalibrate_x`](@ref)
- [`trim_yields`](@ref) / [`restrict_x`](@ref) / [`drop_peaks`](@ref)

**Chimeric spectra**
- [`abundance_matrix`](@ref) / [`partial_correlation`](@ref) / [`cmi_matrix`](@ref)
- [`cluster_ions`](@ref) / [`cluster_spectra`](@ref)

**Ecosystem interoperability**
- [`featurize`](@ref) / [`select_spectra`](@ref) / [`from_matrix`](@ref) / [`spectra_table`](@ref)

**Uncertainties**
- [`nscans`](@ref) / [`stdev`](@ref) / [`sem`](@ref)
- [`measurements`](@ref) / [`withunits`](@ref)

**Exported types**
- [`MassJ.MSscans`](@ref) / [`MassJ.IonCurrent`](@ref) / [`MassJ.MSrun`](@ref)
- [`MassJ.AbstractPeak`](@ref) / [`MassJ.Peak`](@ref) / [`MassJ.TargetPeak`](@ref)
- [`MassJ.YieldCurve`](@ref) / [`MassJ.ChromPeak`](@ref) / [`MassJ.Feature`](@ref)
- [`MassJ.Adduct`](@ref) / [`MassJ.Calibration`](@ref)
- [`MassJ.FormulaCandidate`](@ref)
