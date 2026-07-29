# Changes in version 1.1.9

* Added the canonical `computeQScore()`, `computeQScoreFlags()`,
  `computeOutliersQScore()`, and `applyQScoreModel()` APIs.
* Retained `computeQCScore()`, `computeQCScoreFlags()`,
  `computeOutliersQCScore()`, and `applyQCScoreModel()` as deprecated
  compatibility APIs. They preserve historical arguments, `QC_score`,
  `low_qcscore`, `low_threshold_qcscore`, and `QCScore_model*` metadata.
* Restored historical plotting arguments (`size`, `alpha`, `alphaNumbers`,
  `csize`, `calpha`, `mapNumbersCol`, and `mapAlphaNumbers`) as aliases for
  the canonical plotting arguments. Canonical arguments take precedence when
  both forms are supplied.
* Restored `getModelFormula(formulaVars, verbose=FALSE)` compatibility while
  accepting `metricList` as a named replacement.
* Restored the public `computeLambda(trainDF, modelFormula)` signature,
  complete-case filtering, and support for both `qcscore_train` and
  `QScore_train`.
* User-supplied `modelFormula` values are now fitted without adding, removing,
  or rebuilding terms. Supported subsets, additive formulas, and selected
  interactions are preserved.
* Formula predictors are limited to `log2SignalDensity`, `Area_um`,
  `log2AspectRatio`, and `log2Ctrl_total_ratio`. Unsupported predictors and
  transformations now produce informative errors, and CosMx-only border terms
  are checked against dataset technology.
* Model transfer now has a canonical public API, preserves training matrix
  column order, supports validated custom formulas, and retains the historical
  deprecated interface.
* Fixed the SpaceTrooper utilities vignette name and added the bioRxiv
  citation.
* Follow-up: intercept handling remains unchanged in this release and should be
  reviewed consistently across training, lambda selection, prediction, and
  model transfer.

# Changes in version 1.1.7

* fixing author name typo

# Changes in version 1.1.6

* adding new vignette about SpaceTrooper utilities

# Changes in version 1.1.5

* minor enhancement of cosmx input reading, polygonsCol added
* minor documentation clarification on QCScore computation
* cleaned LICENSE stub and creating LICENSE.md file

# Changes in version 1.1.4

* updating documentation along multiple functions
* adding graphical abstract
* updating README
* updating Vignettes and providing better data examples

# Changes in version 1.1.3

* fixing verbose message in QC functions.

# Changes in version 1.1.2

* adding volume instead of area for MERFISH.
* refactoring log2CountArea to log2DensitySignal
* added size and alpha arguments to plotCellsFovs
* added scaleBar argument to plotCellsFovs, plotCentroids, plotPolygons and plotZoomFovsMap

# Changes in version 1.0.1

* fixing minor bugs in AspectRatio internal computation for technology missing it.
* fixing merfish reading where colData where not properly sorted and sync with assay cells.

# Changes in version 0.99.0

* adding unit tests on QC steps
* adding datasets descriptions
