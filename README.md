# The Coors Effect

What Coors Field does to pitch movement, and which pitchers keep their value there.

**Read the report:** https://charliewiebe.github.io/coors-field/analysis.html

**Portfolio:** https://charliewiebe.github.io

Using 2.8 million public Statcast pitches (2023–26), this project compares the same pitcher, pitch and season at Coors Field and everywhere else. It splits each pitch's acceleration into drag, Magnus lift and seam-shifted wake, rebuilds every pitch as if thrown at Coors, and grades both versions with the [Pitch Shaper](https://github.com/charliewiebe/pitch-shaper) grader.

## Files

```
analysis.html             the report (rendered from coors_field.Rmd)
coors_field.Rmd           report source; reads coors_out/
04_coors_movement.R       full analysis: physics, paired comparisons, bootstrap, Coors simulation, grading
coors_out/                outputs of 04_coors_movement.R used by the report
```

`04_coors_movement.R` runs after the Pitch Shaper pipeline (`01`–`03`), which pulls the Statcast data and builds the grader. Public data only.
