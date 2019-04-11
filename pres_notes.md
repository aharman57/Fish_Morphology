## presentation notes

### bmb

- the pairs plots are nice (you might consider the 'enhanced' versions that are shown in Ian Dworkin's notes: see `car::scatterplotMatrix`). For correlation plots, I like `corrplot::corrplot.mixed(..., lower="ellipse", upper="number")`.
- there seems to be some heteroscedasticity left in your data: try Box-Cox transformation?  Or use a permutation-test approach (as indicated in Ian's notes: `vegan::adonis` and functions in the `geomorph` package) for p-value calculations?
- you might want to think about *Tobit regression* (e.g. see the `AER` package). The Tobit model assumes that the data are normally distributed, but that values less than zero are truncated (https://en.wikipedia.org/wiki/Tobit_model) ...

---

Need to think some more about allometries.  What we have before essentially assumes isometry with respect to other traits, but the trait ratios vary with environmental values.  If we want the normal kind of allometric relationships we might need to add log(total length) as a predictor in the regression ...

$$
\begin{split}
y_1 & = \beta_{0,1} + \beta_{1,1} x_1 + \beta_{2,1} x_2 + \\
y_2 & = \beta_{0,2} + \beta_{1,2} x_1 + \beta_{2,2} x_2 + \\
\end{split}
$$
here the second index indicates which morphometric variable we're looking at (i.e. \beta_{i,j} is the effect of the i'th predictor variable on the j'th response variable). Let's say that all of the $y_j$ values are the logs of the raw morphological variables, and let's say we take the first variable (e.g., total length) as the reference trait. Now we define $a_j = y_j-y_1$ as the allometry of trait $j$ with respect to length.  So

$$
a_j & = (\beta_{0,j}-\beta_{0,1}) + (\beta_{1,j}-\beta_{1,1}) x_1 + ...
$$

so this implies a *constant* ratio between traits
