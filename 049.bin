# EFA, step by step, for your SEC-SAXS data

EFA here means **evolving factor analysis**. The input is the entire family of scattering curves across SEC frames, not just the integrated SEC trace or Rg time series. Integrating over q would discard the shape information needed for separation.

## 1. Arrange the measured curves into a matrix

Let i index q and j index the SEC frame:

\[
D_{ij}=I(q_i,t_j),\qquad D\in\mathbb{R}^{N_q\times N_t}.
\]

A column is one measured SAXS curve. A row follows a single q value through time. After background correction, assume independent dilute particles and K components whose ensemble-averaged scattering shapes remain fixed during elution:

\[
D_{ij}=\sum_{k=1}^{K}S_{ik}C_{jk}+E_{ij},\qquad D=SC^T+E.
\]

S has shape Nq×K and contains the component curves. C has shape Nt×K and contains their elution amplitudes. E includes measurement noise and deviations from the assumptions. A single IDP can have many conformations yet be one scattering component if its ensemble-average curve is constant over time. A changing ensemble or concentration-dependent interaction can add a mathematical component without adding an oligomeric species.

The scale is arbitrary: multiplying one S column by a and dividing the corresponding C column by a gives exactly the same data. Here each C column sums to one, so C is a relative elution coefficient, not a molar concentration or mole fraction.

## 2. Make q points comparable, without destroying rank

Low-q intensities and errors are generally larger than high-q values. Following RAW's approach, calculate the mean error for each q row:

\[
\bar\sigma_i=\frac1{N_t}\sum_j\sigma_{ij},\qquad A_{ij}=D_{ij}/\bar\sigma_i.
\]

This is multiplication by a diagonal matrix W on the left: A=WD. It preserves the noiseless rank and elution factorization, A=(WS)Cᵀ. Element-by-element division by time-varying errors need not preserve a rank-K factorization. The code therefore uses q-row weights for decomposition and the individual propagated errors for residual inspection. We do not subtract the temporal mean, which would remove part of the physical scattering signal.

Residual-background correction adds correlated uncertainty. For a mean background B from n independent frames:

\[
D_j=Y_j-B,\quad \operatorname{Var}(D_j)=\sigma_j^2+\operatorname{Var}(B),
\quad \operatorname{Cov}(D_j,D_l)=\operatorname{Var}(B)\quad(j\ne l),
\]

when the sample frames do not overlap the buffer interval. The Monte Carlo procedure models this additional shared correction. Covariance inherited from earlier beamline subtraction remains unknown.

## 3. Use SVD to measure independent patterns

\[
A=U\Sigma V^T=\sum_r s_r\,u_r v_r^T.
\]

Each term is an outer product: one pattern across q times one pattern across frames. Singular values s1≥s2≥… measure the strength of those patterns in the weighted data. If one fixed species simply changes concentration, every column is proportional and noiseless rank is one. Two distinguishable curves with independently changing proportions can give rank two.

A singular vector is a mathematical basis vector. It can contain negative values and is not itself a purified species curve. Noise makes almost every measured matrix formally full rank, so the question is which modes exceed a plausible noise model and show coherent q/time structure.

Our null calculation fits a rank-one mean, adds simulated measurement noise, repeats the background correction, and asks how large s2 becomes. The clean peaks' s2 values are not unusual using supplied errors; smaller empirically estimated errors make weak deviations detectable. This is why the result is not proof of monodispersity.

## 4. Make the factor analysis evolve through the SEC frames

At each frame j, compute:

\[
A^{\mathrm{forward}}_j=A[:,1:j],\qquad
A^{\mathrm{backward}}_j=A[:,j:N_t].
\]

Run SVD on each growing submatrix and plot its singular values against j. In the forward scan, a new independent signal becomes detectable as a component enters. In the backward scan, adding earlier frames reveals components that are absent from later portions. Jointly the scans suggest component support intervals [ak,bk].

Noise singular values also rise as more frames are included. EFA curves are not concentration curves, and their crossings do not supply exact molecular start/stop times. In WT, the second mode becomes prominent as the later main peak joins the earlier shoulder. The broad dataset therefore contains two distinguishable patterns even though the clean peak alone lacks robust separation constraints.

## 5. Rotate the SVD basis toward physical components

The truncated approximation is A≈UKΣKVKᵀ. For any invertible K×K matrix T:

\[
C=V_K T,\qquad S_w=U_K\Sigma_KT^{-T},\qquad A\approx S_wC^T.
\]

There are infinitely many such rotations. EFA reduces that ambiguity by imposing time support and usually nonnegative elution amplitudes:

\[
C_{jk}=0\ \text{outside }[a_k,b_k],\qquad C_{jk}\ge0.
\]

With M as the binary support mask and \(\widetilde C=M\odot C\), alternating least-squares updates take the form:

\[
S_w=A(\widetilde C^T)^+,\qquad C_{\mathrm{new}}=(S_w^+A)^T,
\]

followed by the concentration constraints and normalization. The superscript + denotes the Moore–Penrose pseudoinverse. RAW's hybrid method initializes from an explicit rotation and then iterates. The original-unit curve is reconstructed by S=D(Ĉᵀ)⁺.

Convergence only means the iteration stopped changing. It does not establish unique components or physical validity. Check sensitivity to time windows, negative curves, baseline, q range, and residuals:

\[
R_{ij}=\frac{D_{ij}-\sum_k S_{ik}C_{jk}}{\sigma_{ij}}.
\]

The report lists mean R² as a diagnostic. It is not a formally reduced chi-square because fitted degrees of freedom and correlations are not fully accounted for.

Our WT shoulder rotation converges, but its early recovered curve and residuals fail these stronger checks. We therefore keep it as exploratory evidence of heterogeneity.

## 6. Estimate Rg from a defensible recovered curve

In the Guinier region:

\[
I(q)\approx I(0)e^{-q^2R_g^2/3},\qquad
\ln I(q)=\ln I(0)-\frac{R_g^2}{3}q^2.
\]

If the slope of ln I against q² is m, then Rg=√(−3m). We fit conservatively with qmax Rg≤1 and examine qmin sensitivity. The dominant clean profiles give about 25.6 Å for WT and 32.4 Å for AA under the primary correction. Those values are apparent sizes of the dominant scattering patterns, not independently established monomer sizes.

For an unresolved dilute mixture, the small-q expansion gives:

\[
R_{g,\mathrm{app}}^2=\frac{\sum_k I_k(0)R_{g,k}^2}{\sum_k I_k(0)}.
\]

Thus the measured Rg is weighted by forward scattering, not by molecule counts. Larger particles can contribute strongly even at small molar fractions.

## 7. Connect components to oligomerization

EFA asks whether distinct scattering patterns can be separated. Oligomer assignment asks how many protein chains each pattern represents. Those are different inference steps.

Your nominal WT monomer mass is 12.1 kDa; a dimer is about 24.2 kDa. Under the same contrast and calibration, I(0) is proportional to number concentration times molecular mass squared, or mass concentration times molecular mass. Without the relevant concentration and intensity calibration, EFA amplitudes cannot supply that molecular mass.

Rg is also sensitive to shape and disorder. An expanded monomer can overlap a compact oligomer in size. Even in an idealized dimer of unchanged identical subunits, Rg,dimer²=Rg,monomer²+d²/4, where d is the separation of subunit centers. There is no universal dimer-to-monomer Rg multiplier for an IDP.

The current evidence supports a dominant clean-peak pattern and additional WT shoulder heterogeneity. It does not resolve an additional IDP-sized oligomer in the clean regions, exclude co-eluting oligomers, or quantify an oligomer fraction.

References: [RAW EFA tutorial](https://bioxtas-raw.readthedocs.io/en/latest/tutorial/s2_efa.html), [Meisburger et al. 2016](https://doi.org/10.1021/jacs.6b01563), and the archived numerical source/provenance supplied with this analysis.
