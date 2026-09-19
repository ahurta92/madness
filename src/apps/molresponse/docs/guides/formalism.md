# Response formalism

The shared machinery every property in this engine specializes. The individual
feature guides ([polarizability](polarizability.md),
[hyperpolarizability](hyperpolarizability.md), [Raman](raman.md),
[two-photon absorption](two_photon_absorption.md)) point back here rather than
re-deriving it. The equations are those of the release report 2026-09-09
(`madness-workspace/reports/2026-09-09_release_report/main.tex`, §1); the same
equations are stated in the Doxygen at the object that implements each one
(`formalism.dox`, page *Response formalism*).

## Objects and conventions

Occupied orbitals $\phi_i$ ($i=1..N_{\rm occ}$, closed shell), ground-state Fock operator
$F^{0}$ with $F^{0}\phi_i=\epsilon_i\phi_i$, projector onto the virtual space
$\hat Q=1-\sum_i|\phi_i\rangle\langle\phi_i|$. A pair density $\gamma(r,r')$ has the diagonal
$\rho_\gamma(r)=\gamma(r,r)$ (carrying the closed-shell factor 2) and enters through the first
derivative of the two-electron operator,

$$
g'[\gamma]f = J[\rho_\gamma]f - c_x K[\gamma]f,\qquad
(K[\gamma]f)(r) = \int \frac{\gamma(r,r')\,f(r')}{|r-r'|}\,dr',\qquad
g'[\gamma]^{\dagger} = g'[\gamma^{\dagger}],
$$

with $c_x=1$ for Hartree–Fock. Every convention that matters is fixed by one statement:
*a response density at $+\omega$ is $\gamma^{B}=|x^{B}\rangle\langle\phi|+|\phi\rangle\langle y^{B}|$,
and its dagger is the same density at $-\omega$.* In the code: `gamma_legs`,
`gamma_dagger_legs`, `ketbra` in `kernels/source_spec.hpp`.

## First order: the linear response in MRA form

For a one-electron perturbation $v^{B}$ at frequency $+\omega_B$ the response is the pair
$(x^{B},y^{B})$, both $\hat Q$-projected, with

$$
\gamma^{B}=\sum_i\bigl[|x_i^{B}\rangle\langle\phi_i|+|\phi_i\rangle\langle y_i^{B}|\bigr],\qquad
F^{B}=v^{B}+g'[\gamma^{B}],\qquad \bar F^{B}=(F^{B})^{\dagger}=v^{B}+g'[\gamma^{B\dagger}],
$$

and the coupled equations

$$
(F^{0}-\epsilon_i-\omega_B)\,x_i^{B}=-\hat QF^{B}\phi_i,\qquad
(F^{0}-\epsilon_i+\omega_B)\,y_i^{B}=-\hat Q\bar F^{B}\phi_i .
$$

MADNESS never forms $F^{0}$ as a matrix. Writing $F^{0}=-\tfrac12\nabla^2+V^{0}$ and moving the
potential to the right, each equation is inverted with the bound-state Helmholtz (BSH) Green's
function,

$$
x_i^{B}=-2\,\hat G_{\mu_i^{-}}\Bigl[V^{0}x_i^{B}+\hat QF^{B}\phi_i\Bigr],\qquad
\hat G_{\mu}=(-\nabla^2+\mu^2)^{-1},\qquad
\mu_i^{\mp}=\sqrt{-2(\epsilon_i\pm\omega_B)},
$$

(and the same for $y$ with $\mu_i^{+}$), iterated to self-consistency with KAIN acceleration;
convergence is judged on the BSH residual and on the change of $\rho_{\gamma^B}$, at each rung
of the threshold ladder (protocol $10^{-4}\to10^{-6}$, wavelet order $k$ following the rung).
This linear solve (`solve_fd_protocol`) is the single most reused component in the engine.

At $\omega_B=0$, $x^{B}=y^{B}$ and only one function per orbital is solved (`Static`); at finite
frequency both channels are kept (`Full`). The polarizability is the trace
$\alpha_{AB}(\omega)=-\mathrm{Tr}(v^{A}\gamma^{B})=-\sum_i\bigl[\langle\phi_i|v^{A}|x_i^{B}\rangle+\langle y_i^{B}|v^{A}|\phi_i\rangle\bigr]$.

**Excited states** are the same operator with $v\equiv0$: $(x^{f},y^{f})$ at $\omega_f$ solve the
equations above with right-hand sides $-\hat Qg'[\gamma^{f}]\phi_i$ and
$-\hat Qg'[\gamma^{f\dagger}]\phi_i$ (RPA / full TDHF; dropping the $y$ channel gives TDA).
Eigenvectors are normalized as $\langle x^f|x^f\rangle-\langle y^f|y^f\rangle=1$.

## Second order: the one source

With both photon frequencies positive, the Fock operator that meets $\gamma^{C}$ is $F^{B}$ at
$+\omega_B$, *undaggered*. Idempotency fixes the occupied–occupied and virtual–virtual blocks of
the second-order density without any solve,

$$
\gamma_L^{BC}=\sum_i\Bigl[|x_i^{B}\rangle\langle y_i^{C}|+|x_i^{C}\rangle\langle y_i^{B}|
  -|\phi_i\rangle\langle\zeta_i^{BC}|-|\phi_i\rangle\langle\zeta_i^{CB}|\Bigr],\qquad
\zeta_i^{BC}=\sum_j\phi_j\,\langle y_i^{B}|x_j^{C}\rangle ,
$$

and the $e^{-i\omega_\sigma t}$ component of the equation of motion, $\hat Q$-projected, is the
second-order linear-response equation with the source

$$
\begin{aligned}
P_p^{BC}&=(1+\mathcal P^{BC})\Bigl[\underbrace{\textstyle\sum_kx_k^{C}F^{B}_{kp}}_{[M]}
   \;\underbrace{-\;\hat QF^{B}x_p^{C}}_{[A]}\Bigr]
   \;\underbrace{-\;g'[\gamma_L^{BC}]\phi_p}_{[L]}
   \;\underbrace{-\;g''[\gamma^{B}\gamma^{C}+\gamma^{C}\gamma^{B}]\phi_p}_{[G],\ \mathrm{HF}:\,0},
   \qquad F^{B}_{kp}=\langle\phi_k|F^{B}|\phi_p\rangle,\\
Q_p^{BC}&=P_p^{BC}\big|_{x\leftrightarrow y\ \text{on every leg}}
   =(1+\mathcal P^{BC})\Bigl[\textstyle\sum_ky_k^{C}\bar F^{B}_{kp}-\hat Q\bar F^{B}y_p^{C}\Bigr]
   -g'[\gamma_L^{BC\dagger}]\phi_p-\dots
\end{aligned}
$$

where $\mathcal P^{BC}$ swaps $B\leftrightarrow C$. This is the *only* second-order source
(`quadratic_source` in `kernels/tpa_source_spec.hpp`; `compute_vbc_spec` in `kernels/vbc.hpp` is
the same equation term by term): β, Raman and 2PA all contract it; they differ only in what they
contract it with.

## β: the 2n+1 contraction

$$
\beta_{ABC}=-2\,(b_1+b_2+b_3),\qquad
b_1=-\bigl[\langle x^{A}|P^{BC}\rangle+\langle y^{A}|Q^{BC}\rangle\bigr],\qquad
b_2=\langle v^{A}y^{C}|x^{B}\rangle+\langle v^{A}\zeta^{BC}|\phi\rangle,\qquad
b_3=b_2|_{B\leftrightarrow C}
$$

with $(x^{A},y^{A})$ the driven response at $+\omega_\sigma=\omega_B+\omega_C$; no second-order
vector is ever solved for. SHG is $\omega_B=\omega_C=\omega$; in the static limit $x=y$, $P=Q$
and Kleinman symmetry holds exactly. See [hyperpolarizability](hyperpolarizability.md).

## Raman: the nuclear-displacement leg

$$
\frac{\partial\alpha_{AB}(\omega)}{\partial Q}=\beta_{ABQ}(-\omega;\omega,0),\qquad
v^{Q}=\frac{\partial V_{\rm nuc}}{\partial Q}=\sum_{\alpha}Z_\alpha\frac{\partial}{\partial Q}\frac{-1}{|r-R_\alpha|},
$$

so the β machinery is reused with $C\to Q$: the $C$ leg is the static response to $v^{Q}$
(`nuclear_operator`, `MolecularDerivativeFunctor`), the $B$ leg the dipole response at $\omega$,
the $A$ leg the dipole response at $\omega_\sigma=\omega$. See [Raman](raman.md).

## Two-photon absorption: the residue of the quadratic response

Write the coupled first-order equations as $(\Lambda-\omega\Delta)|X,Y\rangle=-|P,Q\rangle$ with
$\Lambda$ the electronic Hessian and $\Delta=\mathrm{diag}(1,-1)$. Completeness of the paired
eigenvectors in the $\Delta$ metric gives

$$
(\Lambda-\omega\Delta)^{-1}
= \sum_K\left[\frac{|X^K,Y^K\rangle\langle X^K,Y^K|}{\Omega_K-\omega}
 - \frac{|Y^K,X^K\rangle\langle Y^K,X^K|}{\Omega_K+\omega}\right],
$$

and the residue of $\beta_{ABC}(-\omega';\omega_B,\omega'-\omega_B)$ at $\omega'\to\Omega_N$
identifies the two-photon transition moment as the eigenvector contracted with the *same* source
$(P,Q)$ — no dagger, nothing changed inside $(P,Q)$. With degenerate photons
$\omega_B=\omega_C=\omega_f/2$,

$$
S_{BC}=\sqrt2\,\bigl[\langle x^{f}|P^{BC}\rangle+\langle y^{f}|Q^{BC}\rangle\bigr],\qquad
\delta^{\rm 2PA}=\tfrac{1}{30}\textstyle\sum_{bc}\bigl[F\,S_{bb}S_{cc}+(G{+}H)S_{bc}S_{bc}\bigr],
$$

with $F=G=H=2$ for linearly polarized parallel photons and $\sqrt2$ the translation between
this solver's eigenvector normalization and DALTON's (`tpa_moment_residue` in `kernels/tpa.hpp`).
See [two-photon absorption](two_photon_absorption.md).

## References

- Release report 2026-09-09, §1 "Formalism" (`madness-workspace/reports/2026-09-09_release_report`).
- First-principles derivation and orientation derivation
  (`reports/2026-09-08_first_principles_derivation`, `reports/2026-09-09_orientation_derivation`);
  working-equation form in `reports/2026-09-09_beta_tpa_working_equations`.
- Parker *et al.* (2018) and Sałek *et al.* (2002), mirrored line by line in the derivations above.
