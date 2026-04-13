[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]()

# Vibrational Spectroscopy
__What is it?__: Vibrational spectroscopy detects transitions between the quantised vibrational energy levels associated with bond stretching and/or bond angle bending in molecules.

__How do we do it?__: Transitions are observed by measuring the amount of infrared radiation that is absorbed or emitted by vibrating molecules in solid, liquid, or gas phases.

__Why do we do it?__: A knowledge of the vibrational level spacings gives us the value of the stretching (or bending) force constants which characterise the stiffness of a bond, allows us to estimate the bond dissociation energy, and guves us a means of identfying characteristic functional groups of atoms within large molecules.

# Quantum approach
In quantum mechanics and in other branches of physics, it is common to approach physical problems algebraics and analytic methods. Examples include the use of differential equations for many interesting models, the use of quantum groups in quantum physics, and of differential geometry in relativity theory. In this work, we discuss the Hermite polynomials, some of their properties and a brief description of their applications to the Quantum Harmonic Oscillator.

The Harmonic Oscillator's Quantum Mechanical solution involves Hermite Polynomials, which are introduced here. The wavefunctions for the quantum harmonic oscillator contain the Gaussian form, which allows them to satisfy the necessary boundary conditions at infinity. In the wavefunction associated with a given value of the quantum number $n$, the Gaussian is multiplied by a polynomial of order $n$ (the Hermite polynomials above) and the constants necessary to normalize the wavefunctions.

# Hermite Polynomials
Hermite polynomials, named after the French mathematician [Charles Hermite](https://en.wikipedia.org/wiki/Charles_Hermite), are orthogonal polynomials, in a sense to be described below, of the form

<p align='center'>
    $$H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2} \quad (1)$$
</p>

for $n=0, \ 1, \ 2, \ 3, \ldots$
The first few Hermite polynomials are
- for $n=0$ we have $H_0 (x) = 1$
- for $n=1$ we have $H_1 (x) = 2x$
- for $n=2$ we have $H_2 (x) = 4x^2 - 2$
- for $n=3$ we have $H_3 (x) = 8x^3 - 12x$
- for $n=4$ we have $H_4 (x) = 16x^4 - 48x^2 + 12$
- for $n=5$ we have $H_5 (x) = 32x^5 - 160x^3 + 120x$
- for $n=6$ we have $H_6 (x) = 64x^6 - 480x^4 + 720x^2 -120$
- for $n=7$ we have $H_7 (x) = 128x^7 - 1344x^5 + 3360x^3 - 1680x$

For $n \in \mathbb{N}$, we define Hermite polynomials $H_n (x)$ by
<p align='center'>
    $$\sum_{n=0}^{\infty} \frac{H_n (x)}{n!} r^n = e^{2xr - r^2}, \quad for |r| < \infty \quad (2) $$
</p>

To find $H_n (x)$, expand the right-hand side of eq.(2) as a Maclaurin series in $r$ and equate coefficients. From eq.(2) we drive the closed expression
<p align='center'>
    $$H_n (x) = \sum_{k=0}^{[n/2]} \frac{(-1)^k n!}{k! (n-2k)!} (2x)^{n-2k} \quad (3) $$
</p>

where $[x]$ denotes the largest integer less than or equal to $x$. Checking with $n = 0, \ 1, \ 2, \ldots$ we find that eq.(3) yields the expected Hermite polynomials. To prove that eq.(3) holds in general.


# Connection with Harmonic Oscillator
In this final part, we will show the connection of Hermite Polynomials with the Quantum Harmonic Oscillator. First of all, the analogue of the classical Harmonic Oscillator in Quantum Mechanics is described by the $Schr\ddot{o}dinger$ equation

<p align='center'>
    $$ -\frac{\hbar^2}{2m} \frac{d^2 \psi}{dx^2} + \frac{1}{2} m \omega^2 x^2 \psi = E \psi \quad (4)$$
</p>

There are a bunch of constants sitting in eq.(4) and life is simpler if we can just get rid of them. To this end, define
<p align='center'>
    $$ y = \sqrt{\frac{m \omega}{\hbar} x} \quad \mbox{and} \quad \tilde{E} = \frac{2E}{\hbar \omega} \quad (5)$$
</p>

Then the $Schr\ddot{o}dinger$ equation takes the cleaner form
<p align='center'>
    $$\frac{d^2 \psi}{dy^2} - y^2 \psi = - \tilde{E} \psi \quad (6)$$
</p>

The derivatives are $\psi' = -y \psi$ and $\psi" = y^2\psi - \psi$, so we see that this obeys the $Schr\ddot{o}inger$ with (rescaled) energy $\tilde{E} = 1$.

Furthermore, it's simple to see that all normalisable solutions should fall off in the same exponential fashion, with $\psi \thicksim e^{-y^2 /2}$ as $y \rightarrow \pm \infty$. This follows from looking at the large $y$ behaviour of (eq.5), where the $\tilde{E} \psi$ term is necessarily neglible compared to the $y^2 \psi$. This motivated the general ansatz

<p align='center'>
    $$\psi (y) = h(y) e^{-y^2 /2} \quad (7)$$
</p>

In general, the functions $h(y)$ are known as _Hermite Polynomials_ and have a number of nice properties.

## The Wavefunctions
The wavefunctions for the quantum harmonic oscillator contain the Gaussian form, which allows them to satisfy the necessary boundary conditions at infinity. In the wavefunction associated with a given value of the quantum number $n$, the Gaussian is multiplied by a polynomial of order $n$ called a __Hermite polynomial__. The expressions are simplified by making the substitution,

<p align='center'>
    $$\psi (y) = e^{-y^2 / 2} \quad (8)$$
</p>

Because of the association of the wavefunction with a probability density, it is necessary for the wavefunction to include a normalization constant, $N_n$.

<p align= 'center'>
    $$N_n = \frac{1}{(2^n n! \sqrt{\pi})^{1/2}} \quad (9)$$
</p>

The final form of the harmonic oscillator wavefunction is this
<p align='center'>
    $$\psi_n (y) = N_n H_n (y) e^{-y^2 /2} \quad (10) $$
</p>

where $y = \sqrt{\alpha} x$ and $\alpha = \frac{m \omega}{\hbar}$.

The general formula for the normalized wavefunctions is 
<p align='center'>
    $$\psi_n (y) = \left(\frac{\alpha}{\pi} \right)^{1/4} \frac{1}{\sqrt{2^n n!}} H_n (y)e^{-y^2 /2} \quad (11)$$
</p>

where $H_n$ is the Hermite polynomial. First four harmonic oscillator normalized wavefunctions,
<p align='center'>
    $$\begin{align} \psi_0 &= \left(\frac{\alpha}{\pi} \right)^{1/4} e^{-y^2 /2} \\ \psi_1 &= \left(\frac{\alpha}{\pi} \right)^{1/4} \sqrt{2} y e^{-y^2 /2} \\ \psi_2 &= \left(\frac{\alpha}{\pi} \right)^{1/4} \frac{1}{\sqrt{2}} (2y^2 -1) e^{-y^2 /2} \\ \psi_3 &= \left(\frac{\alpha}{\pi} \right)^{1/4} \frac{1}{\sqrt{3}} (2y^3 - 3y) e^{-y^2 /2} \end{align}$$
</p>

All energies are proportional to $\hbar \omega$, with $\omega$ the frequency of the harmonic oscillator. The energies are

<p align='center'>
    $$E_n = \hbar \omega \left (\frac{1}{2} + n \right) \quad \mbox{with} \quad n = 0, \ 1, \ 2, \ldots \quad (12)$$
</p>

When the $Schr\ddot{o}dinger$ equation for the harmonic oscillator is solved by a series method, the solutions contain this set of polynomials, named the Hermite polynomials.

|n | $H_n (y)$ | $E_n$ |
|--|-----------|-----|
|0 | $1$         |$\frac{1}{2} \hbar \omega$ |
|1 | $2y$        |$\frac{3}{2} \hbar \omega$ |
|2 | $4y^2 -2$ |$\frac{5}{2} \hbar \omega$ |
|3| $8y^3 - 12y$ |$\frac{7}{2} \hbar \omega$ |
|4| $16y^4-48y^2 +12$ | $\frac{9}{2} \hbar \omega$ |
|5| $32y^5 - 160y^3 +120y $| $\frac{11}{2} \hbar \omega$|

The wavefunctions for the quantum harmonic oscillator contain the Gaussian form, which allows them to satisfy the necessary boundary conditions at infinity. In the wavefunction associated with a given value of the quantum number $n$, the Gaussian is multiplied by a polynomial of order $n$ (the Hermite polynomials above) and the constants necessary to normalize the wavefunctions. 
