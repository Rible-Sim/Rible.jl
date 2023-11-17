# 自然坐标
```@autodocs
Modules = [Rible.NCF]
```


# Unifying rigid bodies and rigid bars using natural coordinates

In this section, the natural coordinates
\cite{dejalonKinematicDynamicSimulation1994,pappalardoNaturalAbsoluteCoordinate2015}
are adapted for unifying the non-minimal descriptions of rigid bodies and rigid bars, which are collectively called rigid members, and indistinguishably labeled
by circled numbers \CircledSmall{1}, \CircledSmall{2}, \dots, or circled capital
letters \CircledSmall{$I$}, \CircledSmall{$J$}, \dots, etc. Thus, a quantity
with a capital subscript, such as $()_{I}$, indicates the quantity belongs to
the $I$th rigid member. 


## Rigid bodies of arbitrary shapes

### 3D rigid bodies

\begin{figure}[tbh]
	\centering
	\includegraphics[width=251pt]{3D_NC_251}
	\caption{A 3D rigid body described by four types of natural coordinates. Rigid bodies are drawn by red lines. Basic points and base vectors are colored in green.}
	\label{fig:3D_NC}
\end{figure}

Consider a tetrahedron which exemplifies an arbitrary 3D rigid body, as shown in
\cref{fig:3D_NC}, where basic points
$\bm{r}_{I,i},\bm{r}_{I,j},\bm{r}_{I,k},\bm{r}_{I,l}\in \mathbb{R}^3$ and base
vectors $\bm{u}_I,\bm{v}_I,\bm{w}_I\in \mathbb{R}^3$ are fixed on the rigid body
and expressed in the global inertial frame $Oxyz$. Four types of natural
coordinates, i.e.
$$
	\begin{aligned}
		\bm{q}_{I,\mathrm{ruvw}}&=[\bm{r}^\mathrm{T}_{I,i},\bm{u}_I^\mathrm{T},\bm{v}_I^\mathrm{T},\bm{w}_I^\mathrm{T}]^\mathrm{T},\ \  \bm{q}_{I,\mathrm{rrvw}}=[\bm{r}^\mathrm{T}_{I,i},\bm{r}^\mathrm{T}_{I,j},\bm{v}_I^\mathrm{T},\bm{w}_I^\mathrm{T}]^\mathrm{T},\\
		\bm{q}_{I,\mathrm{rrrw}}&=[\bm{r}^\mathrm{T}_{I,i},\bm{r}^\mathrm{T}_{I,j},\bm{r}^\mathrm{T}_{I,k},\bm{w}_I^\mathrm{T}]^\mathrm{T},\ \  \mathrm{and}\ \ 
		\bm{q}_{I,\mathrm{rrrr}}=[\bm{r}^\mathrm{T}_{I,i},\bm{r}^\mathrm{T}_{I,j},\bm{r}^\mathrm{T}_{I,k},\bm{r}^\mathrm{T}_{I,l}]^\mathrm{T}\in \mathbb{R}^{12},
	\end{aligned}
$$
can be used to describe a 3D rigid body, corresponding to \cref{fig:3D_NC} (a)
to (d), respectively, where $()_\mathrm{ruvw}$, etc, denote the type of natural
coordinates. For the latter three types of natural coordinates, we can formally
define $\bm{u}_I=\bm{r}_{I,j}{-}\bm{r}_{I,i}$,
$\bm{v}_I=\bm{r}_{I,k}{-}\bm{r}_{I,i}$, and
$\bm{w}_I=\bm{r}_{I,l}{-}\bm{r}_{I,i}$, so that they can be converted to the
first type by 
$$
	\bm{q}_{I,\mathrm{ruvw}}=
	\bm{Y}_{\mathrm{ruvw}}\bm{q}_{I,\mathrm{ruvw}}=
	\bm{Y}_{\mathrm{rrvw}}\bm{q}_{I,\mathrm{rrvw}}=
	\bm{Y}_{\mathrm{rrrw}}\bm{q}_{I,\mathrm{rrrw}}=
	\bm{Y}_{\mathrm{rrrr}}\bm{q}_{I,\mathrm{rrrr}},
$$
where the conversion matrices are defined as, respectively,
$$
\bm{Y}_{\mathrm{ruvw}}=\begin{bmatrix*}[r]
			1& 0& 0& 0\\
			0& 1& 0& 0\\
			0& 0& 1& 0\\
			0& 0& 0& 1\\
		\end{bmatrix*}
$$
$$
	\begin{aligned}
		\bm{Y}_{\mathrm{ruvw}}=\begin{bmatrix*}[r]
			1& 0& 0& 0\\
			0& 1& 0& 0\\
			0& 0& 1& 0\\
			0& 0& 0& 1\\
		\end{bmatrix*}\otimes\mathbf{I}_3,\,
		\bm{Y}_{\mathrm{rrvw}}=\begin{bmatrix*}[r]
			1& 0& 0& 0\\
			-1& 1& 0& 0\\
			0& 0& 1& 0\\
			0& 0& 0& 1\\
		\end{bmatrix*}\otimes\mathbf{I}_3,\, 
		\bm{Y}_{\mathrm{rrrw}}=\begin{bmatrix*}[r]
			1& 0& 0& 0\\
			-1& 1& 0& 0\\
			-1& 0& 1& 0\\
			0& 0& 0& 1\\
		\end{bmatrix*}\otimes\mathbf{I}_3,\, \mathrm{and}\,
		\bm{Y}_{\mathrm{rrrr}}=\begin{bmatrix*}[r]
			1& 0& 0& 0\\
			-1& 1& 0& 0\\
			-1& 0& 1& 0\\
			-1& 0& 0& 1\\
		\end{bmatrix*}\otimes\mathbf{I}_3
	\end{aligned}
$$
where $\mathbf{I}_3$ is a $3\times3$ identity matrix, and $\otimes$ denotes the
Kronecker product.

Note that the base vectors are assumed to be non-coplanar, thus the natural
coordinates actually form an affine frame attached to the 3D rigid body.
Consequently, the position vector of a generic point on the 3D rigid body can be
expressed by
$$
	\bm{r}=\bm{r}_{I,i}+c_{I,1}\bm{u}_I+c_{I,2}\bm{v}_I+c_{I,3}\bm{w}_I=
	\bm{C}_{I,\mathrm{body}}\bm{q}_{I,\mathrm{body}},
$$
where $c_{I,1}$, $c_{I,2}$ and $c_{I,3}$ are the affine coordinates;
$\bm{C}_{I,\mathrm{body}}=\left(\left[1, c_{I,1},
c_{I,2}, c_{I,3}\right] \otimes \mathbf{I}_3\right)\bm{Y}_{\mathrm{body}}$ is a transformation matrix for
$\bm{q}_{I,\mathrm{body}}$; 
$()_{\mathrm{body}}$ can be any of 
$()_{\mathrm{ruvw}}$, $()_{\mathrm{rrvw}}$, $()_{\mathrm{rrrw}}$, or $()_{\mathrm{rrrr}}$.

To ensure rigidity of the body, the natural coordinates
$\bm{q}_{I,\mathrm{body}}$ must satisfy six intrinsic constraints
$$
	\begin{aligned}
		\bm{\varPhi}_I( \bm{q}_{I,\mathrm{body}} ) =
		\begin{pmatrix}	
			\bm{u}_I^\mathrm{T}\bm{u}_I-\bar{\bm{u}}_I^\mathrm{T}\bar{\bm{u}}_I\\	
			\bm{v}_I^\mathrm{T}\bm{v}_I-\bar{\bm{v}}_I^\mathrm{T}\bar{\bm{v}}_I\\	
			\bm{w}_I^\mathrm{T}\bm{w}_I-\bar{\bm{w}}_I^\mathrm{T}\bar{\bm{w}}_I\\	
			\bm{v}_I^\mathrm{T}\bm{w}_I-\bar{\bm{v}}_I^\mathrm{T}\bar{\bm{w}}_I\\	
			\bm{u}_I^\mathrm{T}\bm{w}_I-\bar{\bm{u}}_I^\mathrm{T}\bar{\bm{w}}_I\\	
			\bm{u}_I^\mathrm{T}\bm{v}_I-\bar{\bm{u}}_I^\mathrm{T}\bar{\bm{v}}_I\\
		\end{pmatrix}=\bm{0}
		%		\ \  \mathrm{and}\ \ 
		%		\bm{A}_I( \bm{q}_{I,\mathrm{body}} ) =\begin{bmatrix}
			%			\cdot&		2\bm{u}_I^\mathrm{T}&		\cdot&		\cdot\\
			%			\cdot&		\cdot&		2\bm{v}_I^\mathrm{T}&		\cdot\\
			%			\cdot&		\cdot&		\cdot&		2\bm{w}_I^\mathrm{T}\\	
			%			\cdot&		\cdot&		\bm{w}_I^\mathrm{T}&		\bm{v}_I^\mathrm{T}\\
			%			\cdot&		\bm{w}_I^\mathrm{T}&		\cdot&		\bm{u}_I^\mathrm{T}\\
			%			\cdot&		\bm{v}_I^\mathrm{T}&		\bm{u}_I^\mathrm{T}&		\cdot\\
			%		\end{bmatrix}.
	\end{aligned}	
$$
where $\bar{\bm{u}}_I$, $\bar{\bm{v}}_I$ and $\bar{\bm{w}}_I$ are constant
vectors in a local frame, which is fixed on the rigid member (See also
\cref{sec:uni_mass}). Then, the position and orientation of a 6-DoF 3D rigid body
can be defined by twelve coordinates (any type in \cref{eq:3D_natural}) and six
constraints \cref{eq:3D_intrinsic}.
%done remove 2D bar
### 3D rigid bars
\begin{figure}[tbh]
	\centering
	\includegraphics[width=162pt]{Bar_162}
	\caption{A 3D rigid bar described by two types of natural coordinates.}
	\label{fig:Bar}
\end{figure}
Two types of natural coordinates, i.e. 
$\bm{q}_{I,\mathrm{ru}}=[\bm{r}^\mathrm{T}_{I,i},\bm{u}_I^\mathrm{T}]^\mathrm{T}$
and 
$\bm{q}_{I,\mathrm{rr}}=[\bm{r}_{I,i}^\mathrm{T},\bm{r}_{I,j}^\mathrm{T}]^\mathrm{T}\in\mathbb{R}^{6}$,
can describe a 3D rigid bar, corresponding to \cref{fig:Bar} (a) and (b), respectively.
Define conversion matrices 
$$	
	\bm{Y}_\mathrm{ru}=\begin{bmatrix}[r]
		1& 0\\
		0& 1\\
	\end{bmatrix}\otimes \bm{\mathrm{I}}_3

	\bm{Y}_\mathrm{rr}=\begin{bmatrix}[r]
		1& 0\\
		-1& 1\\
	\end{bmatrix}\otimes \bm{\mathrm{I}}_3
$$
Then, the position vector of a generic point along the longitudinal axis of the rigid
bar is given by
$$
	\bm{r}=\bm{r}_{I,i}+c_I\bm{u}_I=\bm{C}_{I,\mathrm{bar}}\bm{q}_{I,\mathrm{bar}}
$$
where the coefficient $c_I$ depends on the relative position of the generic point;
$\bm{C}_{I,bar}=\left(\left[1, c_I\right] \otimes \mathbf{I}_3\right)\bm{Y}_{\mathrm{bar}}$ is the
transformation matrix for $\bm{q}_{I,\mathrm{bar}}$; 
$()_{\mathrm{bar}}$ can be either $()_{\mathrm{ru}}$ or $()_{\mathrm{rr}}$.
And the intrinsic constraint to preserve the bar length is
$$
	\varPhi_I( \bm{q}_{I,\mathrm{bar}} )=\bm{u}_I^\mathrm{T}\bm{u}_I-\bar{\bm{u}}_I^\mathrm{T}\bar{\bm{u}}_I=0
	%	\ \ \mathrm{and}\ \ 
	%	\bm{A}_I( \bm{q}_{I,\mathrm{bar}} ) = \frac{\partial\varPhi_I}{\partial\bm{q}_{I,\mathrm{bar}}}=\left[
	%	\bm{0}^\mathrm{T},\ 		2\bm{u}_I^\mathrm{T}\\
	%	\right]
$$

Hence, the position and orientation of a 5-DoF 3D rigid bar can be
defined by six coordinates and one constraint \cref{eq:bar_intrinsic}.
\subsection{Unified formulations and mass matrices}\label{sec:uni_mass}

\begin{table}[tbh]
	\caption{Polymorphism of natural coordinates for rigid bodies and rigid bars}
	\label{tab:polymorphsim}
	\centering
	\setlength\tabcolsep{4pt}
	\renewcommand{\arraystretch}{1.0}
	\begin{tabular}{c|cccc}
	\toprule
	&
	{\shortstack{Degrees of \\ freedom}} & 
	{\shortstack{Number of \\ coordinates}} & 
	{\shortstack{Number of \\ constraints}} & 
	{\shortstack{Types of \\ natural coordinates}} \\
	\midrule
	3D Rigid Body & 6 & 12 & 6 & ruvw rrvw rrrw rrrr \\
	3D Rigid Bar  & 5 & 6  & 1 & ru rr \\
	\bottomrule
	\end{tabular}
\end{table}

The transformation relations \cref{eq:bar_trans,eq:3D_trans} for the
standard types of natural coordinates can be put into a unifying form
$$
	\bm{r}=\bm{C}_{I}\bm{q}_{I},
$$
which is a polymorphic expression, meaning that the formulations of
$\bm{C}_I$ and $\bm{Y}_I$ vary with the type of $\bm{q}_I$, as summarized in
\cref{tab:polymorphsim}. However, note that $\bm{C}_{I}$ is not a function of
$\bm{q}_I$. Consequently, the velocity of a generic point is given by
$\dot{\bm{r}}=\bm{C}_{I}\dot{\bm{q}}_{I}$, which can be used to derive the mass
matrix. Let $\rho_I$ denote the longitudinal or volume density of the rigid
member $\CircledSmall{I}$. Then, the kinetic energy can be computed by an
integral over its entire domain $\Omega$ as
$$\textstyle
	\begin{aligned}\textstyle
		T_I
		&=\frac{1}{2}{\textstyle\int\limits_{\Omega}}{\rho_I\bm{\dot{r}}^\mathrm{T}\bm{\dot{r}}\mathrm{d}\Omega}
		=\frac{1}{2}{\textstyle\int\limits_{\Omega}}{\rho_I \dot{\bm{q}}_{I}^\mathrm{T}\bm{C}_{I}^\mathrm{T}\bm{C}_I\dot{\bm{q}}_I\mathrm{d}\Omega}=\frac{1}{2}\dot{\bm{q}}_{I}^\mathrm{T}\bm{M}_I\dot{\bm{q}}_I
		%		\\
		%		&\textstyle=\frac{1}{2}{ \dot{\check{\bm{q}}}^{\mathrm{T}}\check{\bm{T}}_{I}^{\mathrm{T}}  \bm{M}_I \check{\bm{T}}_I\dot{\check{\bm{q}}} }+\frac{1}{2}{ \dot{\tilde{\bm{q}}}^{\mathrm{T}}\tilde{\bm{T}}_{I}^{\mathrm{T}}  \bm{M}_I \tilde{\bm{T}}_I\dot{\tilde{\bm{q}}} }
		%		+{ \dot{\check{\bm{q}}}^{\mathrm{T}}\check{\bm{T}}_{I}^{\mathrm{T}}  \bm{M}_I \tilde{\bm{T}}_I\dot{\tilde{\bm{q}}} }
	\end{aligned}
$$
where $\bm{M}_I$ is a constant mass matrix with polymorphism defined by
$$
\begin{aligned}
	\bm{M}_{I}&
	=\int _{\Omega}\rho _I\bm{C}_{I}^{\mathrm{T}}\bm{C}_{I}\mathrm{d}\Omega
	=\bm{Y}_I^{\mathrm{T}}\left(\int _{\Omega}\left( \rho _I \begin{bmatrix}
		1&		\bm{c}_{I}^{\mathrm{T}}\\
		\bm{c}_I&		\bm{c}_I\bm{c}_{I}^{\mathrm{T}}
	\end{bmatrix} \right) \mathrm{d}\Omega \otimes \mathbf{I}_3\right)\bm{Y}_I\\
	&=\bm{Y}_I^{\mathrm{T}}\left(\begin{bmatrix}
		\smallint _{\Omega}\rho _I\mathrm{d}\Omega&		\smallint _{\Omega}\rho _I\bm{c}_{I}^{\mathrm{T}}\mathrm{d}\Omega\\
		\smallint _{\Omega}\rho _I\bm{c}_I\mathrm{d}\Omega&		\smallint _{\Omega}\rho _I\bm{c}_I\bm{c}_{I}^{\mathrm{T}}\mathrm{d}\Omega\\
	\end{bmatrix} \otimes \mathbf{I}_3\right)\bm{Y}_I
\end{aligned}
$$


It is possible to express the mass matrix by conventional inertia properties, 
such as the mass, the center of mass, 
and the moments of inertia of a rigid member. 
To this end, let's introduce a local Cartesian frame $\bar{O}\bar{x}\bar{y}\bar{z}$
which is fixed on the rigid member $\CircledSmall{I}$, as shown in \cref{fig:local}. 
Quantities expressed in this local frame are denoted by an overline $\bar{()}$.
Without loss of generality, let its origin $\bar{O}$ coincide with the mass center, 
such that $\bar{\bm{r}}_{I,g}=\bm{0}$. 
For a 3D rigid body, let its axes align along the principal axes of inertia. 
For a 3D rigid bar, let its $\bar{x}$ axis aligns along the longitudinal direction.

%done remove 2D
\begin{figure}[tbh]
	\centering
	\includegraphics[width=244pt]{local_244}
	\caption{The basic point $\bar{\bm{r}}_{I,i}$, 
	the base vectors $\bar{\bm{u}}_{I}$, $\bar{\bm{v}}_{I}$ and $\bar{\bm{w}}_{I}$, 
	the mass center $\bar{\bm{r}}_{I,g}$, 
	and a generic point $\bar{\bm{r}}_{I}$ in the local Cartesian frame of 
	(a) a 3D rigid body or (b) a 3D rigid bar.}\label{fig:local}
\end{figure}

Because the basic points and base vectors are fixed on the rigid members, their
coordinates in the local frame are constant. Let's define a polymorphic matrix
$$
    \bar{\bm{X}}_I=\left\{
	\begin{aligned}
		&\left[\bar{\bm{u}},\bar{\bm{v}},\bar{\bm{w}}\right] &&\mathrm{for~a~rigid~body}\\
		&\left[\bar{\bm{u}}\right] &&\mathrm{for~a~rigid~bar}
	\end{aligned}\right.
$$
Then, according to \cref{eq:trans}, the position vector of a generic point 
in the local frame can be expressed by
$\bar{\bm{r}} = \bar{\bm{r}}_{I,i} + \bar{\bm{X}}_I \bm{c}_I$, which gives
$$
	\bm{c}_I=\bar{\bm{X}}^{+}_I( \bar{\bm{r}}-\bar{\bm{r}}_{I,i} )
$$
where $()^{+}$ denotes the Moore-Penrose pseudoinverse. For \cref{eq:X_uvw}, because the columns are linearly independent, i.e. $\bar{\bm{X}}$ has full rank, the pseudoinverse is equal to the matrix inverse.

Using \cref{eq:c}, the following expressions for use in \cref{eq:M} can be derived:
$$
	\begin{aligned}
		\int _{\Omega}\rho _I\mathrm{d}\Omega &=m_I\\
		\int _{\Omega}\rho _I\bm{c}_I\mathrm{d}\Omega &= m_I\bar{\bm{X}}^+\left( \bar{\bm{r}}_{I,g}-\bar{\bm{r}}_{I,i} \right)=-m_I\bar{\bm{X}}^+\bar{\bm{r}}_{I,i}\\
		\int _{\Omega}\rho _I\bm{c}_I\bm{c}_{I}^{\mathrm{T}}\mathrm{d}\Omega &=
		\bar{\bm{X}}^+\left( \bar{\bm{J}}_I -m_I\bar{\bm{r}}_{I,i}\bar{\bm{r}}_{I,g}^{\mathrm{T}}-m_I\bar{\bm{r}}_{I,g}\bar{\bm{r}}_{I,i}^{\mathrm{T}}+m_I\bar{\bm{r}}_{I,i}\bar{\bm{r}}_{I,i}^{\mathrm{T}} \right) \bar{\bm{X}}^{+\mathrm{T}}\\
		&=\bar{\bm{X}}^+\left( \bar{\bm{J}}_I +m_I\bar{\bm{r}}_{I,i}\bar{\bm{r}}_{I,i}^{\mathrm{T}} \right) \bar{\bm{X}}^{+\mathrm{T}}
	\end{aligned}
$$
where $m_I$ is the mass of the rigid member $\CircledSmall{I}$; 
$\bar{\bm{J}}_I$ contains the moments of inertia and necessitates some discussions:

For a 3D rigid body, 
$\bar{\bm{J}}_I$ is given by 
$$
	\bar{\bm{J}}_I 
	=\int_{\Omega}\rho_I\bar{\bm{r}}\bar{\bm{r}}^{\mathrm{T}}\mathrm{d}\Omega
	=\int_{\Omega}\rho_I\begin{bmatrix}
		\bar{x}^2     & \bar{y}\bar{x}& \bar{z}\bar{x}\\
		\bar{x}\bar{y}&	     \bar{y}^2& \bar{z}\bar{y}\\
		\bar{x}\bar{z}& \bar{y}\bar{z}&      \bar{z}^2\\
	\end{bmatrix} \mathrm{d}\Omega,
$$
while the conventional inertia matrix is given by 
$$
	\bar{\bm{I}}_I=\int_{\Omega}\rho_I\begin{bmatrix}
		\bar{y}^2+\bar{z}^2     & -\bar{y}\bar{x}& -\bar{z}\bar{x}\\
		-\bar{x}\bar{y}&	     \bar{x}^2+\bar{z}^2& -\bar{z}\bar{y}\\
		-\bar{x}\bar{z}& -\bar{y}\bar{z}&      \bar{x}^2+\bar{y}^2\\
	\end{bmatrix} \mathrm{d}\Omega
$$ 
Hence, we have $\bar{\bm{J}}_I =\tfrac{1}{2}\mathrm{trace}\left( \bar{\bm{I}}_I \right) \mathbf{I}_3-\bar{\bm{I}}_I$.

For a 3D rigid bar, the expression of $\bar{\bm{J}}_I$ is the same as \cref{eq:J_3D}, except that only the element $\bar{x}^2$ is nonzero. And the pseudoinverse of $\bar{\bm{X}}_{I}=\left[ \bar{u}_x,0,0 \right] ^{\mathrm{T}}$ is $\bar{\bm{X}}_{I}^{+}=\left[{1}/{\bar{u}_x},0,0\right]$. Therefore, we have $\bar{\bm{X}}_{I}^{+}\bar{\bm{J}}_I\bar{\bm{X}}_{I}^{+\mathrm{T}}=\left(\int_{\Omega}\rho_I\bar{x}^2\mathrm{d}\Omega\right)/{\bar{u}^2_x}$.

For an advanced treatment of the inertia representation for rigid multibody systems in terms of natural coordinates, 
we refer the interested readers to our previous paper \cite{xuGeneralizedInertiaRepresentation2021}. 

