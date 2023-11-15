## Unit quaternion
Quaternion
$$\bm{q}=\left[ \begin{matrix}	q_0&
	q_1&
	q_2&
	q_3\\\end{matrix} \right] ^{\mathrm{T}}$$

Unit quaternion constraint
$$\varPhi \left( \bm{q} \right) 
=\frac{1}{2}\left( \bm{q}^{\mathrm{T}}\bm{q}-1 \right) $$

Constraint Jacobian
$$\bm{A}\left( \bm{q} \right) =\varPhi 
_{\bm{q}}\left( \bm{q} \right) =\bm{q}^{\mathrm{T}}$$

Constraint force Jacobian
$$\frac{\partial \bm{A}^{\mathrm{T}}\lambda}{\partial 
\bm{q}}=\frac{\partial}{\partial \bm{q}}\left( \lambda 
\bm{q} \right) =\lambda \mathbf{I}$$

Quaternion velocity must satisfy
$$\bm{A}\left( \bm{q} \right) 
\dot{\bm{q}}=\bm{q}^{\mathrm{T}}\dot{\bm{q}}=0$$

Orthonormal basis of nullspace of constraint Jacobian
$$\mathcal{N} \left( \bm{A} \right) 
=\mathrm{span}\left( \bm{N}\left( \bm{q} \right) \right) 
\\\bm{N}\left( \bm{q} \right) 
=\bm{L}^{\mathrm{T}}\left( \bm{q} \right) $$

where
$$\bm{L}\left( \bm{q} \right) =\left[ \begin{array}{r}
	-q_1&	
	q_0&	
	q_3&	
	-q_2\\	-q_2&
		-q_3&
		q_0&
		q_1\\
	-q_3&	
	q_2&	
	-q_1&	
	q_0\\\end{array} \right] $$
	
By orthonormality of nullspace
$$\bm{L}\left( \bm{q} \right) 
\bm{L}^{\mathrm{T}}\left( \bm{q} \right) =\mathbf{I}_3$$

By definition of nullspace
$$\bm{A}\left( \bm{q} \right) \bm{N}\left( \bm{q} 
\right) 
=\bm{q}^{\mathrm{T}}\bm{L}^{\mathrm{T}}\left( \bm{q} 
\right) =\mathbf{0}$$

or
$$\bm{L}\left( \bm{q} \right) \bm{q}=\mathbf{0}$$

Useful identity
$$\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{L}\left( \bm{q} \right) 
=\mathbf{I}_4-\bm{qq}^{\mathrm{T}}$$

Expressions (1.8) to (1.11) are valid for ANY 4-vector $\bm{q}$, including 
$\dot{\bm{q}}$.

Derivative of (1.10) leads to a useful identity

$$\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}+\bm{L}\left( \dot{\bm{q}} \right) 
\bm{q}=\mathbf{0}\\\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}=-\bm{L}\left( \dot{\bm{q}} \right) 
\bm{q}$$

Local angular velocity $\bm{\varOmega }$ as independent velocity

$$\dot{\bm{q}}=\frac{1}{2}\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{\varOmega }$$

Using (1.8), (1.13) leads to 
$$\bm{\varOmega }=2\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}$$

By use of (1.11), (1.14) leads to another useful identity

$$\bm{\varOmega }^{\mathrm{T}}\bm{\varOmega }=4\dot{\bm
{q}}^{\mathrm{T}}\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}\\=4\dot{\bm{q}}^{\mathrm{T}}\left( \mathbf{I}_4-\bm{qq}^{\mathrm{T}} \right) 
\dot{\bm{q}}\\=4\dot{\bm{q}}^{\mathrm{T}}\dot{\bm{q}}$$

## Modified inerita representation

Splitted inertia with nonzero parameter $\gamma$
$$\bm{J}=\bm{J}-\gamma \mathbf{I}_3+\gamma 
\mathbf{I}_3=\bm{J}_{\gamma}+\gamma \mathbf{I}_3$$

where
$$ \bm{J}_{\gamma} = \bm{J}-\gamma \mathbf{I}_3$$


Modified rotational kinetic energy, using
$$T_{\mathrm{rot},\gamma}=\frac{1}{2}\bm{\varOmega }^{\mathrm{T}}\bm{J\varOmega }\\=\frac{1}{2}\bm{\varOmega }^{\mathrm{T}}\bm{J}_{\gamma}\bm{\varOmega }+\gamma 
\frac{1}{2}\bm{\varOmega }^{\mathrm{T}}\bm{\varOmega }\\=2\dot{
\bm{q}}^{\mathrm{T}}\bm{L}^{\mathrm{T}}\left( \bm{q} 
\right) \bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}+2\gamma 
\dot{\bm{q}}^{\mathrm{T}}\dot{\bm{q}}$$	

Modified rotational momentum
$$\bm{p}_{\mathrm{rot},\gamma}=\frac{\partial 
T_{\mathrm{rot},\gamma}}{\partial 
\dot{\bm{q}}^{\mathrm{T}}}=\bm{M}_{\mathrm{rot},\gamma}\left( \bm{q} \right) \dot{\bm{q}}$$

where modified rotational mass matrix
$$\bm{M}_{\mathrm{rot},\gamma}=4\bm{L}^{\mathrm{T}}\left( \bm{q} \right) \bm{J}_{\gamma}\bm{L}\left( \bm{q} 
\right) +4\gamma \mathbf{I}_4$$

$$\bm{M}_{\mathrm{rot},\gamma}\dot{\bm{q}}=4\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}+4\gamma 
\mathbf{I}_4\dot{\bm{q}}\\=-4\bm{L}^{\mathrm{T}}\left( \bm{q} \right) \bm{J}_{\gamma}\bm{L}\left( \dot{\bm{q}} 
\right) \bm{q}+4\gamma \mathbf{I}_4\dot{\bm{q}}\\\frac{\partial 
\left( \bm{M}_{\mathrm{rot},\gamma}\dot{\bm{q}} 
\right)}{\partial 
\bm{q}}_1=-4\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{J}_{\gamma}\bm{L}\left( \dot{\bm{q}} \right) 
\\\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}=\bm{\eta }\\4\bm{L}^{\mathrm{T}}\left( \bm{q} \right) \bm{\eta }=4\left[ \begin{array}{r}	-q_1&
		-q_2&
		-q_3\\
	\mathrm{ }q_0&	
	-q_3&	
	\mathrm{ }q_2\\	\mathrm{ }q_3&
		\mathrm{ }q_0&
		-q_1\\
	-q_2&	
	\mathrm{ }q_1&	
	\mathrm{ }q_0\\\end{array} \right] \left[ \begin{array}{c}	\eta _1\\
	\eta _2\\	\eta _3\\\end{array} \right] 
=4\left[ \begin{array}{r}	-q_1\eta _1-q_2\eta _2-q_3\eta _3\\	q_0\eta _1-q_3\eta 
_2+q_2\eta _3\\	q_3\eta _1+q_0\eta _2-q_1\eta _3\\	-q_2\eta _1+q_1\eta 
_2+q_0\eta _3\\\end{array} \right] \\\frac{\partial 
\left( \bm{M}_{\mathrm{rot},\gamma}\dot{\bm{q}} 
\right)}{\partial \bm{q}}_2=4\left[ \begin{array}{r}	0&
		-\eta _1&
		-\eta _2&
		-\eta _3\\
	\eta _1&	
	0&	
	\eta _3&	
	-\eta _2\\	\eta _2&
		-\eta _3&
		0&
		\eta _1\\
	\eta _3&	
	\eta _2&	
	-\eta _1&	
	0\\\end{array} \right] \\\frac{\partial 
\left( \bm{M}_{\mathrm{rot},\gamma}\dot{\bm{q}} 
\right)}{\partial \bm{q}}=\frac{\partial 
\left( \bm{M}_{\mathrm{rot},\gamma}\dot{\bm{q}} 
\right)}{\partial \bm{q}}_1+\frac{\partial 
\left( \bm{M}_{\mathrm{rot},\gamma}\dot{\bm{q}} 
\right)}{\partial 
\bm{q}}_2=-4\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{J}_{\gamma}\bm{L}\left( \dot{\bm{q}} \right) 
+4\left[ \begin{array}{r}	0&	
	-\eta _1&	
	-\eta _2&	
	-\eta _3\\	\eta _1&
		0&
		\eta _3&
		-\eta _2\\
	\eta _2&	
	-\eta _3&	
	0&	
	\eta _1\\	\eta _3&
		\eta _2&
		-\eta _1&
		0\\\end{array} \right] $$
		
The derivative of mass matrix
$$\dot{\bm{M}}_{\mathrm{rot},\gamma}=4\bm{L}^{\mathrm{T}}\left(
 \bm{q} \right) 
\bm{J}_{\gamma}\bm{L}\left( \dot{\bm{q}} \right) 
+4\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} \right) 
\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) $$

The partial derivative of rotational kinetic energy
$$\frac{\partial T_{\mathrm{rot},\gamma}}{\partial 
\bm{q}^{\mathrm{T}}}=\frac{\partial}{\partial 
\bm{q}}\left( 2\bm{q}^{\mathrm{T}}\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} \right) 
\bm{J}_{\gamma}\bm{L}\left( \dot{\bm{q}} \right) 
\bm{q} \right) \\=4\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} 
\right) \bm{J}_{\gamma}\bm{L}\left( \dot{\bm{q}} \right) 
\bm{q}\\=-4\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} \right) 
\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}$$

Quadratic velocity term for Centrifugal and Coriolis forces, using

$$\dot{\bm{M}}_{\mathrm{rot},\gamma}\dot{\bm{q}}-\frac{\partial 
T_{\gamma}}{\partial 
\bm{q}^{\mathrm{T}}}=\left( 4\bm{L}^{\mathrm{T}}\left( \bm{q} \right) \bm{J}_{\gamma}\bm{L}\left( \dot{\bm{q}} 
\right) +4\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} \right) 
\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) \right) 
\dot{\bm{q}}+4\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} 
\right) \bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}\\=8\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} 
\right) \bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}$$

## rigid body
$$\bm{x}=\left[ \begin{matrix}
	\bm{r}_{O}^{\mathrm{T}}&	
	\bm{q}^{\mathrm{T}}\\\end{matrix} \right] ^{\mathrm{T}}$$
	
Generic point

$$\bm{r}=\bm{r}_O+\bm{R}\left( \bm{q} \right) 
\bm{c}$$

where rotation matrix

$$\bm{R}\left( \bm{q} \right) =\left[ \begin{matrix}
	q_{0}^{2}+q_{1}^{2}-q_{2}^{2}-q_{3}^{2}&	
	2\left( q_1q_2-q_0q_3 \right)&	
	2\left( q_1q_3+q_0q_2 \right)\\	2\left( q_1q_2+q_0q_3 
\right)&	
	q_{0}^{2}-q_{1}^{2}+q_{2}^{2}-q_{3}^{2}&	
	2\left( q_2q_3-q_0q_1 \right)\\	2\left( q_1q_3-q_0q_2 
\right)&		2\left( q_2q_3+q_0q_1 \right)&	
	q_{0}^{2}-q_{1}^{2}-q_{2}^{2}+q_{3}^{2}\\\end{matrix} \right] 
\\\bm{R}\left( \bm{q} \right) =\left[ \begin{matrix}
	q_{0}^{2}+q_{1}^{2}-q_{2}^{2}-q_{3}^{2}&	
	2\left( q_1q_2 \right)&	
	2\left( q_1q_3 \right)\\	2\left( q_1q_2 \right)&	
	q_{0}^{2}-q_{1}^{2}+q_{2}^{2}-q_{3}^{2}&	
	2\left( q_2q_3 \right)\\	2\left( q_1q_3 \right)&	
	2\left( q_2q_3 \right)&	
	q_{0}^{2}-q_{1}^{2}-q_{2}^{2}+q_{3}^{2}\\\end{matrix} \right] 
+2\left[ \begin{array}{r}	&	
	-q_0q_3&	
	q_0q_2\\	q_0q_3&
		&
		-q_0q_1\\
	-q_0q_2&	
	q_0q_1&	
	\\\end{array} \right] \\$$
	
Generic point velocity and variation
$$\dot{\bm{r}}=\dot{\bm{r}}_O+\dot{\bm{R}}\left( \bm{q} \right) 
\bm{c}\\=\dot{\bm{r}}_O+\bm{R}\left( \bm{q} 
\right) \hat{\bm{\varOmega}}\left( \bm{q} \right) 
\bm{c}\\=\dot{\bm{r}}_O-\bm{R}\left( \bm{q} 
\right) \hat{\bm{c}}\bm{\varOmega }\left( \bm{q} \right) 
\\=\dot{\bm{r}}_O-\bm{R}\left( \bm{q} \right) 
\hat{\bm{c}}2\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}\\=\bm{C}\left( \bm{x} \right) 
\left( \begin{array}{c}	\dot{\bm{r}}_0\\
	\dot{\bm{q}}\\\end{array} \right) $$
	
$$\delta \bm{r}=\bm{C}\left( \bm{x} \right) \delta 
\bm{x}$$

where

$$\bm{C}\left( \bm{x};\bm{c} \right) 
=\left[ \begin{matrix}	\mathbf{I}&	
	-2\bm{R}\left( \bm{q} \right) 
\hat{\bm{c}}\bm{L}\left( \bm{q} \right)\\\end{matrix} 
\right] $$

$$\bm{R}\left( \bm{q} \right) 
\bm{\eta }=\left( \begin{array}{c}
	\left( q_{0}^{2}+q_{1}^{2}-q_{2}^{2}-q_{3}^{2} \right) \eta 
_1+2\left( q_1q_2-q_0q_3 \right) \eta _2+2\left( q_1q_3+q_0q_2 \right) \eta _3\\
	2\left( q_1q_2+q_0q_3 \right) \eta 
_1+\left( q_{0}^{2}-q_{1}^{2}+q_{2}^{2}-q_{3}^{2} \right) \eta 
_2+2\left( q_2q_3-q_0q_1 \right) \eta _3\\	2\left( q_1q_3-q_0q_2 \right) \eta 
_1+2\left( q_2q_3+q_0q_1 \right) \eta 
_2+\left( q_{0}^{2}-q_{1}^{2}-q_{2}^{2}+q_{3}^{2} \right) \eta _3\\\end{array} 
\right) \\\frac{\partial \bm{R}\left( \bm{q} \right) 
\underline{\bm{\eta }}}{\partial \bm{q}}=\left[ \begin{array}{r}
	2q_0\eta _1-2q_3\eta _2+2q_2\eta _3&	 {2q_1\eta _1+2q_2\eta _2+2q_3\eta _3}&	 {-2q_2\eta _1+2q_1\eta _2+2q_0\eta _3}&	
 {2q_3\eta _1-2q_0\eta _2+2q_1\eta _3}\\	2q_3\eta 
_1+2q_0\eta _2-2q_1\eta _3& {2q_2\eta 
_1-2q_1\eta _2-2q_0\eta _3}& {2q_1\eta 
_1+2q_2\eta _2+2q_3\eta _3}& {2q_0\eta 
_1-2q_3\eta _2+2q_2\eta _3}\\	-2q_2\eta _1+2q_1\eta _2+2q_0\eta _3&	
{2q_3\eta _1+2q_0\eta _2-2q_1\eta _3}&	
{-2q_0\eta _1+2q_3\eta _2-2q_2\eta _3}&	
{2q_1\eta _1+2q_2\eta _2+2q_3\eta _3}\\\end{array} \right] 
\\2\left( q_0\mathbf{I}+\left[ \begin{array}{r}	&
		-q_3&
		q_2\\
	q_3&	
	&	
	-q_1\\	-q_2&
		q_1&
		\\\end{array} \right] \right) 
\left( \begin{array}{c}	\eta _1\\	\eta _2\\
	\eta _3\\\end{array} \right) \\2\left( q_1\eta _1+q_2\eta _2+q_3\eta _3 \right) 
\mathbf{I}-2\left( \left( q_0\mathbf{I}+\left[ \begin{array}{r}	&
		-q_3&
		q_2\\
	q_3&	
	&	
	-q_1\\	-q_2&
		q_1&
		\\\end{array} \right] \right) 
\left( \begin{array}{c}	\eta _1\\	\eta _2\\
	\eta _3\\\end{array} \right) \right) _{\times}\\\frac{\partial}{\partial 
\bm{q}}\left( \bm{R}\left( \bm{q} \right) 
\hat{\bm{c}}2\bm{L}\left( \dot{\bm{q}} \right) 
\underline{\bm{q}} \right) =\\\frac{\partial}{\partial 
\bm{q}}\left( -\bm{R}\left( \bm{q} \right) 
\hat{\bm{c}}2\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}} \right) =2\frac{\partial}{\partial 
\bm{q}}\left( \bm{R}\left( \bm{q} \right) 
\hat{\bm{c}}\bm{L}\left( \dot{\bm{q}} \right) 
\bm{q} \right) \\=2\left( \bm{R}\left( \bm{q} \right) 
\hat{\bm{c}}\bm{L}\left( \dot{\bm{q}} \right) 
-\frac{\partial}{\partial 
\bm{q}}\left( \bm{R}\left( \bm{q} \right) 
\hat{\bm{c}}\bm{L}\left( \underline{\bm{q}} \right) 
\dot{\bm{q}} \right) \right) \\\frac{\partial 
\dot{\bm{r}}}{\partial \bm{x}}=\frac{\partial 
\bm{C}\left( \bm{x} \right) 
\underline{\dot{\bm{x}}}}{\partial 
\bm{x}}\\=\frac{\partial}{\partial 
\bm{x}}\left( \dot{\bm{r}}_O-\bm{R}\left( \bm{q} 
\right) \hat{\bm{c}}2\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}} \right) \\=\left[ \begin{matrix}	\bm{0}&
		\frac{\partial}{\partial 
\bm{q}}\left( -\bm{R}\left( \bm{q} \right) 
\hat{\bm{c}}2\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}} \right)\\\end{matrix} \right] $$

Virtual work of a concentrated force
$$\delta W\,\,=\,\,\bm{f}^{\mathrm{T}}\delta 
\bm{r}=\bm{f}^{\mathrm{T}}\left[ \begin{matrix}	\mathbf{I}&
		-2\bm{R}\left( \bm{q} 
\right) \hat{\bm{c}}\bm{L}\left( \bm{q} 
\right)\\\end{matrix} \right] \delta 
\bm{x}\\=\bm{F}^{\mathrm{T}}\delta \bm{x}$$

where the generalized concentrated force is
$$\bm{F}=\left[ \begin{array}{c}	\mathbf{I}\\
	-2\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\hat{\bm{c}}^{\mathrm{T}}\bm{R}^{\mathrm{T}}\left( \bm{
q} \right)\\\end{array} \right] \bm{f}$$

$$\frac{\partial}{\partial 
\bm{q}}\left( \bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{c} \right) =\frac{\partial}{\partial 
\bm{q}}\left( \begin{array}{r}	-q_1c_1-q_2c_2-q_3c_3\\
	q_0c_1-q_3c_2+q_2c_3\\	q_3c_1+q_0c_2-q_1c_3\\
	-q_2c_1+q_1c_2+q_0c_3\\\end{array} \right) =\left[ \begin{array}{r}	0&
		-c_1&
		-c_2&
		-c_3\\
	c_1&	
	0&	
	c_3&	
	-c_2\\	c_2&
		-c_3&
		0&
		c_1\\
	c_3&	
	c_2&	
	-c_1&	
	0\\\end{array} \right] \\\bm{R}\left( \bm{q} \right) 
\bm{f}=\left( \begin{array}{c}
	\left( q_{0}^{2}+q_{1}^{2}-q_{2}^{2}-q_{3}^{2} \right) \eta 
_1+2\left( q_1q_2-q_0q_3 \right) \eta _2+2\left( q_1q_3+q_0q_2 \right) \eta _3\\
	2\left( q_1q_2+q_0q_3 \right) \eta 
_1+\left( q_{0}^{2}-q_{1}^{2}+q_{2}^{2}-q_{3}^{2} \right) \eta 
_2+2\left( q_2q_3-q_0q_1 \right) \eta _3\\	2\left( q_1q_3-q_0q_2 \right) \eta 
_1+2\left( q_2q_3+q_0q_1 \right) \eta 
_2+\left( q_{0}^{2}-q_{1}^{2}-q_{2}^{2}+q_{3}^{2} \right) \eta _3\\\end{array} 
\right) \\\bm{R}\left( \bm{q} \right) =\left[ \begin{matrix}
	\left( q_{0}^{2}+q_{1}^{2}-q_{2}^{2}-q_{3}^{2} \right)&	
	2q_1q_2&	
	2q_1q_3\\	2q_1q_2&
	
	\left( q_{0}^{2}-q_{1}^{2}+q_{2}^{2}-q_{3}^{2} \right)&	
	2q_2q_3\\	2q_1q_3&
		2q_2q_3&
	
	\left( q_{0}^{2}-q_{1}^{2}-q_{2}^{2}+q_{3}^{2} \right)\\\end{matrix} \right] 
+2\left[ \begin{array}{r}	0&	
	-q_0q_3&	
	q_0q_2\\	q_0q_3&
		0&
		-q_0q_1\\
	-q_0q_2&	
	q_0q_1&	
	0\\\end{array} \right] \\\bm{R}^{\mathrm{T}}\left( \bm{q} \right) 
=\left[ \begin{matrix}	\left( q_{0}^{2}+q_{1}^{2}-q_{2}^{2}-q_{3}^{2} \right)&	
	2q_1q_2&	
	2q_1q_3\\	2q_1q_2&
	
	\left( q_{0}^{2}-q_{1}^{2}+q_{2}^{2}-q_{3}^{2} \right)&	
	2q_2q_3\\	2q_1q_3&
		2q_2q_3&
	
	\left( q_{0}^{2}-q_{1}^{2}-q_{2}^{2}+q_{3}^{2} \right)\\\end{matrix} \right] 
+2\left[ \begin{array}{r}	0&	
	q_0q_3&	
	-q_0q_2\\	-q_0q_3&
		0&
		q_0q_1\\
	q_0q_2&	
	-q_0q_1&	
	0\\\end{array} \right] \\\bm{R}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{f}=\left[ \begin{matrix}
	\left( q_{0}^{2}+q_{1}^{2}-q_{2}^{2}-q_{3}^{2} \right)&	
	2q_1q_2&	
	2q_1q_3\\	2q_1q_2&
	
	\left( q_{0}^{2}-q_{1}^{2}+q_{2}^{2}-q_{3}^{2} \right)&	
	2q_2q_3\\	2q_1q_3&
		2q_2q_3&
	
	\left( q_{0}^{2}-q_{1}^{2}-q_{2}^{2}+q_{3}^{2} \right)\\\end{matrix} \right] 
\bm{f}+2\left[ \begin{array}{r}	0&	
	q_0q_3&	
	-q_0q_2\\	-q_0q_3&
		0&
		q_0q_1\\
	q_0q_2&	
	-q_0q_1&	
	0\\\end{array} \right] \bm{f}\\\frac{\partial}{\partial 
\bm{q}}q_0\left[ \begin{array}{r}	0&
		q_3&
		-q_2\\
	-q_3&	
	0&	
	q_1\\	q_2&
		-q_1&
		0\\\end{array} \right] 
\bm{f}=\frac{\partial}{\partial \bm{q}}q_0\left( \begin{array}{c}
	q_3f_2-q_2f_3\\	-q_3f_1+q_1f_3\\
	q_2f_1-q_1f_2\\\end{array} \right) =\left[ -2\begin{matrix}
	\bm{v}_{\times}\bm{f}&	
	2q_0\mathbf{I}\bm{f}_{\times}\\\end{matrix} \right] 
\\\frac{\partial}{\partial 
\bm{q}}\left( \bm{R}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{f} \right) =2\left( q_0\mathbf{I}-\bm{v}_{\times} \right) 
\bm{f}\\2\left( q_1\eta _1+q_2\eta _2+q_3\eta _3 \right) 
\mathbf{I}+2\left( \left( q_0\mathbf{I}-\bm{v}_{\times} \right) 
\bm{f} \right) _{\times}\\\frac{\partial 
\bm{R}\left( \bm{q} \right) 
\underline{\bm{\eta }}}{\partial \bm{q}}=\left[ \begin{array}{r}
	2q_0\eta _1-2q_3\eta _2+2q_2\eta _3&	
{2q_1\eta _1+2q_2\eta _2+2q_3\eta _3}&	
{-2q_2\eta _1+2q_1\eta _2+2q_0\eta _3}&	
	{2q_3\eta _1-2q_0\eta _2+2q_1\eta _3}\\	2q_3\eta 
_1+2q_0\eta _2-2q_1\eta _3&	{2q_2\eta 
_1-2q_1\eta _2-2q_0\eta _3}&	{2q_1\eta 
_1+2q_2\eta _2+2q_3\eta _3}&	{2q_0\eta 
_1-2q_3\eta _2+2q_2\eta _3}\\	-2q_2\eta _1+2q_1\eta _2+2q_0\eta _3&	
	{2q_3\eta _1+2q_0\eta _2-2q_1\eta _3}&	
	{-2q_0\eta _1+2q_3\eta _2-2q_2\eta _3}&	
	{2q_1\eta _1+2q_2\eta _2+2q_3\eta _3}\\\end{array} \right] 
\\2\left( q_0\mathbf{I}+\bm{v}_{\times} \right) 
\bm{\eta }\\2\left( q_1\eta _1+q_2\eta _2+q_3\eta _3 \right) 
\mathbf{I}-2\left( \left( q_0\mathbf{I}+\bm{v}_{\times} \right) 
\bm{\eta } \right) \\\frac{\partial}{\partial 
\bm{q}}\left( -2\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\hat{\bm{c}}^{\mathrm{T}}\bm{R}^{\mathrm{T}}\left( \bm{
q} \right) \bm{f} \right) 
=-2\left( \bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\hat{\bm{c}}^{\mathrm{T}}\frac{\partial}{\partial 
\bm{q}}\left( \bm{R}^{\mathrm{T}}\left( \bm{q} \right) 
\bm{f} \right) +\frac{\partial}{\partial 
\bm{q}}\left( \bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\underline{\hat{\bm{c}}^{\mathrm{T}}\bm{R}^{\mathrm{T}}\left( \bm{q} \right) \bm{f}} \right) \right) \\\frac{\partial 
\left( \bm{C}^{\mathrm{T}}\bm{f} \right)}{\partial 
\bm{q}}=\frac{\partial}{\partial \bm{q}}\left( \begin{array}{c}
	\bm{f}\\
	-2\bm{L}^{\mathrm{T}}\left( \bm{q} \right) 
\hat{\bm{c}}^{\mathrm{T}}\bm{R}^{\mathrm{T}}\left( \bm{
q} \right) \bm{f}\\\end{array} \right) \\=\left[ \begin{array}{c}
	\bm{0}\\	\\\end{array} \right] $$
	
Total kinetic energy
$$T_{\gamma}=T_{\mathrm{tra}}+T_{\mathrm{rot},\gamma}$$

Where translational kinetic energy
$$T_{\mathrm{tra}}=m\dot{\bm{r}}_{O}^{\mathrm{T}}\dot{\bm{r}}_O
$$

Total momentum
$$\bm{p}=\frac{\partial T_{\gamma}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}=\bm{M}_{\gamma}\dot{\bm{x}}$$

Where total mass matrix
$$\bm{M}_{\gamma}=\left[ \begin{matrix}	m\mathbf{I}_3&
		\bm{0}\\
	\bm{0}&	
	\bm{M}_{\mathrm{rot},\gamma}\\\end{matrix} \right] $$
	
Quadratic velocity term
$$\frac{\partial T_{\gamma}}{\partial 
\bm{x}^{\mathrm{T}}}=\left[ \begin{matrix}	\bm{0}&
		\bm{0}\\
	\bm{0}&	
	-4\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} \right) 
\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}\\\end{matrix} \right] $$

(not needed for variational integrator)

$$\dot{\bm{M}}_{\gamma}\dot{\bm{x}}=\left[ \begin{matrix}
	\bm{0}&	
	\bm{0}\\	\bm{0}&
	
	4\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} \right) 
\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}\\\end{matrix} \right] $$

$$\dot{\bm{M}}_{\gamma}\dot{\bm{x}}-\frac{\partial 
T_{\gamma}}{\partial \bm{x}^{\mathrm{T}}}=\left[ \begin{matrix}
	\bm{0}&	
	\bm{0}\\	\bm{0}&
	
	8\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} \right) 
\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}\\\end{matrix} \right] $$

### Dynamics
The Lagrangian 
$$L_{\gamma}=T_{\gamma}-V$$

Action
$$S=\int_{t_0}^{t_{\mathrm{f}}}{\left( L_{\gamma}\left( \bm{x},\dot{\bm{x}} \right) +\bm{\varPhi }^{\mathrm{T}}\bm{\lambda } 
\right) \mathrm{d}t}$$

Lagrange’s equations of the first kind

Differential-algebraic equations (not needed for variational integrator)

$$\left\{ \begin{array}{r}	\left[ \begin{matrix}
	m\mathbf{I}_3&	
	\bm{0}\\	\bm{0}&
	
	\bm{M}_{\mathrm{rot},\gamma}\\\end{matrix} \right] 
\ddot{\bm{x}}+\left[ \begin{matrix}	\bm{0}&
		\bm{0}\\
	\bm{0}&	
	8\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}} \right) 
\bm{J}_{\gamma}\bm{L}\left( \bm{q} \right) 
\dot{\bm{q}}\\\end{matrix} \right] +\frac{\partial V}{\partial 
\bm{x}^{\mathrm{T}}}-\bm{F}+\bm{A}^{\mathrm{T}}\bm{\lambda }=\bm{0}\\	\bm{\varPhi }\left( \bm{q} \right) 
=\bm{0}\\\end{array} \right. $$	

Discrete action
$$S_{\left( k-{{1}/{2}} \right)}=\left( L_{\gamma ,\left( k-{{1}/{2}} 
\right)}+\frac{1}{2}\left( \bm{\varPhi }^{\mathrm{T}}\left( \bm
{x}_{\left( k-1 \right)} \right) 
+\bm{\varPhi }^{\mathrm{T}}\left( \bm{x}_{\left( k \right)} 
\right) \right) \bm{\lambda }_{\left( k-{{1}/{2}} \right)} \right) h$$

where 
$$L_{\gamma ,\left( k-{{1}/{2}} \right)}\coloneqq 
L_{\gamma}\left( \bm{x}_{\left( k-{{1}/{2}} 
\right)},\dot{\bm{x}}_{\left( k-{{1}/{2}} \right)} \right) $$

With central finite differences
$$\dot{\bm{x}}_{\left( k-{{1}/{2}} \right)}\coloneqq 
\frac{\bm{x}_{\left( k \right)}-\bm{x}_{\left( k-1 
\right)}}{h}\\\bm{x}_{\left( k-{{1}/{2}} \right)}\coloneqq 
\frac{\bm{x}_{\left( k \right)}+\bm{x}_{\left( k-1 \right)}}{2}$$

Discrete Lagrange-d’Alembert principle
$$\delta S_{\left( k-{{1}/{2}} \right)}=\delta \bm{q}_{\left( k 
\right)}^{\mathrm{T}}\frac{\partial S_{\left( k-{{1}/{2}} 
\right)}}{\partial {\bm{x}^{\mathrm{T}}}_{\left( k \right)}}+\delta 
\bm{q}_{\left( k-1 \right)}^{\mathrm{T}}\frac{\partial 
S_{\left( k-{{1}/{2}} \right)}}{\partial 
{\bm{x}^{\mathrm{T}}}_{\left( k-1 \right)}}=0$$

$$\frac{\partial S_{\left( k-{{1}/{2}} \right)}}{\partial 
{\bm{x}^{\mathrm{T}}}_{\left( k-1 
\right)}}=\left( \frac{1}{2}\frac{\partial L_{\gamma ,\left( k-{{1}/{2}} 
\right)}}{\partial \bm{x}^{\mathrm{T}}}-\frac{1}{h}\frac{\partial 
L_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}+\frac{1}{2}\bm{A}^{\mathrm{T}}\left(
 \bm{x}_{\left( k-1 \right)} \right) 
\bm{\lambda }_{\left( k-{{1}/{2}} \right)} \right) h\\\frac{\partial 
S_{\left( k-{{1}/{2}} \right)}}{\partial 
{\bm{x}^{\mathrm{T}}}_{\left( k 
\right)}}=\left( \frac{1}{2}\frac{\partial L_{\gamma ,\left( k-{{1}/{2}} 
\right)}}{\partial \bm{x}^{\mathrm{T}}}+\frac{1}{h}\frac{\partial 
L_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}+\frac{1}{2}\bm{A}^{\mathrm{T}}\left(
 \bm{x}_{\left( k \right)} \right) 
\bm{\lambda }_{\left( k-{{1}/{2}} \right)} \right) h$$	

$$\frac{\partial L_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}=\bm{p}_{\left( k-{{1}/{2}} 
\right)}=\bm{M}_{\gamma}\left( \bm{x}_{\left( k-{{1}/{2}} 
\right)} \right) \dot{\bm{x}}_{\left( k-{{1}/{2}} 
\right)}\\\frac{\partial T_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\bm{x}^{\mathrm{T}}}=\left[ \begin{matrix}	\bm{0}&
		\bm{0}\\
	\bm{0}&	
	4\bm{L}^{\mathrm{T}}\left( \dot{\bm{q}}_{\left( k-{{1}/{2}} 
\right)} \right) 
\bm{J}_{\gamma}\bm{L}\left( \dot{\bm{q}}_{\left( k-{{1}
/{2}} \right)} \right) \bm{q}_{\left( k-{{1}/{2}} 
\right)}\\\end{matrix} \right] $$

Introduce discrete conjugate momentum
$$\bm{p}_{\left( k-1 
\right)}=\bm{M}_{\gamma}\left( \bm{x}_{\left( k-1 \right)} \right) 
\dot{\bm{x}}_{\left( k-1 \right)}\\\bm{p}_{\left( k 
\right)}=\bm{M}_{\gamma}\left( \bm{x}_{\left( k \right)} \right) 
\dot{\bm{x}}_{\left( k \right)}$$	

$$\delta \bm{q}_{k-1}^{\mathrm{T}}\bm{p}_{k-1}-\delta 
\bm{q}_{k}^{\mathrm{T}}\bm{p}_k=0$$

$$\left\{ \begin{aligned}	\bm{p}_{\left( k-1 
\right)}&=-\frac{\partial S_{\left( k-{{1}/{2}} \right)}}{\partial 
\bm{x}_{\left( k-1 
\right)}^{\mathrm{T}}}-\frac{h}{2}\bm{F}_{\left( k-{{1}/{2}} 
\right)}\\	\bm{p}_{\left( k \right)}&=\frac{\partial 
S_{\left( k-{{1}/{2}} \right)}}{\partial \bm{x}_{\left( k 
\right)}^{\mathrm{T}}}+\frac{h}{2}\bm{F}_{\left( k-{{1}/{2}} 
\right)}\\	\bm{\varPhi }\left( \bm{q}_{\left( k \right)} \right) 
&=\bm{0}\\\end{aligned} \right. $$	

$$\left\{ \begin{aligned}	\bm{p}_{\left( k-1 
\right)}&=-\left( \frac{h}{2}\frac{\partial L_{\gamma ,\left( k-{{1}/{2}} 
\right)}}{\partial \bm{x}^{\mathrm{T}}}-\frac{\partial 
L_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}+\frac{h}{2}\bm{A}^{\mathrm{T}}\left(
 \bm{x}_{\left( k-1 \right)} \right) 
\bm{\lambda }_{\left( k-{{1}/{2}} \right)} \right) 
-\frac{h}{2}\bm{F}_{\left( k-{{1}/{2}} \right)}\\
	\bm{p}_{\left( k \right)}&=\left( \frac{h}{2}\frac{\partial 
L_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\bm{x}^{\mathrm{T}}}+\frac{\partial L_{\gamma ,\left( k-{{1}/{2}} 
\right)}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}+\frac{h}{2}\bm{A}^{\mathrm{T}}\left(
 \bm{x}_{\left( k \right)} \right) 
\bm{\lambda }_{\left( k-{{1}/{2}} \right)} \right) 
+\frac{h}{2}\bm{F}_{\left( k-{{1}/{2}} \right)}\\
	\bm{\varPhi }\left( \bm{q}_{\left( k \right)} \right) 
&=\bm{0}\\\end{aligned} \right. $$

$$\bm{p}_{\left( k-1 \right)}+\bm{p}_{\left( k 
\right)}=2\frac{\partial L_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}+\frac{h}{2}\left( \bm{A}^{\mathrm{T}
}\left( \bm{x}_{\left( k \right)} \right) 
-\bm{A}^{\mathrm{T}}\left( \bm{x}_{\left( k-1 \right)} \right) 
\right) \bm{\lambda }_{\left( k-{{1}/{2}} 
\right)}\\=2\bm{M}_{\gamma}\left( \bm{x}_{\left( k-{{1}/{2}} 
\right)} \right) \dot{\bm{x}}_{\left( k-{{1}/{2}} 
\right)}+\frac{h}{2}\left( \bm{A}^{\mathrm{T}}\left( \bm{x}_{\left( k \right)} \right) 
-\bm{A}^{\mathrm{T}}\left( \bm{x}_{\left( k-1 \right)} \right) 
\right) \bm{\lambda }_{\left( k-{{1}/{2}} \right)}$$

$$\begin{aligned}
	\bm{M}_{\gamma ,\left( k-1 \right)}^{-1}\bm{p}_{\left( k-1 
\right)}&=-\bm{M}_{\gamma ,\left( k-1 
\right)}^{-1}\left( \frac{h}{2}\frac{\partial L_{\gamma ,\left( k-{{1}/{2}} 
\right)}}{\partial \bm{x}^{\mathrm{T}}}+\frac{\partial 
L_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}+\frac{h}{2}\bm{A}^{\mathrm{T}}\left(
 \bm{x}_{\left( k \right)} \right) 
\bm{\lambda }_{\left( k-{{1}/{2}} \right)} \right) 
-\frac{h}{2}\bm{M}_{\gamma ,\left( k-1 
\right)}^{-1}\bm{F}_{\left( k-{{1}/{2}} \right)}\\
	\bm{M}_{\gamma ,\left( k \right)}^{-1}\bm{p}_{\left( k 
\right)}&=\bm{M}_{\gamma ,\left( k 
\right)}^{-1}\left( \frac{h}{2}\frac{\partial L_{\gamma ,\left( k-{{1}/{2}} 
\right)}}{\partial \bm{x}^{\mathrm{T}}}-\frac{\partial 
L_{\gamma ,\left( k-{{1}/{2}} \right)}}{\partial 
\dot{\bm{x}}^{\mathrm{T}}}+\frac{h}{2}\bm{A}^{\mathrm{T}}\left(
 \bm{x}_{\left( k-1 \right)} \right) 
\bm{\lambda }_{\left( k-{{1}/{2}} \right)} \right) 
+\frac{h}{2}\bm{M}_{\gamma ,\left( k 
\right)}^{-1}\bm{F}_{\left( k-{{1}/{2}} \right)}\\\end{aligned}$$

$$\left\{ \begin{array}{r}
	\bm{p}_{\left( k-1/2 \right)}-\bm{p}_{\left( k-1 
\right)}-\tfrac{h}{2}\bm{F}_{\left( k-1/2 
\right)}-\tfrac{h}{2}\bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 \right)}=0\\
	\bm{p}_{\left( k \right)}-\bm{p}_{\left( k-1/2 
\right)}-\tfrac{h}{2}\bm{F}_{\left( k-1/2 
\right)}-\tfrac{h}{2}\bm{A}\left( \bm{q}_{\left( k \right)} 
\right) ^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 \right)}=0\\\end{array} 
\right. \\\bm{p}_{\left( k \right)}+\bm{p}_{\left( k-1 
\right)}=2\bm{p}_{\left( k-1/2 
\right)}+\tfrac{h}{2}\left( \bm{A}\left( \bm{q}_{\left( k \right)} 
\right) ^{\mathrm{T}}-\bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}} \right) \bm{\lambda }_{\left( k-1/2 
\right)}\\\bm{p}_{\left( k-1/2 
\right)}=\bm{M}\left( \bm{q}_{k-1/2} \right) 
\dot{\bm{q}}_{\left( k-1/2 \right)}\\\bm{F}_{\left( k-1/2 
\right)}=\frac{\partial T_{\gamma ,\left( k-1/2 \right)}}{\partial 
\bm{x}^{\mathrm{T}}}-\frac{\partial V_{\left( k-1/2 \right)}}{\partial 
\bm{x}^{\mathrm{T}}}+\bm{F}_{\left( k-1/2 
\right)}^{\mathrm{ex}}$$	

MSI
$$\bm{M}\mathrm{d}\bm{v}-\bm{F}\mathrm{d}t-\bm{
A}^{\mathrm{T}}\mathrm{d}\bm{\lambda }_b-\bm{D}^{\mathrm{T}}\mathrm{d}\bm{\lambda }_u=0\\\bm{p}_{\left( k 
\right)}-\bm{p}_{\left( k-1/2 \right)}+\bm{p}_{\left( k-1/2 
\right)}-\bm{p}_{\left( k-1 
\right)}-h\bm{F}\left( \frac{\bm{q}_{\left( k 
\right)}+\bm{q}_{\left( k-1 \right)}}{2},\frac{\bm{q}_{\left( k 
\right)}-\bm{q}_{\left( k-1 \right)}}{h},t_{\left( k-{{1}/{2}} 
\right)} \right) 
-\tfrac{h}{2}\left( \bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}}+\bm{A}\left( \bm{q}_{\left( k \right)} 
\right) ^{\mathrm{T}} \right) \bm{\lambda }_{\left( k-1/2 
\right)}-\tfrac{h}{2}\left( \bm{D}\left( \bm{q}_{\left( k-1 
\right)} \right) ^{\mathrm{T}}+\bm{D}\left( \bm{q}_{\left( k 
\right)} \right) ^{\mathrm{T}} \right) \bm{\varLambda }_{\left( k-1/2 
\right)}=0\\\bm{p}_{\left( k \right)}-\bm{p}_{\left( k-1/2 
\right)}+\bm{p}_{\left( k-1/2 \right)}-\bm{p}_{\left( k-1 
\right)}-h\bm{F}\left( \frac{\bm{q}_{\left( k 
\right)}+\bm{q}_{\left( k-1 \right)}}{2},\frac{\bm{q}_{\left( k 
\right)}-\bm{q}_{\left( k-1 \right)}}{h},t_{\left( k-{{1}/{2}} 
\right)} \right) 
-\tfrac{h}{2}\left( \bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}}+\bm{A}\left( \bm{q}_{\left( k \right)} 
\right) ^{\mathrm{T}} \right) \bm{\lambda }_{\left( k-1/2 
\right)}=0\\\bm{p}_{\left( k-1/2 
\right)}^{+}-h\bm{F}_{\left( k-{{1}/{2}} 
\right)}-\tfrac{h}{2}\left( \bm{A}\left( \bm{q}_{\left( k-{{1}\
Bigg/{2}} \right)} \right) ^{\mathrm{T}} \right) 
\bm{\lambda }_{\left( k-{{1}/{2}} 
\right)}=0\\\bm{q}_{\left( k \right)}\rightarrow 
\bm{q}_{\left( k-1/2 \right)}\\\bm{p}_{\left( k 
\right)}-\bm{p}_{\left( k-1/2 
\right)}-\frac{h}{2}\bm{F}_{\left( k-1/2 
\right)}-\frac{h}{2}\bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 
\right)}-\frac{h}{2}\bm{D}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}}\bm{\varLambda }_{\left( k-1/2 
\right)}=0\\\bm{p}_{\left( k-1/2 \right)}-\bm{p}_{\left( k-1 
\right)}-\frac{h}{2}\bm{F}_{\left( k-1/2 
\right)}-\frac{h}{2}\bm{A}\left( \bm{q}_{\left( k \right)} \right) 
^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 
\right)}-\frac{h}{2}\bm{D}\left( \bm{q}_{\left( k \right)} \right) 
^{\mathrm{T}}\bm{\varLambda }_{\left( k-1/2 
\right)}=0\\\\\\\\\bm{p}_{\left( k \right)}+\bm{p}_{\left( k-1 
\right)}=2\bm{p}_{\left( k-1/2 
\right)}+\tfrac{h}{2}\left( \bm{A}\left( \bm{q}_{\left( k \right)} 
\right) ^{\mathrm{T}}-\bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}} \right) \bm{\lambda }_{\left( k-1/2 
\right)}\\\bm{p}_{\left( k-1/2 
\right)}=\bm{M}\left( \bm{q}_{k-1/2} \right) 
\dot{\bm{q}}_{\left( k-1/2 \right)}\\\bm{F}_{\left( k-1/2 
\right)}=\frac{\partial T_{\gamma ,\left( k-1/2 \right)}}{\partial 
\bm{x}^{\mathrm{T}}}-\frac{\partial V_{\left( k-1/2 \right)}}{\partial 
\bm{x}^{\mathrm{T}}}+\bm{F}_{\left( k-1/2 
\right)}^{\mathrm{ex}}$$	


Integration scheme
$$\left\{ \begin{array}{r}
	h\bm{M}\left( \bm{q}_{k-1/2} \right) 
\dot{\bm{q}}_{\left( k-1/2 \right)}-h\bm{p}_{\left( k-1 
\right)}-\tfrac{h^2}{2}\bm{F}_{\left( k-1/2 
\right)}-\tfrac{h^2}{2}\bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 \right)}=0\\
	\bm{\varPhi }\left( \bm{q}_{\left( k \right)} \right) 
=0\\\end{array} \right. \\\bm{p}_{\left( k 
\right)}=-\bm{p}_{\left( k-1 \right)}+2\bm{p}_{\left( k-1/2 
\right)}+\tfrac{h}{2}\left( \bm{A}\left( \bm{q}_{\left( k \right)} 
\right) ^{\mathrm{T}}-\bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}} \right) \bm{\lambda }_{\left( k-1/2 
\right)}\\\dot{\bm{q}}_{\left( k-1/2 
\right)}=\frac{\bm{q}_{\left( k \right)}-\bm{q}_{\left( k-1 
\right)}}{h},\quad \frac{\partial \dot{\bm{q}}_{\left( k-1/2 
\right)}}{\partial \bm{q}_{\left( k 
\right)}}=\frac{1}{h}\mathbf{I}\\\bm{q}_{\left( k-1/2 
\right)}=\frac{\bm{q}_{\left( k \right)}+\bm{q}_{\left( k-1 
\right)}}{2},\quad \frac{\partial \bm{q}_{\left( k-1/2 \right)}}{\partial 
\bm{q}_{\left( k \right)}}=\frac{1}{2}\mathbf{I}\\\frac{\partial 
\left( h\bm{M}\left( \bm{q}_{k-1/2} \right) 
\dot{\bm{q}}_{\left( k-1/2 \right)} \right)}{\partial 
\bm{q}_{\left( k \right)}}=\frac{\partial 
\left( \bm{M}\left( \bm{q}_{k-1/2} \right) 
\left( \bm{q}_{\left( k \right)}-\bm{q}_{\left( k-1 \right)} 
\right) \right)}{\partial \bm{q}_{\left( k \right)}}\\=\frac{\partial 
\left( \bm{M}\left( \bm{q}_{k-1/2} \right) 
\bm{q}_{\left( k \right)} \right)}{\partial \bm{q}_{\left( k 
\right)}}\\=\bm{M}\left( \bm{q}_{k-1/2} \right) +\frac{\partial 
\left( \bm{M}\left( \bm{q}_{k-1/2} \right) {\bm{q}_{\left( k \right)}} \right)}{\partial 
\bm{q}_{\left( k \right)}}\\\\\bm{F}_{\left( k-1/2 
\right)}=\frac{\partial T_{\gamma ,\left( k-1/2 \right)}}{\partial 
\bm{x}^{\mathrm{T}}}-\frac{\partial V_{\left( k-1/2 \right)}}{\partial 
\bm{x}^{\mathrm{T}}}+\bm{F}_{\left( k-1/2 
\right)}^{\mathrm{ex}}$$	

### Contact 
Denote the $i$th contact point’s post-contact velocity by $\dot{\bm{r}}_{i}^{+}$, then 
the post-contact relative velocity reads
$$\acute{\bm{v}}_{i}^{+}=\left[ \begin{array}{c}
	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] 
\dot{\bm{r}}_{i}^{+}$$

The $i$th contact point is distinguished as impact or persistent.

a)	Continuous Newton’s impact law 
$$\acute{v}_{i,n,\mathrm{imp}}^{+}+e_i\acute{v}_{i,n,\mathrm{imp}}^{-}\ge 0$$

b)	Continuous non-compliance at velocity level
$$\acute{v}_{i,n,\mathrm{per}}^{+}\ge 0$$	

Let
$$\hat{\bm{\nu}}_i=\acute{\bm{v}}_{i}^{+}+\bm{b}_i$$

where 

$$\bm{b}_i=\left( \begin{array}{c}	w_i+\acute{v}_{i,t}^{+}\\	0\\
	0\\\end{array} \right) $$
, and 
$w_i=\begin{cases}	e_i\acute{v}_{i,n,\mathrm{imp}}^{-}\\
	0\\\end{cases}$
Cone complementarity condition for a frictional contact point reads 

$$\mathcal{K} \ni \hat{\bm{\varLambda}}_i\bot \hat{\mathbf{\nu}}_i\in 
\mathcal{K} $$

Time-stepping process
(2)	Predict position $\bm{q}_{\left( k-{{1}/{2}} \right)}^{*}\coloneqq 
\bm{q}_{\left( k-1 \right)}+\frac{h}{2}\dot{\bm{q}}_{\left( k-1 \right)}$ 
to find contact points’ local coordinates $\bm{c}_{I,i}$and contact directions 
$$\left[ \begin{array}{c}	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} 
\right] $$

 both of which are held constant during timestep $k$.
(3)	Discrete Newton’s impact law 
$$\acute{\bm{v}}_{i,\left( k \right) ,\mathrm{imp}}\coloneqq 
\left[ \begin{array}{c}	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] 
\dot{\bm{r}}_{i,\left( k \right)}\\\bm{b}_{i,\left( k 
\right) ,\mathrm{imp}}\coloneqq \left( \begin{array}{c}
	e_i\bm{n}_{i}^{\mathrm{T}}\dot{\bm{r}}_{i,\left( k-1 
\right)}+\acute{v}_{i,t,\left( k \right) ,\mathrm{imp}}\\	0\\
	0\\\end{array} \right) $$
	
$$\mathcal{K} \ni \hat{\bm{\varLambda}}_{i,\mathrm{imp}}\bot 
\acute{\bm{v}}_{i,\left( k 
\right) ,\mathrm{imp}}+\bm{b}_{i,\left( k \right) ,\mathrm{imp}}\in 
\mathcal{K} $$	

(4)	Discrete non-compliance at velocity level 
$$\acute{\bm{v}}_{i,\left( k \right) ,\mathrm{per}}\coloneqq 
\left[ \begin{array}{c}	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] 
\dot{\bm{r}}_{i,\left( k-{{1}/{2}} \right)}=\left[ \begin{array}{c}
	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] 
\frac{1}{h}\left( \bm{r}_{\left( k \right)}-\bm{r}_{\left( k-1 
\right)} \right) \\\bm{b}_{i,\left( k \right) ,\mathrm{per}}\coloneqq 
\left( \begin{array}{c}	\acute{v}_{i,t,\left( k-{{1}/{2}} 
\right) ,\mathrm{imp}}\\	0\\	0\\\end{array} \right) $$

$$\mathcal{K} \ni \hat{\bm{\varLambda}}_{i,\mathrm{per}}\bot 
\acute{\bm{v}}_{i,\left( k 
\right) ,\mathrm{per}}+\bm{b}_{i,\left( k \right) ,\mathrm{per}}\in 
\mathcal{K} $$

(5)	Central approximations of the generalized directions
$$\bm{D}_{i,\left( k-{{1}/{2}} 
\right)}\,\,=\frac{1}{2}\left( \bm{D}_{i,\left( k-1 
\right)}+\bm{D}_{i,\left( k \right)} \right) 
\\=\frac{1}{2}\left[ \begin{array}{c}	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] 
\left( \bm{C}_i\left( \bm{q}_{\left( k-1 
\right)},\bm{c}_{I,i} \right) 
+\bm{C}_i\left( \bm{q}_{\left( k \right)},\bm{c}_{I,i} 
\right) \right) $$	

(6)	To be consistent with the symplectic scheme, the generalized contact force has the form
$$\bm{D}_{i,\left( k-1 
\right)}^{\mathrm{T}}\bm{\varLambda }_{i,\left( k-{{1}/{2}} 
\right)}$$	

(7)	Jacobian of the nonholonomic constraint without friction 
$$\bm{D}_{i,n,\left( k \right) ,\mathrm{imp}}=\frac{\partial 
\left( \acute{\bm{v}}_{i,n,\left( k \right) ,\mathrm{imp}} 
\right)}{\partial \dot{\bm{q}}_{\left( k 
\right)}}=\left[ \bm{n}_{i}^{\mathrm{T}} \right] \frac{\partial 
\left( \dot{\bm{r}}_{i,\left( k \right)} \right)}{\partial 
\dot{\bm{q}}_{\left( k \right)}}=\left[ \bm{n}_{i}^{\mathrm{T}} 
\right] \bm{C}_i\left( \bm{q}_{\left( k 
\right)},\bm{c}_{I,i} \right) $$

(8)	Jacobian of the holonomic constraint without friction
$$\bm{D}_{i,n,\left( k \right) ,\mathrm{per}}=\frac{\partial 
\left( \acute{\bm{v}}_{i,n,\left( k \right) ,\mathrm{per}} 
\right)}{\partial \bm{q}_{\left( k 
\right)}}=\frac{1}{h}\left[ \bm{n}_{i}^{\mathrm{T}} \right] 
\frac{\partial \left( \bm{r}_{i,\left( k \right)} \right)}{\partial 
\bm{q}_{\left( k 
\right)}}=\frac{1}{h}\left[ \bm{n}_{i}^{\mathrm{T}} \right] 
\bm{C}_i\left( \bm{q}_{\left( k \right)},\bm{c}_{I,i} 
\right) $$

Non-smooth integration scheme
$$\left\{ \begin{array}{r}	h\bm{M}\left( \bm{q}_{k-1/2} \right) 
\dot{\bm{q}}_{\left( k-1/2 \right)}-h\bm{p}_{\left( k-1 
\right)}-\tfrac{h^2}{2}\bm{F}_{\left( k-1/2 
\right)}-\bm{A}\left( \bm{q}_{\left( k-1 \right)} \right) 
^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 
\right)}-\bm{D}_{\left( k-1 
\right)}^{\mathrm{T}}\bm{H}_{\left( k 
\right)}\hat{\bm{\varLambda}}_{\left( k-{{1}/{2}} \right)}=0\\
	\bm{\varPhi }\left( \bm{q}_{\left( k \right)} \right) =0\\
	\hat{\mathbf{\nu}}_{\left( k 
\right) ,\mathrm{imp}}-\left( \acute{\bm{v}}_{i,\left( k 
\right) ,\mathrm{imp}}+\bm{b}_{i,\left( k \right) ,\mathrm{imp}} \right) 
=0\\	\hat{\mathbf{\nu}}_{\left( k 
\right) ,\mathrm{per}}-\left( \acute{\bm{v}}_{i,\left( k 
\right) ,\mathrm{per}}+\bm{b}_{i,\left( k \right) ,\mathrm{per}} \right) 
=0\\	\mathcal{C} \ni \hat{\bm{\varLambda}}_{\left( k-{{1}/{2}} 
\right)}\bot \hat{\mathbf{\nu}}_{\left( k \right)}\in \mathcal{C}\\\end{array} 
\right. \\\bm{p}_{\left( k \right)}=-\bm{p}_{\left( k-1 
\right)}+2\bm{M}\left( \bm{q}_{k-1/2} \right) 
\dot{\bm{q}}_{\left( k-1/2 
\right)}+\frac{1}{h}\left( \bm{A}\left( \bm{q}_{\left( k \right)} 
\right) ^{\mathrm{T}}-\bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}} \right) \bm{\lambda }_{\left( k-1/2 
\right)}+\frac{1}{h}\left( \bm{D}_{\left( k 
\right)}^{\mathrm{T}}-\bm{D}_{\left( k-1 \right)}^{\mathrm{T}} \right) 
\bm{H}_{\left( k 
\right)}\hat{\bm{\varLambda}}_{\left( k-{{1}/{2}} 
\right)}\\\frac{\partial \bm{q}_{k-1/2}}{\partial 
\bm{q}_{\left( k \right)}}=\frac{1}{2}\mathbf{I}\\\frac{\partial 
\dot{\bm{q}}_{\left( k-1/2 \right)}}{\partial \bm{q}_{\left( k 
\right)}}=\frac{1}{h}\mathbf{I}\\\frac{\partial 
\bm{M}\left( \bm{q}_{k-1/2} \right) 
\underline{\dot{\bm{q}}_{\left( k-1/2 \right)}}}{\partial 
\bm{q}_{\left( k \right)}}=\frac{\partial 
\bm{M}\left( \bm{q}_{k-1/2} \right) 
\underline{\dot{\bm{q}}_{\left( k-1/2 \right)}}}{\partial 
\bm{q}_{\left( k-{{1}/{2}} \right)}}\frac{\partial 
\bm{q}_{\left( k-{{1}/{2}} \right)}}{\partial 
\bm{q}_{\left( k \right)}}=\frac{1}{2}\frac{\partial 
\bm{M}\left( \bm{q}_{k-1/2} \right) 
\underline{\dot{\bm{q}}_{\left( k-1/2 \right)}}}{\partial 
\bm{q}_{\left( k-{{1}/{2}} \right)}}\\\frac{\partial 
\bm{p}_{\left( k \right)}}{\partial \bm{q}_{\left( k 
\right)}}=2\bm{M}\left( \bm{q}_{k-1/2} \right) \frac{\partial 
\dot{\bm{q}}_{\left( k-1/2 \right)}}{\partial \bm{q}_{\left( k 
\right)}}+2\frac{\partial \bm{M}\left( \bm{q}_{k-1/2} \right) 
\underline{\dot{\bm{q}}_{\left( k-1/2 \right)}}}{\partial 
\bm{q}_{\left( k \right)}}+\frac{1}{h}\frac{\partial 
\bm{A}\left( \bm{q}_{\left( k \right)} \right) 
^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 \right)}}{\partial 
\bm{q}_{\left( k \right)}}+\frac{1}{h}\frac{\partial 
\bm{D}_{\left( k \right)}^{\mathrm{T}}\bm{H}_{\left( k 
\right)}\hat{\bm{\varLambda}}_{\left( k-{{1}/{2}} 
\right)}}{\partial \bm{q}_{\left( k 
\right)}}\\=\frac{2}{h}\bm{M}\left( \bm{q}_{k-1/2} \right) 
+\frac{\partial \bm{M}\left( \bm{q}_{k-1/2} \right) 
\underline{\dot{\bm{q}}_{\left( k-1/2 \right)}}}{\partial 
\bm{q}_{\left( k-{{1}/{2}} \right)}}+\frac{1}{h}\frac{\partial 
\bm{A}\left( \bm{q}_{\left( k \right)} \right) 
^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 \right)}}{\partial 
\bm{q}_{\left( k \right)}}+\frac{1}{h}\frac{\partial 
\bm{D}_{\left( k \right)}^{\mathrm{T}}\bm{H}_{\left( k 
\right)}\hat{\bm{\varLambda}}_{\left( k-{{1}/{2}} 
\right)}}{\partial \bm{q}_{\left( k 
\right)}}\\\bm{M}^{-1}\left( \bm{q}_{k-1/2} \right) 
\frac{\partial \bm{p}_{\left( k \right)}}{\partial 
\bm{q}_{\left( k 
\right)}}=\frac{2}{h}\mathbf{I}+\bm{M}^{-1}\left( \bm{q}_{k-1/2} 
\right) \left( \frac{\partial \bm{M}\left( \bm{q}_{k-1/2} \right) 
\underline{\dot{\bm{q}}_{\left( k-1/2 \right)}}}{\partial 
\bm{q}_{\left( k-{{1}/{2}} \right)}}+\frac{1}{h}\frac{\partial 
\bm{A}\left( \bm{q}_{\left( k \right)} \right) 
^{\mathrm{T}}\bm{\lambda }_{\left( k-1/2 \right)}}{\partial 
\bm{q}_{\left( k \right)}}+\frac{1}{h}\frac{\partial 
\bm{D}_{\left( k \right)}^{\mathrm{T}}\bm{H}_{\left( k 
\right)}\hat{\bm{\varLambda}}_{\left( k-{{1}/{2}} 
\right)}}{\partial \bm{q}_{\left( k \right)}} \right) 
\\\dot{\bm{q}}_{\left( k 
\right)}=\bm{M}^{-1}\left( \bm{q}_{k-1/2} \right) 
\bm{p}_{\left( k \right)}\\\frac{\partial \dot{\bm{q}}_{\left( k 
\right)}}{\partial \bm{q}_{\left( k 
\right)}}=\bm{M}^{-1}\left( \bm{q}_{k-1/2} \right) 
\frac{\partial \bm{p}_{\left( k \right)}}{\partial 
\bm{q}_{\left( k \right)}}+\frac{\partial 
\bm{M}^{-1}\left( \bm{q}_{k-1/2} \right) 
\underline{\bm{p}_{\left( k \right)}}}{\partial \bm{q}_{\left( k 
\right)}}\\\frac{\partial \dot{\bm{q}}_{\left( k \right)}}{\partial 
\bm{\lambda }_{\left( k-1/2 
\right)}}=\bm{M}^{-1}\left( \bm{q}_{k-1/2} \right) 
\frac{\partial \bm{p}_{\left( k \right)}}{\partial 
\bm{\lambda }_{\left( k-1/2 
\right)}}\\=\bm{M}^{-1}\left( \bm{q}_{k-1/2} \right) 
\left( \frac{1}{h}\left( \bm{A}\left( \bm{q}_{\left( k \right)} 
\right) ^{\mathrm{T}}-\bm{A}\left( \bm{q}_{\left( k-1 \right)} 
\right) ^{\mathrm{T}} \right) \right) $$	

$$\frac{\partial \bm{C}\left( \bm{q}_{I,\left( k 
\right)};\bm{c}_{I,i} \right) 
\underline{\dot{\bm{q}}_{I,\left( k \right)}}}{\partial 
\bm{q}_{\left( k \right)}}=\frac{\partial 
\bm{C}\left( \bm{q}_{I,\left( k \right)};\bm{c}_{I,i} 
\right) \underline{\dot{\bm{q}}_{I,\left( k \right)}}}{\partial 
\bm{q}_{I,\left( k \right)}}\frac{\partial \bm{q}_{I,\left( k 
\right)}}{\partial \bm{q}_{\left( k \right)}}\\=\frac{\partial 
\bm{C}\left( \bm{q}_{I,\left( k \right)};\bm{c}_{I,i} 
\right) \underline{\bm{T}_I\dot{\bm{q}}_{\left( k 
\right)}}}{\partial \bm{q}_{I,\left( k 
\right)}}\bm{T}_I\\\frac{\partial \dot{\bm{r}}_{i,\left( k 
\right)}}{\partial \bm{q}_{\left( k \right)}}=\frac{\partial 
\bm{C}\left( \bm{q}_{\left( k \right)};\bm{c}_{I,i} 
\right) \dot{\bm{q}}_{\left( k \right)}}{\partial 
\bm{q}_{\left( k 
\right)}}=\bm{C}_I\left( \bm{q}_{\left( k 
\right)};\bm{c}_{I,i} \right) \frac{\partial 
\dot{\bm{q}}_{\left( k \right)}}{\partial \bm{q}_{\left( k 
\right)}}+\frac{\partial \bm{C}\left( \bm{q}_{\left( k 
\right)};\bm{c}_{I,i} \right) \underline{\dot{\bm{q}}_{\left( k 
\right)}}}{\partial \bm{q}_{\left( k \right)}}\\\\\frac{\partial 
\acute{\bm{v}}_{i,\left( k \right) ,\mathrm{imp}}}{\partial 
\bm{q}_{\left( k \right)}}=\left[ \begin{array}{c}
	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] \frac{\partial 
\dot{\bm{r}}_{i,\left( k \right)}}{\partial \bm{q}_{\left( k 
\right)}}\\=\left[ \begin{array}{c}	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] 
\bm{C}_I\left( \bm{q}_{\left( k \right)};\bm{c}_{I,i} 
\right) \frac{\partial \dot{\bm{q}}_{\left( k \right)}}{\partial 
\bm{q}_{\left( k \right)}}+\left[ \begin{array}{c}
	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] \frac{\partial 
\bm{C}\left( \bm{q}_{\left( k \right)};\bm{c}_{I,i} 
\right) \underline{\dot{\bm{q}}_{\left( k \right)}}}{\partial 
\bm{q}_{\left( k \right)}}\\=\bm{D}_i\frac{\partial 
\dot{\bm{q}}_{\left( k \right)}}{\partial \bm{q}_{\left( k 
\right)}}+\frac{\partial 
\bm{D}_i\underline{\dot{\bm{q}}_{\left( k \right)}}}{\partial 
\bm{q}_{\left( k \right)}}\\\frac{\partial \hat{\mathbf{\nu}}_{\left( k 
\right)}}{\partial \bm{q}_{\left( k \right)}}=\frac{\partial}{\partial 
\bm{q}_{\left( k \right)}}\left( \acute{\bm{v}}_{i,\left( k 
\right) ,\mathrm{imp}}+\bm{b}_{i,\left( k \right) ,\mathrm{imp}} \right) 
\\\\\frac{\partial \acute{\bm{v}}_{i,\left( k 
\right) ,\mathrm{per}}}{\partial \bm{q}_{\left( k 
\right)}}=\frac{1}{h}\left[ \begin{array}{c}	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] \frac{\partial 
\bm{r}_{i,\left( k \right)}}{\partial \bm{q}_{\left( k 
\right)}}\\=\frac{1}{h}\bm{D}_i$$	

$$\bm{D}_{i,\left( k \right)}=\left[ \begin{array}{c}
	\bm{n}_{i}^{\mathrm{T}}\\
	\bm{t}_{i,1}^{\mathrm{T}}\\
	\bm{t}_{i,2}^{\mathrm{T}}\\\end{array} \right] 
\bm{C}_i\left( \bm{q}_{I,\left( k \right)},\bm{c}_{I,i} 
\right) \bm{T}_I\\\bm{f}_i=\left[ \begin{matrix}
	\bm{n}_i&	
	\bm{t}_{i,1}&	
	\bm{t}_{i,2}\\\end{matrix} \right] 
\bm{\varLambda }_i\\\bm{D}_{i,\left( k 
\right)}^{\mathrm{T}}\bm{\varLambda }_i=\bm{T}_{I}^{\mathrm{T}}
\bm{C}_{i}^{\mathrm{T}}\left( \bm{q}_{I,\left( k 
\right)},\bm{c}_{I,i} \right) \left[ \begin{matrix}	\bm{n}_i&
		\bm{t}_{i,1}&
		\bm{t}_{i,2}\\\end{matrix} 
\right] 
\bm{\varLambda }_i\\=\bm{T}_{I}^{\mathrm{T}}\bm{C}_{i}^
{\mathrm{T}}\left( \bm{q}_{I,\left( k \right)},\bm{c}_{I,i} 
\right) \bm{f}_i\\\frac{\partial \left( \bm{D}_{i,\left( k 
\right)}^{\mathrm{T}}\bm{\varLambda }_i \right)}{\partial 
\bm{q}_{I,\left( k \right)}}=$$


