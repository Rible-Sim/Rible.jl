$$
\begin{aligned}
	&\frac{\partial \kappa}{\partial e_i}=\frac{\partial}{\partial e_i}\frac{\left| \mathbf{r}_x\times \mathbf{r}_{xx} \right|}{\left| \mathbf{r}_x \right|^3}=\frac{1}{g^2}\left( g\frac{\partial f}{\partial e_i}-f\frac{\partial g}{\partial e_i} \right)\\
	&f=\left| \mathbf{r}_x\times \mathbf{r}_{xx} \right|,\quad g=\left| \mathbf{r}_x \right|^3\\
	&\frac{\partial f}{\partial e_i}=\frac{\partial}{\partial e_i}\left| \mathbf{r}_x\times \mathbf{r}_{xx} \right|=\frac{\partial}{\partial e_i}\left( \left( \mathbf{r}_x\times \mathbf{r}_{xx} \right) ^{\mathrm{T}}\left( \mathbf{r}_x\times \mathbf{r}_{xx} \right) \right) ^{1/2}\\
	&=\left( \left( \mathbf{r}_x\times \mathbf{r}_{xx} \right) ^{\mathrm{T}}\left( \mathbf{r}_x\times \mathbf{r}_{xx} \right) \right) ^{-1/2}\left( \left( \mathbf{r}_x\times \mathbf{r}_{xx} \right) ^{\mathrm{T}}\left( \frac{\partial}{\partial e_i}\mathbf{r}_x\times \mathbf{r}_{xx}+\mathbf{r}_x\times \frac{\partial}{\partial e_i}\mathbf{r}_{xx} \right) \right)\\
	&\frac{\partial g}{\partial e_i}=\frac{\partial}{\partial e_i}\left( \mathbf{r}_{x}^{\mathrm{T}}\mathbf{r}_x \right) ^{3/2}=3\left( \mathbf{r}_{x}^{\mathrm{T}}\mathbf{r}_x \right) ^{1/2}\left( \mathbf{r}_{x}^{\mathrm{T}}\frac{\partial \mathbf{r}_x}{\partial e_i} \right)\\
\end{aligned}
\\
\mathrm{\delta}\kappa =\frac{\partial \kappa}{\partial \boldsymbol{e}}\mathrm{\delta}\boldsymbol{e}
\\
\mathrm{\delta}\varepsilon _{xx}^{\mathrm{a}}=\frac{\partial \varepsilon _{xx}^{\mathrm{a}}}{\partial \boldsymbol{e}}\mathrm{\delta}\boldsymbol{e}=\boldsymbol{r}_{x}^{\mathrm{T}}\boldsymbol{S}_x\mathrm{\delta}\boldsymbol{e}
\\
\mathrm{\delta}W=\int_0^L{EI\kappa \mathrm{\delta}\kappa +EA\varepsilon _{xx}^{\mathrm{a}}\mathrm{\delta}\varepsilon _{xx}^{\mathrm{a}}\mathrm{d}x}
\\
\mathrm{\delta}W_I=EI\int_0^L{\kappa \mathrm{\delta}\kappa \mathrm{d}x}=EI\int_0^L{\kappa \frac{\partial \kappa}{\partial \boldsymbol{e}}\mathrm{\delta}\boldsymbol{e}\mathrm{d}x}=EI\left( \int_0^L{\kappa \frac{\partial \kappa}{\partial \boldsymbol{e}}\mathrm{d}x} \right) \mathrm{\delta}\boldsymbol{e}, \boldsymbol{Q}_I=EI\left( \int_0^L{\kappa \frac{\partial \kappa}{\partial \boldsymbol{e}}\mathrm{d}x} \right) 
\\
\mathrm{\delta}W_A=EA\int_0^L{\varepsilon _{xx}^{\mathrm{a}}\mathrm{\delta}\varepsilon _{xx}^{\mathrm{a}}\mathrm{d}x}=EA\int_0^L{\varepsilon _{xx}^{\mathrm{a}}\boldsymbol{r}_{x}^{\mathrm{T}}\boldsymbol{S}_x\mathrm{\delta}\boldsymbol{e}\mathrm{d}x}=EA\left( \int_0^L{\varepsilon _{xx}^{\mathrm{a}}\boldsymbol{r}_{x}^{\mathrm{T}}\boldsymbol{S}_x\mathrm{d}x} \right) \mathrm{\delta}\boldsymbol{e}, \boldsymbol{Q}_A=EA\left( \int_0^L{\varepsilon _{xx}^{\mathrm{a}}\boldsymbol{r}_{x}^{\mathrm{T}}\boldsymbol{S}_x\mathrm{d}x} \right) 
\\
\boldsymbol{Q}=E\left( \int_0^L{\left( I\kappa \frac{\partial \kappa}{\partial \boldsymbol{e}}+A\varepsilon _{xx}^{\mathrm{a}}\boldsymbol{r}_{x}^{\mathrm{T}}\boldsymbol{S}_x \right) \mathrm{d}x} \right) 
\\
8.96{{\mathrm{g}}\Bigg/{\mathrm{cm}^3}}=8.96\times 10^{-3}{{\mathrm{kg}}\Bigg/{\left( 10^{-2}\mathrm{m} \right) ^3}}=8.96\times 10^3{{\mathrm{kg}}\Bigg/{\mathrm{m}^3}}
\\
V=gA\int_0^L{\rho r_3\mathrm{d}x}=gA\int_0^L{\rho \left[ \boldsymbol{Se} \right] _3\mathrm{d}x}=gA\int_0^L{\rho \left[ \boldsymbol{e}^{\mathrm{T}}\boldsymbol{S}^{\mathrm{T}} \right] _3\mathrm{d}x}
\\
V=-\int_0^L{A\rho \boldsymbol{r}^{\mathrm{T}}\boldsymbol{g}\mathrm{d}x}=-\int_0^L{A\rho \boldsymbol{r}^{\mathrm{T}}\mathrm{d}x\boldsymbol{g}}=-\int_0^L{A\rho \left( \boldsymbol{e}^{\mathrm{T}}\boldsymbol{S}^{\mathrm{T}} \right) \mathrm{d}x\boldsymbol{g}}
\\
\boldsymbol{G}=-\frac{\partial V}{\partial \boldsymbol{e}^{\mathrm{T}}}=\left( \int_0^L{A\rho \boldsymbol{S}^{\mathrm{T}}\mathrm{d}x} \right) \boldsymbol{g}
\\
\boldsymbol{S}_g=\left( \int_0^L{A\rho \boldsymbol{S}\mathrm{d}x} \right) /m
\\
\boldsymbol{G}=\boldsymbol{S}_{g}^{\mathrm{T}}m\boldsymbol{g}
\\
L=T-V
\\
\frac{\mathrm{d}}{\mathrm{d}t}\left( \frac{\partial T}{\partial \dot{\boldsymbol{q}}^{\mathrm{T}}} \right) +\frac{\partial V}{\partial \boldsymbol{q}}=0
\\
\frac{\mathrm{d}}{\mathrm{d}t}\left( \frac{\partial T}{\partial \dot{\boldsymbol{q}}^{\mathrm{T}}} \right) -\boldsymbol{Q}=0
\\
\boldsymbol{Q}=-\frac{\partial V}{\partial \boldsymbol{q}}
\\
\boldsymbol{r}_x=\boldsymbol{S}_x\boldsymbol{e}
\\
\left( -\frac{3}{4}+\frac{3}{4}\xi ^2 \right) \boldsymbol{r}_i
\\
\left( +\frac{3}{4}-\frac{3}{4}\xi ^2 \right) \boldsymbol{r}_j
\\
\frac{L}{8}\left( -1-2\xi +3\xi ^2 \right) \boldsymbol{r}_{x,i}
\\
\frac{L}{8}\left( -1+2\xi +3\xi ^2 \right) \boldsymbol{r}_{x,j}
\\

\\

\\

$$
