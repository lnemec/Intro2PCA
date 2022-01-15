### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1cdb4e08-99f7-48d9-924b-e105d1b266dc
begin
    using Plots; gr()
	using Printf
	using LaTeXStrings
	using PlutoUI
	using LinearAlgebra
	import Statistics: mean, cov
	using MultivariateStats
	using Random
	
	md"""**Import julia libraries**"""
end

# ╔═╡ 721e5d79-4842-43cf-a2ff-4cf0fdb13f87
begin
    rng = MersenneTwister(220115)
    Base.show(io::IO, f::Float64) = @printf io "%1.2f" f

	md"""Set the random seed and change the format of floating numbers for showing arrays."""
end

# ╔═╡ 969a133a-9861-4248-92f0-757fd12f12fb
md"# Principal Component Analysis (PCA): *A Physically Intuitive Mathematical Introduction*

The principal component analysis (PCA) involves rotating a cloud of data points in Euclidean space such that the variance is maximal along the first axis (the so-called first principal component). The principal axis theorem ensures that the data can be rotated in such away. In mathematical terms, the PCA is a coordinate transformation, or equivalently, a change of basis or orthogonal linear transformation.

The mathematics behind the PCA is found again in the description of rotations of rigid bodies. This physical interpretation is instructive in understanding the PCA."

# ╔═╡ 624fcb2a-8446-42c2-863f-e285f2c1fd86
begin
	N_slider = @bind N html"<input type=range min=10 max=1000 step=10 value = 500>"
	R1_slider = @bind R1 html"<input type=range min=0.1 max=10 step=0.1 value=4.0>"
	R2_slider = @bind R2 html"<input type=range min=0.1 max=10 step=0.1 value=8.5>"
	angle_slider = @bind angle html"<input type=range min=0.0 max=180.0 step=0.01 value=35.0>"
	
	md"""## Define parameters to generate sample data
		
	"""
end

# ╔═╡ d8005539-13d5-4a3a-ad62-861c4f4bb9e0
md"First, we will generate a cloud of $N$ randomly distributed data points in Euclidean space $\mathbb{R}^n$ with coordinates $\{ \vec{x}^{(1)}, \vec{x}^{(2)},\dots, \vec{x}^{(N)} \} = X$. We will demonstrate the concept based on $3$-dimensional data, where $n=3$ and $X \subset \mathbb{R}^3$ with basis vectors $\vec{e}_1$, $\vec{e}_2$ and $\vec{e}_3$ centered around the origin $(0, 0, 0)$. For simplicity, we will set all coordinates along the basis vector $\vec{e}_3$ to zero. It will allow us to visualise the data in the $\vec{e}_1$ and $\vec{e}_2$ plane."

# ╔═╡ 5619b1f9-38d1-42de-bd81-708567e2fa96
md"""
**Define Parameters:**

|       | Min | Slider        | Max   |
|-------|-----|:-------------:|-------|
| N     | 10  | $(N_slider)   | 1000  |
| R1    | 0.1 | $(R1_slider)  | 10    |
| R2    | 0.1 | $(R2_slider)  | 10    |
| angle | 0.0 |$(angle_slider)| 180   |
"""

# ╔═╡ 24f9bdd6-0337-4a23-9373-cfecabc2b029
begin
	Rmax = max(R1, R2)
	
	md"""
	
	N defines the number of points of the sample data. R1 and R2 define the
	variance of the data and angle rotates the sample data in space.

	| Variable | N  | R1  | R2  | angle | Rmax |
	|----------|----|-----|-----|-------|------|
	|**Values**| $(N) | $(R1) | $(R2) | $(round(angle, digits=2)) |$(Rmax) |
	"""
end

# ╔═╡ eacd1340-79af-4985-bf60-b8db8a380a20
begin
	rot = [cos(deg2rad(angle)) sin(deg2rad(angle));
		  -sin(deg2rad(angle)) cos(deg2rad(angle))]
	
	pts = randn(rng, N, 3)
	pts[:, 1] .*= R1
	pts[:, 2] .*= R2
    pts[:, 3] = zeros(N)
	
	pts[:, 1:2] = pts[:, 1:2] * rot
	
	md"""Create $3$-dim array with sample data based on parameters defined in the table above.
	
	---
	"""
	
end

# ╔═╡ 6bcb1adc-633c-4730-b4e2-6b5b6f5237c4
begin
	lim = [-3 * Rmax, 3 * Rmax]
	
	u1 = [0.0 0.0; lim[2] lim[2]*tan(deg2rad(angle+90.0))]
	u1[2,:] = normalize!(u1[2,:])
	ann1 = (u1[2,1].*2.5*Rmax, u1[2,2].*2.5*Rmax, Plots.text(L"$\mathbf{\vec{u}_1}$", :left))
	
	u2 = [0.0 0.0; lim[2] lim[2]*tan(deg2rad(angle))]
	u2[2,:] = normalize!(u2[2,:])
	ann2 = (u2[2,1].*2.5*Rmax, u2[2,2].*2.5*Rmax, Plots.text(L"$\mathbf{\vec{u}_2}$", :left))
	
	l = @layout [  a{1.0w, 0.01h}
	              [b{0.9w, 0.1h} _
	               c{0.9w, 0.9h} d{0.1w, 0.9h}]]

	title = plot(title = "Random Set of Points",
		         framestyle = :none,
		         grid = false,
		         showaxis = false,
		         titlelocation = :center,
		         bottom_margin = -40Plots.px)
	
	h1 =  histogram(pts[:,1],
		            bins=range(lim[1], stop = lim[2], length = 20),
		            normalize = true,
		            fillcolor = "#0044aa",
		            orientation = :vertical,
		            framestyle = :none,
	                title = L"$a)$",
	                titlefontsize = 10,
		            titlelocation=:left)
	
	h2 =  histogram(pts[:,2],
		            bins=range(lim[1], stop = lim[2], length = 20),
		            normalize = true,
		            fillcolor = "#0044aa",
		            orientation = :horizontal,
		            framestyle = :none,
	                title = L"$c)$",
	                titlefontsize = 10,
		            titlelocation=:right )
	
	s = scatter(pts[:,1], pts[:,2],
		     markershape = :circle,
		     markersize = 2,
		     markercolor = "#0044aa",
		     markerstrokecolor = "#0044aa",
		     framestyle = :orign,
		     xlabel = L"\vec{e}_1",
	         ylabel = L"\vec{e}_2",
		     xlims = (lim[1], lim[2]),
		     ylims = (lim[1], lim[2]),
		     aspect_ratio = :equal,
		     title = L"$b)$",
	         titlefontsize = 10,
		     titlelocation=:left)
	
	plot(title, h1, s, h2,
		 layout = l,
		 legend = false,
		 size = (600, 610))

	plot!(lim, lim*tan(deg2rad(angle)),
		  color="#a6a6a6",
		  label="",
		  annotations = ann1,
		  subplot = 3 )
	
	plot!(lim, lim*tan(deg2rad(angle+90.0)),
		  color="#a6a6a6",
		  label="",
		  annotations = ann2,
		  subplot = 3)
end

# ╔═╡ 075b42b8-476b-4436-a886-dcbb082b0ad7
md"*Figure 1. a) Distribution of random points along $\vec{e}_2$. b) Cloud of randomly distributed data points in Euclidean space $\mathbb{R}^3$ with coordinates $\{ \vec{x}^{(1)}, \vec{x}^{(2)},\dots, \vec{x}^{(N)} \}$ with the basis vectors (axis in the plot b)) $\vec{e}_1$ and $\vec{e}_2$. The grey lines in plot b) indicate the principal axis $\vec{u}_1$ and $\vec{u}_2$. c) Distribution of random points along $\vec{e}_1$. All values along $\vec{e}_3$ are equal zero.*

---
"

# ╔═╡ 35f8cd6a-a87d-4af9-8590-3bcac79ef3ea
md"## Moment of Inertia $\boldsymbol{J}$

The moment of inertia $\boldsymbol{J}$ of a rigid body, also called rotational inertia, determines the torque required for a desired angular acceleration around an axis of rotation. It depends on the mass distribution of the body and the selected axis of rotation. A body with a larger moment of inertia $\boldsymbol{J}$ requires more torque to change the body's rate of rotation. For the same rigid body, different axes of rotation will have different moments of inertia associated with them. In other words, it depends on the body's mass distribution and the axis chosen, with larger moments requiring more torque to change the body's rotation rate.

All moments of inertia of a rigid body can be summarised by a matrix. In general, it can be determined with respect to any point in space. For simplicity, we will calculate the moment of inertia with respect to the center of mass.

The principal axes of the body, also known as figure axes, and the principal moments of inertia can be found by rotating the cloud of point masses. In mathematical terms, the principal moments are calculated by a coordinate transformation, or equivalently, a change of basis or orthogonal linear transformation.

The figure axis corresponding to the largest principal moment of inertia is the area vector of the plane with the maximal spread of mass points.

| ![Rolling cylinder with different moments of inertia](https://upload.wikimedia.org/wikipedia/commons/9/95/RollingVsInertia.gif) | 
|:--| 
| *Figure 2. The six cylinders have the same mass but different moments of inertia $\boldsymbol{J}$. As they roll down the slope, cylinders with lower moments of inertia accelerate more quickly.  ([Image taken from wikipedia](https://en.wikipedia.org/wiki/File:RollingVsInertia.gif))* |
"

# ╔═╡ 9ee5e06b-318b-4750-82bb-9db2999f2d1e
md"## Visual Comparison of the Principal Axis in PCA and Moment of Inertia

In the following, we will interpret the set of random data points from above in two ways. First, we interpret the data points $X$ as statistically distributed data with Covariance matrix $\boldsymbol{C}$. Second, $X$ represents a set of point masses representing a rigid body with the Moment of Inertia matrix $\boldsymbol{J}$. 

---"

# ╔═╡ 5a400a58-4bfb-49fe-8286-2b8e56803710
begin
	l2 = @layout [c{0.5w, 1.0h} d{0.5w, 1.0h}]

	s2 = scatter(pts[:,1], pts[:,2],
		     markershape = :circle,
		     markersize = 2,
		     markercolor = "#0044aa",
		     markerstrokecolor = "#0044aa",
		     framestyle = :orign,
		     xlabel = L"\vec{e}_1",
	         ylabel = L"\vec{e}_2",
		     xlims = (lim[1], lim[2]),
		     ylims = (lim[1], lim[2]),
		     aspect_ratio = :equal)
	
	plot(s2, s2,
		 layout = l2,
		 legend = false,
		 size = (600, 310))

	plot!(lim, lim*tan(deg2rad(angle+90.0)),
		  color="#a6a6a6",
		  label="",
		  annotations = ann1,
		  subplot = 1,
		  title = "a) PCA",
	      titlefontsize = 12,
		  titlelocation=:left )

	plot!( u1[:,1].*Rmax,
		   u1[:,2].*Rmax,
		   subplot = 1,
		   arrow=true,color="#6a0000",linewidth=3,label="")
	
	
	plot!(lim, lim*tan(deg2rad(angle)),
		  color="#a6a6a6",
		  label="",
		  annotations = ann2,
		  subplot = 2,
		  title = "b) Moment of Inertia",
	      titlefontsize = 12,
		  titlelocation=:left)

	plot!( u2[:,1].*Rmax,
		   u2[:,2].*Rmax,
		   subplot = 2,
		   arrow=true, color="#01493e",linewidth=3,label="")
	
end

# ╔═╡ 07371ed0-72e9-47d9-82d8-6f3cb0e16ec8
md" *Figure 3. a) PCA: The first principal component $\vec{u}_1$ is along the axis where the variance is maximally indicated by a red arrow b) Moment of Inertia: The figure axis $\vec{u}_2$ corresponding to the (actually second) largest principal moment of inertia indicated by a green arrow.*

---
"

# ╔═╡ 2f20596e-efa9-4d00-9769-6e8689c99f05
md"With the visual support of *Figure 1* and *3*, we expect that the principal axes of the PCA and Moment of Inertia are the same. However, the value of the largest principal component and principal moment of inertia will differ for most sets of data points.

> *Note:* In physics, the moment of inertia is defined for a $3$-dimensional rigid body. For simplicity, we projected the data into the plane spanned by $\vec{e}_1$ and $\vec{e}_2$. In our example, the plane with the maximal spread of mass points is then spaned by $\vec{e}_1 \times \vec{e}_2$. The principal axis corresponding to the largest moment is pointing out of the plane along $\vec{e}_3$ and is orthogonal to $\vec{e}_1$ and $\vec{e}_2$.

The line along $\vec{u}_1$ is equivalent to the direction where the variance is maximal. In the following, we will support our visual understanding by exploring the PCA and moment of inertia mathematically.
"

# ╔═╡ 03664f5c-d45c-11ea-21b6-91cd647a07aa
md"## Definition of the Moment of Inertia of a Ridgid Body

For a rigid object of $N$ point masses $m_i$ in $\mathbb{R}^n$, the moment of inertia  $\boldsymbol{J}$ is given by

$$\begin{equation*}
    \boldsymbol{J} = 
        \begin{bmatrix}
            J_{1,1} & \cdots & J_{1,n} \\
            \vdots  & \ddots & \vdots  \\
            J_{n,1} & \cdots & J_{n,n} 
        \end{bmatrix}
\end{equation*}$$

Its components are defined by Eq. (1) as

$$\begin{equation}
    J_{j,j'} =
        \frac{1}{M}
        \sum_{i=1}^{N} m_{i}
        \left( ||\vec{x}^{(i)}||^2 \delta_{j,j'} - x_j^{(i)} x_{j'}^{(i)} \right)
\end{equation}$$

where $\delta_{j,j'}$ is the Kronecker delta and $M = \sum_{i}^{N} m_i$ is the total mass.

> *Note:* Here, we normalise the moment of inertia by the total mass. In physics, the moment of inertia would not commonly be normalised like this.

Looking closer, we see that $\boldsymbol{J}$ is symmetric with $J_{j,j'} = J_{j',j}$. The spectral theorem tells us that $\boldsymbol{J}$ has real eigenvalues $\lambda$ and is diagonalisable by an orthogonal matrix (orthogonally diagonalizable).

"

# ╔═╡ 2e1aa912-e7bd-4ffd-8060-05556a77bc8b
md"## Definition of the Covariance Matrix

The covariance matrix $\boldsymbol{C}$ for a cloud of points in Euclidean space centered around the mean is given by

$$\begin{equation*}
    \boldsymbol{C} = 
        \begin{bmatrix}
            C_{1,1} & \cdots & C_{1,n} \\
            \vdots  & \ddots & \vdots  \\
            C_{n,1} & \cdots & C_{n,n} 
        \end{bmatrix}
\end{equation*}$$

Its components are defined by Eq. (2) as

$$\begin{equation}
    C_{j,j'} =
        \frac{1}{N}
        \sum_{i=1}^{N}
        \left( x_j^{(i)} x_{j'}^{(i)} \right)
\end{equation}$$"

# ╔═╡ 8c8167ca-8679-4d2a-b15b-4c8f32cfb8b6
md"## Solving the Eigenvalue Problem

The principal axes of the PCA and moment of inertia can be determined by rotating the data points in space. More precisely, the principal components and axes are calculated by a coordinate transformation, or equivalently, a change of basis or orthogonal linear transformation.

A real symmetric matrix (like $\boldsymbol{C}$ and $\boldsymbol{J}$ ) has the eigendecomposition into the product of a rotation matrix $\boldsymbol{R}$ and a diagonal matrix $\boldsymbol{\Lambda}$

$$\begin{equation*}
    \boldsymbol{\Lambda} = 
        \begin{bmatrix}
            \lambda_{1} & \cdots & 0 \\
            \vdots  & \ddots & \vdots  \\
            0 & \cdots & \lambda_{n} 
        \end{bmatrix}
\end{equation*}$$

given by

$\boldsymbol{J} = \boldsymbol{R} \boldsymbol{\Lambda } \boldsymbol{R}^T$

The columns of the rotation matrix $\boldsymbol{R}$ define the directions of the principal axes, and the constants $\lambda_{1}, \dots, \lambda_{n}$ are the diagonal elements of the matrix $\boldsymbol{\Lambda}$ and
are called the principal moments.

The structure of the matrices $\boldsymbol{J}$ and $\boldsymbol{C}$ is the same except for the sign of the off-diagonal elements. We will see in the following that the eigenvectors will be the same for $\boldsymbol{C}$ and $\boldsymbol{J}$. In addition, we will see how the eigenvalues $\boldsymbol{\Lambda}$ of $\boldsymbol{C}$ relate to the eigenvalues of $\boldsymbol{J}$. 
"

# ╔═╡ 1d1d7052-977e-4180-af75-f1171f1492b2
md"## Showing the Equivalency of the Eigenvectors of $\boldsymbol{C}$ and $\boldsymbol{J}$

Let's rewrite the moment of inertia matrix $\boldsymbol{J}$ defined by Eq. (1) in terms of the covariance matrix $\boldsymbol{C}$ in Eq. (2).

$$\begin{equation}
    \boldsymbol{J} = tr(\boldsymbol{C})\boldsymbol{I} - \boldsymbol{C} \qquad \text{Eq. (3)}
\end{equation}$$

where $\boldsymbol{I}$ is the identity matrix.

$$\begin{align}
J_{j,j'} & =  tr(\boldsymbol{C})\boldsymbol{I_{j,j}} - \boldsymbol{C_{j,j'}} \\

& = \frac{1}{N} \sum_{i=1}^{N} \sum_{j=1}^{n} \left( x_j^{(i)} x_{j}^{(i)} \right) -      \frac{1}{N} \sum_{i=1}^{N} \left( x_j^{(i)} x_{j'}^{(i)} \right)\\

& = \frac{1}{N} \sum_{i=1}^{N} \left( \sum_{j=1}^{n} x_j^{(i)} x_{j}^{(i)}  - x_j^{(i)} x_{j'}^{(i)} \right)\\

& = \frac{1}{N} \sum_{i=1}^{N} \left( ||\vec{x}^{(i)}||^2 \delta_{j,j'} - x_j^{(i)} x_{j'}^{(i)} \right) \\

& \text{multiply the expression under the sum by } 1 \\

& = \frac{1}{N} \sum_{i=1}^{N} 1 \cdot \left( ||\vec{x}^{(i)}||^2 \delta_{j,j'} - x_j^{(i)} x_{j'}^{(i)} \right) \\

& \text{assume each data point has mass } m_{i} = 1 \\


& = \frac{1}{\sum_{i}^{N} m_i} \sum_{i=1}^{N} m_{i} \left( ||\vec{x}^{(i)}||^2 \delta_{j,j'} - x_j^{(i)} x_{j'}^{(i)} \right) \\

& = \frac{1}{M} \sum_{i=1}^{N} m_{i} \left( ||\vec{x}^{(i)}||^2 \delta_{j,j'} - x_j^{(i)} x_{j'}^{(i)} \right) \\

\end{align}$$
"

# ╔═╡ 48d9a159-5725-4c8d-8a3f-334c0f6e24cf
md"To obtain the eigenvectors and eigenvalues of $\boldsymbol{C}$, we solve the eigenvalue problem by decomposing $\boldsymbol{C}$ into the product of a rotation matrix $\boldsymbol{R}$ and a diagonal matrix $\boldsymbol{\Lambda}_{cov}$

$\boldsymbol{C} = \boldsymbol{R} \boldsymbol{\Lambda}_{cov} \boldsymbol{R}^T,$

where $\boldsymbol{R}$ is composed of the eigenvectors $\vec{v}^{j}$. In the case of a $3$-dimensional space

$\boldsymbol{R} = \left[ \vec{v}^{1} \; \vec{v}^{2} \; \vec{v}^{2} \right]$

and

$$\begin{equation*}
    \boldsymbol{\boldsymbol{\Lambda}_{cov}} = 
        \begin{bmatrix}
            \lambda^{1}_{cov} & 0 & 0 \\
            0  & \lambda^{2}_{cov} & 0 \\
            0 & 0 & \lambda^{3}_{cov}
        \end{bmatrix}
\end{equation*}$$.


The $j^{th}$ eigenvector $\vec{v}^{j}$ and $j^{th}$ eigenvalue $\lambda^{j}_{cov}$ of the covariance matrix $\boldsymbol{C}$ are given by

$$\begin{equation}
    \boldsymbol{C}\lambda^{j}_{cov} = \lambda^{j}_{cov} \vec{v}^{j}.
\end{equation}$$

In the following, we will drop the index $j$ and rewrite the formula from above to $\boldsymbol{C}\lambda_{cov} = \lambda_{cov} \vec{v}$. We multiply equation Eq. (3) by the eigenvector $\vec{v}$ from the right side.

$$\begin{align}
    \boldsymbol{J}\vec{v} &= \left( tr(\boldsymbol{C})\boldsymbol{I} - \boldsymbol{C} \right) \vec{v}  \\
    &= tr(\boldsymbol{C})\boldsymbol{I}\vec{v} - \boldsymbol{C} \vec{v} \\

    &= tr(\boldsymbol{C})\boldsymbol{I}\vec{v} - \lambda_{cov} \vec{v} \\

    &= tr(\boldsymbol{C})\vec{v} - \lambda_{cov} \vec{v} \\

    &= \left( tr(\boldsymbol{C}) - \lambda_{cov} \right) \vec{v} \\

    &= \lambda_{J} \vec{v} \qquad \qquad \qquad \qquad \text{Eq. (4)}\\
\end{align}$$

From Eq. (4), we see that $\boldsymbol{C}$ and $\boldsymbol{J}$ have the same eigenvectors $\vec{v}$.
"

# ╔═╡ 9e8ce452-78bc-4f1c-95fa-37be90c1006c
md"## Calculating the Eigenvalues of the Moment of Inertia Matrix

Based on Eq. (4), we first need to make a small side calculation:

$tr\left( \boldsymbol{C} \right)= tr \left( \boldsymbol{R} \boldsymbol{\Lambda}_{cov} \boldsymbol{R}^T \right)$

The trace of a matrix is invariant under cyclic permutations

$$
\begin{align}
    tr \left( \boldsymbol{R} \boldsymbol{\Lambda}_{cov} \boldsymbol{R}^T \right) &= 
    tr \left( \boldsymbol{R}^T \boldsymbol{R} \boldsymbol{\Lambda}_{cov} \right)\\
    &= tr \left( \boldsymbol{I} \boldsymbol{\Lambda}_{cov} \right) \\
    &= \sum_{j=1}^{n} \lambda^{j}_{cov} \qquad \quad \text{Eq. (5)}
\end{align}$$

"

# ╔═╡ d3eb6336-9e2e-4503-8fa1-f49385b45b40
md"From Eq. (4) and Eq. (5), we can calculate the eigenvalues $\lambda_{J}$ of the moment of inertia $\boldsymbol{J}$ (Eq. (1)).

$$\begin{align}
    \lambda^{k}_{J} &= tr(\boldsymbol{C}) - \lambda^{k}_{cov} \\
    &= \sum_{j=1}^{n} \lambda^{j}_{cov} - \lambda^{k}_{cov} \qquad \text{Eq. (6)}
  \end{align}$$

We see that the $k^{th}$ eigenvalue $\lambda^{k}_{J}$ can be calculated from the eigenvalues $\Lambda_{cov}$.
"

# ╔═╡ 96b8d373-3ce1-460a-a70c-3f9e026bb783
md"## Calculating eigenvalues and eigenvectors using the dataset $X$"

# ╔═╡ 42c9d5b9-7939-4f1f-b148-c9347acb3b3d
md"In the following, we calculate the covariance matrix and its eigenvalues and eigenvectors. The covariance matrix $\boldsymbol{C}$ of the dataset $X$ is"

# ╔═╡ cbdf4dd7-5710-47d9-91aa-8dbae9712e26
C = cov(pts; corrected=false)

# ╔═╡ be55c1dd-0b85-4192-bbbb-2f397fae1d49
md"For the cavariance matrix $\boldsymbol{C}$, we find the eigenvalues $\lambda_{cov}$ and the eigenvectors $v_{cov}$"

# ╔═╡ d283412a-a3a3-458c-aa8b-63e1c4f3bf15
λ_cov, v_cov = eigen(C; sortby=-)

# ╔═╡ 4b17638d-1465-4b21-a2a3-61f94cbd0bca
md"The eigenvalues $\lambda_{cov}$ are"

# ╔═╡ 19378277-61f3-4d59-b840-452cfdf58a3b
λ_cov

# ╔═╡ 30513f7e-d866-4078-95e2-ab05bfbcf547
md"and the eigenvectors $\vec{v}_{cov}$ are"

# ╔═╡ 32191a9f-14ed-4bf6-a8ea-6586b117e848
v_cov[:,1], v_cov[:,2], v_cov[:,3]

# ╔═╡ 11c37b73-ddfd-41cf-98f4-c461b265d1bf
md"Using Eq. (6), we can calculate the eigenvalues $\lambda_{J}$ of $\boldsymbol{J}$"

# ╔═╡ b20fed20-a678-450b-a689-a51baa337153
sum(λ_cov)-λ_cov[1], sum(λ_cov)-λ_cov[2], sum(λ_cov)-λ_cov[3]

# ╔═╡ abbd82f5-31cc-4f45-9c7c-f401796174cb
md"## Calculating the Eigenvalues $\lambda_{J}$ and Eigenvectors $\vec{v}$ of $\boldsymbol{J}$"

# ╔═╡ eaa487b9-2728-49f8-8b5e-003926f4d544
δ(i, j) = i == j ? 1 : 0

# ╔═╡ 3e6282ce-f3ba-4591-976d-0c37d14755c7
"""
	calcJ(pts::AbstractMatrix)

Function to calculate the moment of inertia matrix ``\\boldsymbol{J}``
"""
function calcJ(pts::AbstractMatrix)
	N, n = size(pts)
	J = zeros(Float64, n, n)
	
	for j1 in 1:n
		for j2 in j1:n
			J[j1, j2] = sum(1:N) do i
				norm(pts[i, :])^2 * δ(j1, j2) - pts[i, j1] * pts[i, j2]
			end

			J[j2, j1] = J[j1, j2]
		end
	end

	return J ./ N
end

# ╔═╡ 6b58bfb2-be89-4cd2-9c67-85ca95d02267
J = calcJ(pts)

# ╔═╡ 51f20ae2-e0b6-410d-8783-689f6de9a025
md"For the moment of inertia matrix $\boldsymbol{J}$, we find the Eigenvalues $\lambda$ and the eigenvectors $v$"

# ╔═╡ 4b47a6af-3c41-4575-a576-09145c535e22
λ, v = eigen(J; sortby=-)

# ╔═╡ 1ad77ac0-bd4c-407f-858f-6e435085d4d8
md"In both interpretations of the cloud of data points, first, as a point cloud-centered around the mean and second as a rigid body rotating around the center of mass, we obtain the same eigenvalues and eigenvectors."

# ╔═╡ 9a7fdec6-b17b-47f1-9c81-c3d8525f227c
scale_factor =  λ_cov[1]/λ_cov[2]

# ╔═╡ 2f593024-6533-489e-928e-ed91ea382396
md"---
**Plotting the Unit Vectors Scaled by the Scale Factor**
"

# ╔═╡ e54b7b64-ebbe-426d-92fa-4959eebba826
begin
    p1 = plot(title = "Random selection of points",
	         xlims = (-3 * Rmax, 3 * Rmax),
	         ylims = (-3 * Rmax, 3 * Rmax) )
	
    scatter!(pts[:,1], pts[:,2],
		     markershape = :circle,
		     markersize = 2,
		     markercolor = "#0044aa",
		     markerstrokecolor = "#0044aa",
		     framestyle = :zerolines,
		     xlabel = L"\vec{e}_1",
	         ylabel = L"\vec{e}_2",
		     aspect_ratio = :equal,
		     label="Points")

	plot!( [0, v_cov[1,1] * sqrt(λ_cov[1])*2],
		   [0, v_cov[2,1] * sqrt(λ_cov[1])*2],
		    arrow=true,color="#6a0000",linewidth=3,label="")
	
	plot!( [0, v_cov[1,2] * sqrt(λ_cov[2])*2],
		   [0, v_cov[2,2] * sqrt(λ_cov[2])*2],
		    arrow=true, color="#01493e",linewidth=3,label="")
end

# ╔═╡ b6c98fa0-8b40-4654-aa0c-b0204ed9e291
md"*Figure 4. Cloud of randomly distributed data points $X$. Overlayed are the scaled eigenvectors (shown as red and green arrows).*

---
"

# ╔═╡ 3c8c50dc-7b25-48d3-84f2-c2e1a9fb85dd
md"As a next step we use the eigenvectors and eigenvalues to rotate the data and represent it in the new basis $v = v_{cov} = [\vec{u}_1, \vec{u}_2, \vec{u}_3]$."

# ╔═╡ fd8ce478-87a8-48aa-b8d1-df80df93254b
pts_rot = pts * v_cov

# ╔═╡ 8d9840d1-1e71-403b-8935-c6aeca8a0ee0
md"
---
"

# ╔═╡ 2be79841-4bfd-47d1-9c50-fd1e33a766c5
begin
    p2 = plot(title = "Random selection of points",
	         xlims = (-3 * Rmax, 3 * Rmax),
	         ylims = (-3 * Rmax, 3 * Rmax) )
	
    scatter!(pts_rot[:,1], pts_rot[:,2],
		     markershape = :circle,
		     markersize = 2,
		     markercolor = :steelblue,
		     markerstrokecolor = :steelblue,
		     xlabel = L"\vec{u}_1",
	         ylabel = L"\vec{u}_2",
		     aspect_ratio = :equal,
		     label="Points")

	plot!( [0, v_cov[1,1]*10.0],
		   [0, v_cov[2,1]*10.0],
		    arrow=true,color=:grey,linewidth=3,label="")
	
	plot!( [0, v_cov[1,2]*10.0],
		   [0, v_cov[2,2]*10.0],
		    arrow=true,color=:grey,linewidth=3,label="")
end

# ╔═╡ 9a90438b-67c2-4024-bf45-696d60b718df
md"*Figure 4. Cloud of randomly distributed data points $X$ represented in the basis spanned by the vectors $[\vec{u}_1 \; \vec{u}_2 \; \vec{u}_3]$. The grey arrows indicate the basis vectors $\vec{e}_1$ and $\vec{e}_2$*

---
"

# ╔═╡ 14bf83f9-9a19-42e0-ba41-c6a16c58a89c
md"""## Final Thoughts
In machine learning and data science, PCA is used for two reasons.
First, the accuracy and numeric stability of some Machine Learning algorithms are sensitive towards correlated input data. In particular, Machine Learning algorithms that perform an inversion of the covariance matrix may experience the singularity problem (Gaussian Mixture Models come to mind). A different example is the application of random forest algorithms to detect interactions between different features because highly correlated data can mask these interactions. Here, first performing a PCA can support and safeguard the interpretability.


Second, it is used to reduce the dimensionality of the data set for example for data compression. In our example, we used a 3-dimension dataset $X$, however, the $\vec{e}_3$ component by construction does not carry any information. We, therefore, could use the PCA, to drop the third dimension and use the projection of the data on the $\vec{u}_1 \times \vec{u}_2$ only. This becomes a powerful tool to handle high-dimensional data and deal with the [curse of dimensionality](https://en.wikipedia.org/wiki/Curse_of_dimensionality).

Both applications of PCA have been extensively discussed. Below you find recommendations for further reading.

"""

# ╔═╡ 3b409a8a-b03a-41e6-9e14-b9f73af1847b
md"## Further Reading
  1. **Bauckhage, Christian & Dong, Tiansi** A Tutorial on Principal Component Analysis Part 1: Motivation (2018).
  2. **Bauckhage, Christian & Dong, Tiansi** A Tutorial on Principal Component Analysis Part 2: Principal Axis via Moment of Inertia (2018).
  3. **Hong, Liang** [The Proof of Equivalence between \" Moment of Inertia Based Method\" and \"PCA Based Method\" in LOS Detection](10.4028/www.scientific.net/AMR.945-949.2071) Advanced Materials Research. 945-949. 2071-2074 (2014).
  4. **Murphy, Kevin P.** [Machine Learning A Probabilistic Perspective](https://isbnsearch.org/isbn/9780262018029), The MIT Press, Chapter 12.2 (2012)
  5. **Marsland, Stephen** [Machine Learning An Algorithmic Perspective](https://isbnsearch.org/isbn/9781466583283), 2nd Edition, Chapman and Hall/CRC (2014)
  6. [Badr, Will](https://medium.com/@will.badr) [Why Feature Correlation Matters …. A Lot!](https://towardsdatascience.com/why-feature-correlation-matters-a-lot-847e8ba439c4) Blog TowardsDataScience (2019)."

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MultivariateStats = "6f286f6a-111f-5878-ab1e-185364afe411"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
LaTeXStrings = "~1.3.0"
MultivariateStats = "~0.8.0"
Plots = "~1.23.6"
PlutoUI = "~0.7.19"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0bc60e3006ad95b4bb7497698dd7c6d649b9bc06"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra"]
git-tree-sha1 = "2ff92b71ba1747c5fdd541f8fc87736d82f40ec9"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.4.0"

[[Arpack_jll]]
deps = ["Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "e214a9b9bd1b4e1b4f15b22c0994862b66af7ff7"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.0+3"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

[[ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "30f2b340c2fff8410d89bfcdc9c0a6dd661ac5f7"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.62.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "8d958ff1854b166003238fe191ec34b9d592860a"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.8.0"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun"]
git-tree-sha1 = "0d185e8c33401084cab546a756b387b15f76720c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.23.6"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "e071adf21e165ea0d904b595544a8e514c8bb42c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.19"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "eb35dcc66558b2dda84079b9a1be17557d32091a"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.12"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─1cdb4e08-99f7-48d9-924b-e105d1b266dc
# ╟─721e5d79-4842-43cf-a2ff-4cf0fdb13f87
# ╟─969a133a-9861-4248-92f0-757fd12f12fb
# ╟─624fcb2a-8446-42c2-863f-e285f2c1fd86
# ╟─d8005539-13d5-4a3a-ad62-861c4f4bb9e0
# ╟─5619b1f9-38d1-42de-bd81-708567e2fa96
# ╟─24f9bdd6-0337-4a23-9373-cfecabc2b029
# ╟─eacd1340-79af-4985-bf60-b8db8a380a20
# ╟─6bcb1adc-633c-4730-b4e2-6b5b6f5237c4
# ╟─075b42b8-476b-4436-a886-dcbb082b0ad7
# ╟─35f8cd6a-a87d-4af9-8590-3bcac79ef3ea
# ╟─9ee5e06b-318b-4750-82bb-9db2999f2d1e
# ╟─5a400a58-4bfb-49fe-8286-2b8e56803710
# ╟─07371ed0-72e9-47d9-82d8-6f3cb0e16ec8
# ╟─2f20596e-efa9-4d00-9769-6e8689c99f05
# ╟─03664f5c-d45c-11ea-21b6-91cd647a07aa
# ╟─2e1aa912-e7bd-4ffd-8060-05556a77bc8b
# ╟─8c8167ca-8679-4d2a-b15b-4c8f32cfb8b6
# ╟─1d1d7052-977e-4180-af75-f1171f1492b2
# ╟─48d9a159-5725-4c8d-8a3f-334c0f6e24cf
# ╟─9e8ce452-78bc-4f1c-95fa-37be90c1006c
# ╟─d3eb6336-9e2e-4503-8fa1-f49385b45b40
# ╟─96b8d373-3ce1-460a-a70c-3f9e026bb783
# ╟─42c9d5b9-7939-4f1f-b148-c9347acb3b3d
# ╠═cbdf4dd7-5710-47d9-91aa-8dbae9712e26
# ╟─be55c1dd-0b85-4192-bbbb-2f397fae1d49
# ╠═d283412a-a3a3-458c-aa8b-63e1c4f3bf15
# ╟─4b17638d-1465-4b21-a2a3-61f94cbd0bca
# ╠═19378277-61f3-4d59-b840-452cfdf58a3b
# ╟─30513f7e-d866-4078-95e2-ab05bfbcf547
# ╠═32191a9f-14ed-4bf6-a8ea-6586b117e848
# ╟─11c37b73-ddfd-41cf-98f4-c461b265d1bf
# ╠═b20fed20-a678-450b-a689-a51baa337153
# ╟─abbd82f5-31cc-4f45-9c7c-f401796174cb
# ╟─eaa487b9-2728-49f8-8b5e-003926f4d544
# ╟─3e6282ce-f3ba-4591-976d-0c37d14755c7
# ╠═6b58bfb2-be89-4cd2-9c67-85ca95d02267
# ╟─51f20ae2-e0b6-410d-8783-689f6de9a025
# ╠═4b47a6af-3c41-4575-a576-09145c535e22
# ╟─1ad77ac0-bd4c-407f-858f-6e435085d4d8
# ╠═9a7fdec6-b17b-47f1-9c81-c3d8525f227c
# ╟─2f593024-6533-489e-928e-ed91ea382396
# ╠═e54b7b64-ebbe-426d-92fa-4959eebba826
# ╟─b6c98fa0-8b40-4654-aa0c-b0204ed9e291
# ╟─3c8c50dc-7b25-48d3-84f2-c2e1a9fb85dd
# ╠═fd8ce478-87a8-48aa-b8d1-df80df93254b
# ╟─8d9840d1-1e71-403b-8935-c6aeca8a0ee0
# ╟─2be79841-4bfd-47d1-9c50-fd1e33a766c5
# ╟─9a90438b-67c2-4024-bf45-696d60b718df
# ╟─14bf83f9-9a19-42e0-ba41-c6a16c58a89c
# ╟─3b409a8a-b03a-41e6-9e14-b9f73af1847b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
