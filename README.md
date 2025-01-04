# The Schwarz-Christoffel transformation

A conformal map $f: \mathbb{C} \rightarrow \mathbb{C}$ is an angle-preserving transformation of points located in the complex plane.
Conformal maps are particularly useful for solving electrostatics problems governed by Laplace's equation since solutions $w$ are harmonic and conformal maps preserve this property.
Hence, if an electrostatics boundary value problem can be solved in one geometry, it can be mapped onto a different boundary geometry.
The Schwarz-Christoffel transformation implements the mapping from a circular boundary to an arbitrary polygon.

This library implements

$$w = f(z) = c \int_0^z \prod_{k=1}^N \left(1-\frac{\tilde{z}}{z_k}\right)^{-\beta_k}\,{\rm d}\tilde{z} \ , $$

where $c$ is a constant for scaling and rotating, and $z_k \in \mathbb{C}$ are so-called prevertices which map to the vertices $w_k = f(z_k)$ of the polygon.
Finding the prevertices to a given polygon is called the Schwarz-Christoffel parameter problem and is also implemented.

<img width="162" alt="{8E118458-8086-4166-A454-618B9A2918E4}" src="https://github.com/user-attachments/assets/9cd893cc-557f-428f-81ae-94cb67fcd33c" />

Vertices at infinity are also supported, enabling open boundary geometries.

<img width="164" alt="{513687E8-3633-4ED5-8128-3E80EC9A6D7B}" src="https://github.com/user-attachments/assets/bc010589-436d-4c0c-b37e-092f92e954d3" />
