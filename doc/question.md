### 问题描述

两点边值问题

$ -\frac{d}{dx}(c(x)\frac{du(x)}{dx}) = f(x), a<x<b$

为保证解的适定性(存在、唯一、稳定)还需要加上适当的边界条件, 如

- Dirichlet边界条件

  $u(a)=g_a, u(b)=g_b.$

- Neumann边界条件

  $u(a)=g_a, u'(b)=r_b$或$u'(a)=r_a, u(b)=g_b.$

- Robin边界条件

  $u(a)=g_a, u'(b)+q_bu(b)=p_b$或$u'(a)+q_au(a)=p_a, u(b)=g_b.$

### 弱形式

对此方程两边同时乘以测试函数$v(x)$, 并在区间$[a,b]$上积分  则有

$-\int^b_a(cu')'vdx = \int^b_afvdx$

利用分部积分公式有

$\int^b_a(cu')'vdx = \int^b_avd(cu' ) = cu'v|^b_a-\int^b_acu'v'dx$ 

于是得到下式

$-c(b)u'(b)v(b)+c(a)u'(a)v(a)+\int^b_acu'v'dx=\int^b_afvdx$

令$a(u,v)=\int^b_acu'v'dx, F(v)=\int^b_afvdx+c(b)u'(b)v(b)-c(a)u'(a)v(a)$

测试函数的选择需要满足最终形成的方程的个数等于未知数的个数, 因此我们如果想要求解该点处的未知量就需要对该点做测试, 如果该点为已知点，则不做测试，即令该点的测试函数的值为0.

- Dirichlet边界条件:$u(a) = g_a, u(b) = g_b$.

   选择测试函数$v(x)$使得$v(a)=v(b)=0$, 则此时

  $a(u,v) = \int^b_acu'v'dx, F(v) = \int^b_afvdx.$


- Neumann边界条件:$u(a) = g_a, u'(b) = r_b$.

  选择测试函数$v(x)$使得$v(a) = 0$, 则此时

  $a(u,v) = \int^b_acu'v'dx, F(v) = \int^b_afvdx+c(b)u'(b)v(b)=\int^b_afvdx+c(b)r(b)v(b)$


- Robin边界条件:$u(a) = g_a, u'(b)+q_bu(b)=p_b$.

  选择测试函数$v(x)$使得$v(a)=0$,则此时

  $a(u,v) = \int^b_acu'v'dx, F(v) = \int^b_afvdx+c(b)u'(b)v(b)=\int^b_afvdx+c(b)(p_b-q_bu(b))v(b)$.

我们称$a(u,v)=F(v)$为上述方程的弱形式. 那么$u, v$应该属于什么函数空间? 为了保证上面的各式是有意的, 则在上面的推导过程中所有的积分应该是有限的, 因此所给的函数$u, v$应该满足$u\in H^1(I), $$v\in H^1(I), v(a)=v(b)=0$. 可以严格写下方程的弱形式:

&emsp;&emsp;求$u\in H^1(I)$使得$a(u,v)=F(v)$对任意$v\in H^1_0(I)$成立, 其中$I=(a,b)$.

### 有限元离散

&emsp;&emsp;上述函数空间均是无穷维空间, 在实际计算中我们无法求解, 有限元Galerkin方法的思想就是用有限维空间去逼近无穷维空间, 假设存在有限维空间$U_h\subset H^1[a,b]$, 则有如下的Galekin形式:

&emsp;&emsp;求$u_h\in U_h$使得$a(u_h,v_h)=F(v_h)$对任意$v_h\in U_h$.

对区间$[a,b]$作一致剖分, 分为$N$等分, 则分割的大小为$h=\frac{(b-a)}{N}$. $x_i=a+(i-1)h,$ $i=1,...,$$N+1$表示网格节点, $E_n=[x_n, x_{n+1}], n=1,...,N.$表示网格单元.

定义一维线性有限元空间:

&emsp;&emsp;$U_h= \{\phi\in C[a,b], \phi(x)在每个网格单元E_n上是线性的, n=1,...,N.\}$

对有限维空间总可以找到它的一组基将它线性表出, 我们可以证明$U_h$是$C[a,b]$的$N+1$维子空间. 证明如下:

&emsp;&emsp;首先, $U_h$是$C[a,b]$的子空间是显然的(真的是显然的).

&emsp;&emsp;其次, 我们可以找到$U_h$的$N+1$个连续的分段线性基函数, 则证明完成.

&emsp;&emsp;考虑$\phi_j(x_i)=\delta_{ij} = \begin{cases} 1, &\text{if $i=j$} \\ 0, &\text{if $i\neq j$} \end{cases}$, 帽子函数(hat function), $$

事实上, 

&emsp;&emsp;$\phi_1(x)=\begin{cases}\frac{x_2-x}{h}, &\text{if  $x_1\leq x\leq x_2$}, \\ 0, &\text{otherwise}.\end{cases}$

&emsp;&emsp;$\phi_j(x)=\begin{cases}\frac{x-x_{j-1}}{h}, &\text{if $x_{j-1}\leq x\leq x_j$},\\ \frac{x_{j+1}-x}{h},&\text{if $x_j\leq x\leq x_{j+1}$}, \\ 0,&\text{otherwise}.\end{cases}$

&emsp;&emsp;$\phi_{N+1}(x)=\begin{cases} \frac{x-x_N}{h}, &\text{if $x_N\leq x\leq x_{N+1}$},\\ 0,&\text{otherwise}.\end{cases}$

下面说明$\phi_i,i=1,...,N+1.$是$U_h$的一组基(只需要说明两点1.线性无关, 2.任意的$f\in U_h$, 均能由其唯一的线性表出):

&emsp;&emsp;考虑$\sum_{j=1}^{N+1}c_j\phi_j(x)=0,\text{for $\forall x \in [a,b].$}$令$x=x_i,(i=1,...,N+1.)$ , 则有$c_i=0,i=1,...,N+1.$所以其线性无关.

&emsp;&emsp;$\forall f \in U_h$, 考虑$g(x)=\sum_{j=1}^{N+1}f(x_i)\phi_j(x).$则有

$g(x_i)=f(x_i),i=1,...,N+1.$且由于$f(x)$和$g(x)$在$[x_i,x_{i+1}],i=1,...,N.$是线性的, 所以

&emsp;&emsp;$f(x)=g(x), [x_i,x_{i+1}], i=1,...,N.$

即, $f(x)=g(x)=\sum_{i=1}^{N+1}f(x_i)\phi_i(x)$.

证毕!

