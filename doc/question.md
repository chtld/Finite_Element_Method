### 问题描述

两点边值问题

$ -\frac{d}{dx}(c(x)\frac{du(x)}{dx}) = f(x), a<x<b$

$u(a) = g_a, u(b) = g_b$

对此方程两边同时乘以测试函数$v(x)$, 并在区间$[a,b]$上积分  则有

$-\int^b_a(cu')'vdx = \int^b_afvdx$

利用分部积分公式有

$\int^b_a(cu')'vdx = \int^b_avd(cu' ) = cu'v|^b_a-\int^b_acu'v'dx$ 

于是得到下式

$-c(b)u'(b)v(b)+c(a)u'(a)v(a)+\int^b_acu'v'dx=\int^b_afvdx$

令$a(u,v)=\int^b_acu'v', F(v)=\int^b_afvdx+c(b)u'(b)v(b)-c(a)u'(a)v(a)$

对Dirichlet边界条件

