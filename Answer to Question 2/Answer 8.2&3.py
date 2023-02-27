import utili
import matplotlib.pyplot as plt

'''
Input: 
(1) The initial concentration of E, S, ES, and P.
(2) The rate of constants: k1, k2, and k3.
(3) The upper limit of calculation: verh.
(4) The length of step: h.

e, s, c and p denote the concentration of E, S, ES, and P.
1,000,000 uM = 1 mol
'''
e = 1
s = 10
c = 0
p = 0

k1 = 100
k2 = 600
k3 = 150

verh = 5
h = 0.001
t = 1


'''
Output:
(1) The concrete value of s, c, e, and p.
(2) Velocity (v) was defined as the quotient of delta S and delta t.
(3) Figure 1: The visualization of s, c, e, and p.
(4) Figure 2: The visualization of the relaitionship between V and c.
(5) Label the maximum of V in the Figure 2.
'''

results = utili.PROCEDURE(verh, h, t, s, c, e, p, k1, k2, k3)
t = results[:, 0]
s = results[:, 1]
c = results[:, 2]
e = results[:, 3]
p = results[:, 4]
v = []

for i in range(p.shape[0]):
    if i != 0:
        v_value = (p[i] - p[i-1]) / (t[i] - t[i-1])
        v.append(v_value)

# Figure 1: The visualization of s, c, e, and p.
plt.figure(dpi=300)
plt.plot(t, s, label='S')
plt.plot(t, c, label='ES')
plt.plot(t, e, label='E')
plt.plot(t, p, label='P')
plt.xlabel('Minutes')
plt.ylabel('Concentration (uM)')
plt.legend()
plt.savefig('Concentration-Minutes.png', dpi=300)

# Figure 2: The visualization of the relaitionship between V and c.
plt.figure(dpi=300)
plt.plot(s[1:], v)
max_y = max(v)
max_x = s[1:][v.index(max(v))]
plt.plot(max_x, max_y, 'o', markersize=2, color='r')
plt.text(max_x, max_y, f'({max_x:.2f}, {max_y:.2f})', ha='center', va='bottom', color='r')
plt.xlabel('the concentration of S')
plt.ylabel('velocity (uM/min)')
plt.grid()
plt.savefig('velocity-S.png', dpi=300)
# plt.show()
