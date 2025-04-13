# Enzyme Kinetics Solution

## Steady State Concentrations

E = 
$$\frac{Et k2 k3 + Et k3 k_{-1} + Et k_{-1} k_{-2} + 2 Et P k2 k_{-3} + 2 Et P k_{-1} k_{-3} + 2 Et P k_{-2} k_{-3} + 2 Et S k1 k2 + 2 Et S k1 k3 + 2 Et S k1 k_{-2}}{k2 k3 + k3 k_{-1} + k_{-1} k_{-2} + P k2 k_{-3} + P k_{-1} k_{-3} + P k_{-2} k_{-3} + S k1 k2 + S k1 k3 + S k1 k_{-2}}$$

ES = 
$$\frac{- Et P k_{-2} k_{-3} - Et S k1 k3 - Et S k1 k_{-2}}{k2 k3 + k3 k_{-1} + k_{-1} k_{-2} + P k2 k_{-3} + P k_{-1} k_{-3} + P k_{-2} k_{-3} + S k1 k2 + S k1 k3 + S k1 k_{-2}}$$

EP = 
$$\frac{- Et P k2 k_{-3} - Et P k_{-1} k_{-3} - Et S k1 k2}{k2 k3 + k3 k_{-1} + k_{-1} k_{-2} + P k2 k_{-3} + P k_{-1} k_{-3} + P k_{-2} k_{-3} + S k1 k2 + S k1 k3 + S k1 k_{-2}}$$

## Reaction Rate

### Original Form (from PDF)

Net reaction rate = 
$$\frac{-Et(P k_{-1} k_{-2} k_{-3} - S k1 k2 k3)}{k2 k3 + k3 k_{-1} + k_{-1} k_{-2} + P k2 k_{-3} + P k_{-1} k_{-3} + P k_{-2} k_{-3} + S k1 k2 + S k1 k3 + S k1 k_{-2}}$$

### Simplified Form with Kinetic Parameters

Net reaction rate = 
$$\frac{Et k_{forward} S - Et k_{reverse} P}{1 + \frac{S}{K_s} + \frac{P}{K_p}}$$

where:

## Kinetic Parameters

### Equilibrium Constants

Kp = $$k2 k3 + k3 k_{-1} + k_{-1} k_{-2}$$

Ks = $$\frac{k2 k3 + k3 k_{-1} + k_{-1} k_{-2}}{k1(k2 + k3 + k_{-2})}$$

### Rate Constants

k_forward = $$\frac{k2 k3}{k2 + k3 + k_{-2}}$$

k_reverse = $$\frac{k_{-1} k_{-2}}{k2 + k_{-1} + k_{-2}}$$

## Special Cases

### Initial Velocity (P = 0)

v_0 = $$\frac{Et \cdot k_{forward} \cdot S}{1 + \frac{S}{K_s}}$$

### Product Inhibition (S = 0)

v_{reverse} = $$\frac{-Et \cdot k_{reverse} \cdot P}{1 + \frac{P}{K_p}}$$

