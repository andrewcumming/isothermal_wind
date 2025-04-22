import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.animation import FuncAnimation
from scipy import optimize

sys.path.insert(0, '../athena/vis/python')
import athena_read

adiabatic = True

# Function to read data
def read_data(i):
    data = athena_read.athdf(f"out/isowind.out1.{i:05d}.athdf")
    x = data['x1v']
    rho = data['rho'][0, 0, :]
    vel = data['vel1'][0, 0, :]
    return x, rho, vel
        
# Initialize the plot
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 8), gridspec_kw={'hspace': 0})
ax1.set_yscale('log')
#ax1.set_xscale('log')
#ax2.set_xscale('log')

ax1.set_ylabel(r'$\rho$')
ax2.set_ylabel(r'$v$')
ax2.set_xlabel('r')

ax2.set_ylim((-3,10))

ax2_right = ax2.twinx()
ax2_right.set_ylabel(r'$\dot{M}$')

line1, = ax1.plot([], [], label='Density')
line2, = ax2.plot([], [], label='Velocity')
line3, = ax2_right.plot([], [], 'r', label=r'$\dot{M}$')
x, rho, vel = read_data(0)
print("stellar surface=", x[0])
Mdot_guess = 4*np.pi * np.exp(1.5) * np.exp(-2/x[0])
ax2_right.plot([x[0],x[-1]], [Mdot_guess, Mdot_guess], 'r:')

def analytic(u,r):
    return u*u - np.log(u*u) - 4*np.log(r) -4/r +3
u_guess = np.array([])
for r in x:
    if len(u_guess) == 0:
        x0 = 0.01
    elif r>1:
        x0 = 2
    else:
        x0 = u_guess[-1]
    u = optimize.fsolve(analytic, x0=x0, args=(r,))     
    u_guess = np.append(u_guess, u)
ax2.plot(x,u_guess,':')

rho_guess = Mdot_guess / (4*np.pi*x**2*u_guess)
ax1.plot(x,rho_guess,':')

# Set up axes ticks and limits (optional)
for ax in [ax1, ax2, ax2_right]:
    ax.tick_params(which='both', direction='in', right=True, top=True)
ax2.tick_params(which='both', direction='in', left=True, right=False, labelleft=True)

# Function to update the plot
def update(i):
    x, rho, vel = read_data(i)
    Mdot = 4*np.pi*x**2 * rho * vel

    # Update data for each line
    line1.set_data(x, rho)
    line2.set_data(x, vel)
    line3.set_data(x, Mdot)

    if 0:
        # Rescale axes
        ax1.relim()
        ax1.autoscale_view()
        ax2.relim()
        ax2.autoscale_view()
        ax2_right.relim()
        ax2_right.autoscale_view()

    titl.set_text(f"Iteration: {i}")
    
    plt.savefig(f'png/frame_{i:05d}.png')
    
    return line1, line2, line3, titl

titl = ax1.text(0.5, 1.05, "Iteration: 0", transform=ax1.transAxes, 
                      ha="center", fontsize=14)

# Create animation
ani = FuncAnimation(fig, update, frames=range(201), interval=10, blit=False, repeat=False)

# Show the animation
plt.tight_layout()
plt.show()
