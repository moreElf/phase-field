import numpy as np
import matplotlib.pyplot as plt
import os

def load_phi_file(filename):
    """
    Load a phi_XXXX.d file into a 2D numpy array.
    """
    data = np.loadtxt(filename)
    Lx = int(np.max(data[:, 0]) + 1)
    Ly = int(np.max(data[:, 1]) + 1)
    phi = np.zeros((Lx, Ly))
    for x, y, val in data:
        phi[int(x), int(y)] = val
    return phi

def plot_phi(phi, title="φ distribution", save_filename=None):
    """
    Display the phi distribution as a heatmap and save as image.
    """
    plt.figure(figsize=(6, 6))
    plt.imshow(phi, cmap='coolwarm', origin='lower')
    plt.colorbar(label="φ")
    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    
    if save_filename:
        plt.savefig(save_filename, dpi=150, bbox_inches='tight')
        print(f"Plot saved as {save_filename}")
    else:
        plt.show()

def plot_multiple_phi(filenames, save_filename="phi_evolution.png"):
    """
    Plot multiple phi distributions in a grid layout.
    """
    n_files = len(filenames)
    cols = min(4, n_files)
    rows = (n_files + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 4*rows))
    if n_files == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes.reshape(1, -1)
    
    for i, filename in enumerate(filenames):
        if os.path.exists(filename):
            phi = load_phi_file(filename)
            row, col = i // cols, i % cols
            ax = axes[row, col] if rows > 1 else axes[col]
            
            im = ax.imshow(phi, cmap='coolwarm', origin='lower')
            ax.set_title(filename)
            ax.axis("off")
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        else:
            print(f"File '{filename}' not found.")
    
    # Hide unused subplots
    for i in range(n_files, rows * cols):
        row, col = i // cols, i % cols
        ax = axes[row, col] if rows > 1 else axes[col]
        ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(save_filename, dpi=150, bbox_inches='tight')
    print(f"Evolution plot saved as {save_filename}")
    plt.close()

if __name__ == "__main__":
    # Single file visualization
    filename = 'phi_0000.d'
    if os.path.exists(filename):
        phi = load_phi_file(filename)
        plot_phi(phi, title=f"Visualization of {filename}", save_filename=f"{filename.replace('.d', '.png')}")
    else:
        print(f"File '{filename}' not found. Please place it in the current directory.")
    
    # Multiple file visualization (evolution over time)
    # Select files at specific time intervals
    evolution_files = ['phi_0000.d', 'phi_0010.d', 'phi_0020.d', 'phi_0030.d', 
                      'phi_0040.d', 'phi_0050.d', 'phi_0060.d', 'phi_0070.d']
    existing_files = [f for f in evolution_files if os.path.exists(f)]
    
    if existing_files:
        plot_multiple_phi(existing_files, "phi_evolution.png")
    else:
        print("No evolution files found for multi-plot visualization.")
