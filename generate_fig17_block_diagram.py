#!/usr/bin/env python3
"""
Generate Fig 17: System Block Diagram for Six-State ODE Model
- Larger fonts, no overlap
- Pathway label colors match arrow colors
- Decay arrow outside Stress box
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Create figure - extra large for better readability
fig, ax = plt.subplots(figsize=(40, 28))
ax.set_xlim(0, 110)
ax.set_ylim(0, 80)
ax.set_aspect('equal')
ax.axis('off')

# Colors
colors = {
    'growth': '#90EE90',      # Light green - LAI, C_buf, X_d
    'stress': '#FF6B6B',      # Light red - ROS, Stress
    'anth': '#DDA0DD',        # Plum - Anthocyanin
    'input': '#FFD700',       # Gold - UV-A input
    'damage': '#FFA500',      # Orange - damage mechanisms
    'arrow_positive': '#228B22',  # Forest green
    'arrow_negative': '#DC143C',  # Crimson
    'arrow_neutral': '#4169E1',   # Royal blue
}

# Font sizes - SIGNIFICANTLY LARGER
title_fontsize = 52
box_fontsize = 40          # State box labels - very large
box_detail_fontsize = 28   # Secondary text in boxes
pathway_fontsize = 26      # Pathway labels
legend_fontsize = 26

def draw_state_box(ax, x, y, width, height, label, sublabel, color, fontsize=box_fontsize):
    """Draw a state variable box"""
    box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                         boxstyle="round,pad=0.03,rounding_size=0.8",
                         facecolor=color, edgecolor='black', linewidth=3)
    ax.add_patch(box)

    if sublabel:
        ax.text(x, y + 1.8, label, ha='center', va='center', fontsize=fontsize, fontweight='bold')
        ax.text(x, y - 2.2, sublabel, ha='center', va='center', fontsize=box_detail_fontsize)
    else:
        ax.text(x, y, label, ha='center', va='center', fontsize=fontsize, fontweight='bold')


def draw_arrow(ax, start, end, color='black', style='->', linewidth=3, connectionstyle="arc3,rad=0"):
    """Draw an arrow between two points"""
    arrow = FancyArrowPatch(start, end,
                           arrowstyle=style,
                           mutation_scale=30,
                           color=color,
                           linewidth=linewidth,
                           connectionstyle=connectionstyle)
    ax.add_patch(arrow)


def add_pathway_label(ax, x, y, text, fontsize=pathway_fontsize, color='#333333'):
    """Add pathway label WITHOUT box, with color matching arrow"""
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            color=color, fontweight='bold', style='italic')


# ============================================================
# LAYOUT - Spread out more to avoid overlap with larger fonts
# ============================================================

# Input (left side)
uva_x, uva_y = 12, 48

# Main state variables (spread out more)
lai_x, lai_y = 42, 64
cbuf_x, cbuf_y = 68, 64
xd_x, xd_y = 94, 64

ros_x, ros_y = 42, 36
stress_x, stress_y = 68, 36
anth_x, anth_y = 94, 36

# Damage mechanisms (bottom) - spread out more
vuln_x, vuln_y = 28, 10
circadian_x, circadian_y = 54, 10
nonlin_x, nonlin_y = 82, 10

# ============================================================
# DRAW STATE BOXES - Larger boxes for larger fonts
# ============================================================

# Input
draw_state_box(ax, uva_x, uva_y, 18, 14, 'UV-A', '(I_UVA)', colors['input'])

# Growth states
draw_state_box(ax, lai_x, lai_y, 16, 12, 'LAI', '', colors['growth'])
draw_state_box(ax, cbuf_x, cbuf_y, 16, 12, 'C_buf', '', colors['growth'])
draw_state_box(ax, xd_x, xd_y, 18, 12, 'X_d', '(Biomass)', colors['growth'])

# Stress states
draw_state_box(ax, ros_x, ros_y, 16, 12, 'ROS', '', colors['stress'])
draw_state_box(ax, stress_x, stress_y, 16, 12, 'Stress', '', colors['stress'])

# Anthocyanin
draw_state_box(ax, anth_x, anth_y, 16, 12, 'Anth', '', colors['anth'])

# Damage mechanism boxes
draw_state_box(ax, vuln_x, vuln_y, 20, 12, 'LAI', 'Vulnerability', colors['damage'])
draw_state_box(ax, circadian_x, circadian_y, 18, 12, 'Circadian', 'Damage', colors['damage'])
draw_state_box(ax, nonlin_x, nonlin_y, 20, 12, 'Nonlinear', 'Damage', colors['damage'])

# ============================================================
# DRAW ARROWS WITH PATHWAY LABELS (colors match arrows)
# ============================================================

# UV-A → LAI (Morphological effect) - curved up - POSITIVE
draw_arrow(ax, (uva_x + 9, uva_y + 5), (lai_x - 8, lai_y - 3),
           colors['arrow_positive'], linewidth=4, connectionstyle="arc3,rad=-0.3")
add_pathway_label(ax, 20, 64, 'Morphological\nEffect\n(SLA↑, LAI↑)', color=colors['arrow_positive'])

# UV-A → ROS (ROS production) - NEGATIVE
draw_arrow(ax, (uva_x + 9, uva_y - 5), (ros_x - 8, ros_y + 3),
           colors['arrow_negative'], linewidth=4, connectionstyle="arc3,rad=0.2")
add_pathway_label(ax, 20, 38, 'ROS Production\n(k_ros × I_UVA)', color=colors['arrow_negative'])

# LAI → C_buf (Light interception) - POSITIVE
draw_arrow(ax, (lai_x + 8, lai_y), (cbuf_x - 8, cbuf_y),
           colors['arrow_positive'], linewidth=4)
add_pathway_label(ax, 55, 72, 'Light\nInterception', color=colors['arrow_positive'])

# C_buf → X_d (Carbon allocation) - POSITIVE
draw_arrow(ax, (cbuf_x + 8, cbuf_y), (xd_x - 9, xd_y),
           colors['arrow_positive'], linewidth=4)
add_pathway_label(ax, 81, 72, 'Carbon\nAllocation', color=colors['arrow_positive'])

# ROS → Stress (Damage accumulation) - NEGATIVE
draw_arrow(ax, (ros_x + 8, ros_y), (stress_x - 8, stress_y),
           colors['arrow_negative'], linewidth=4)
add_pathway_label(ax, 55, 42, 'Damage\nAccumulation', color=colors['arrow_negative'])

# Stress → C_buf (Growth inhibition) - curved - NEGATIVE
draw_arrow(ax, (stress_x - 2, stress_y + 6), (cbuf_x - 2, cbuf_y - 6),
           colors['arrow_negative'], linewidth=4, connectionstyle="arc3,rad=0.4")
add_pathway_label(ax, 60, 50, 'Growth\nInhibition (−)', color=colors['arrow_negative'])

# Stress → Anth (Stress-induced synthesis) - POSITIVE
draw_arrow(ax, (stress_x + 8, stress_y), (anth_x - 8, anth_y),
           colors['arrow_positive'], linewidth=4)
add_pathway_label(ax, 81, 42, 'Stress-Induced\nSynthesis', color=colors['arrow_positive'])

# Anth → X_d (Antioxidant protection) - curved - POSITIVE (protects growth)
draw_arrow(ax, (anth_x, anth_y + 6), (xd_x, xd_y - 6),
           colors['arrow_positive'], linewidth=4, connectionstyle="arc3,rad=-0.4")
add_pathway_label(ax, 100, 50, 'Antioxidant\nProtection (−)', color=colors['arrow_positive'])

# Stress → LAI (Inhibits morphological effect) - curved - NEGATIVE
draw_arrow(ax, (stress_x - 6, stress_y + 6), (lai_x + 6, lai_y - 6),
           colors['arrow_negative'], linewidth=4, connectionstyle="arc3,rad=-0.3")
add_pathway_label(ax, 48, 52, 'Inhibits Morph.\nEffect (−)', color=colors['arrow_negative'])

# Stress decay - arrow from Stress going OUT and curving back (self-loop outside) - NEUTRAL
# Draw as curved arrow going right and down then pointing back
draw_arrow(ax, (stress_x + 8, stress_y - 3), (stress_x + 8, stress_y + 3),
           colors['arrow_neutral'], linewidth=3, connectionstyle="arc3,rad=-1.2")
add_pathway_label(ax, 80, 32, 'decay', color=colors['arrow_neutral'])

# ============================================================
# DAMAGE MECHANISMS CONNECTIONS
# ============================================================

# LAI → LAI Vulnerability (feedback line)
draw_arrow(ax, (lai_x - 6, lai_y - 6), (vuln_x + 6, vuln_y + 6),
           colors['arrow_neutral'], linewidth=3, connectionstyle="arc3,rad=0.2")
add_pathway_label(ax, 28, 40, 'v(LAI) =\nA·exp(−k·LAI)+1', color=colors['arrow_neutral'])

# LAI Vulnerability → ROS (modulates damage)
draw_arrow(ax, (vuln_x + 6, vuln_y + 6), (ros_x - 6, ros_y - 6),
           colors['arrow_neutral'], linewidth=3, connectionstyle="arc3,rad=-0.2")
add_pathway_label(ax, 38, 20, 'k1·ROS·v(LAI)\n(D12 groups)', color=colors['arrow_neutral'])

# Circadian → Stress - NEGATIVE
draw_arrow(ax, (circadian_x, circadian_y + 6), (stress_x - 4, stress_y - 6),
           colors['arrow_negative'], linewidth=3, connectionstyle="arc3,rad=0.15")
add_pathway_label(ax, 54, 24, 'Night UVA\n(L6D6-N)', color=colors['arrow_negative'])

# Nonlinear → Stress - NEGATIVE
draw_arrow(ax, (nonlin_x, nonlin_y + 6), (stress_x + 4, stress_y - 6),
           colors['arrow_negative'], linewidth=3, connectionstyle="arc3,rad=-0.15")
add_pathway_label(ax, 78, 24, 'k2·ROS·nonlin\n(H12D3)', color=colors['arrow_negative'])

# Gompertz formula near Nonlinear Damage box
add_pathway_label(ax, 95, 10, 'Gompertz\n(threshold=10.5h)', color=colors['arrow_neutral'])

# ============================================================
# TITLE AND LEGEND
# ============================================================

ax.set_title('Six-State ODE Model: System Dynamics Block Diagram',
             fontsize=title_fontsize, fontweight='bold', pad=40)

# Legend - positioned at bottom left, horizontal layout
legend_x, legend_y = 3, 2
legend_items = [
    (colors['growth'], 'Growth States (LAI, C_buf, X_d)'),
    (colors['stress'], 'Stress States (ROS, Stress)'),
    (colors['anth'], 'Anthocyanin (Anth)'),
    (colors['damage'], 'Damage Mechanisms'),
]

for i, (color, text) in enumerate(legend_items):
    rect = patches.Rectangle((legend_x + i * 28, legend_y), 4, 3,
                             facecolor=color, edgecolor='black', linewidth=2)
    ax.add_patch(rect)
    ax.text(legend_x + i * 28 + 5, legend_y + 1.5, text,
            fontsize=legend_fontsize, va='center')

# Save
plt.tight_layout()
plt.savefig('paper_figures/Fig17_system_block_diagram.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig('paper_figures/Fig17_system_block_diagram.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("Fig 17 generated:")
print("  - Box labels: 40pt")
print("  - Pathway labels: 26pt (colors match arrows)")
print("  - Title: 52pt")
print("  - Figure size: 40x28 inches")
print("  - Decay arrow moved outside Stress box")
print("Saved to paper_figures/Fig17_system_block_diagram.png/pdf")
