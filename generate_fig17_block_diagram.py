#!/usr/bin/env python3
"""
Generate Fig 17: System Block Diagram for Six-State ODE Model
Font sizes significantly increased - pathway labels same as state boxes, state boxes even larger
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Create figure - extra large for better readability
fig, ax = plt.subplots(figsize=(36, 24))
ax.set_xlim(0, 100)
ax.set_ylim(0, 70)
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
title_fontsize = 48
box_fontsize = 36          # State box labels - very large
box_detail_fontsize = 24   # Secondary text in boxes
pathway_fontsize = 28      # Pathway labels - same as boxes
legend_fontsize = 24

def draw_state_box(ax, x, y, width, height, label, sublabel, color, fontsize=box_fontsize):
    """Draw a state variable box"""
    box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                         boxstyle="round,pad=0.03,rounding_size=0.8",
                         facecolor=color, edgecolor='black', linewidth=3)
    ax.add_patch(box)

    if sublabel:
        ax.text(x, y + 1.5, label, ha='center', va='center', fontsize=fontsize, fontweight='bold')
        ax.text(x, y - 2, sublabel, ha='center', va='center', fontsize=box_detail_fontsize)
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
    """Add pathway label WITHOUT box, with larger font"""
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            color=color, fontweight='bold', style='italic')


# ============================================================
# LAYOUT - Spread out more to avoid overlap with larger fonts
# ============================================================

# Input (left side)
uva_x, uva_y = 12, 42

# Main state variables (spread out more)
lai_x, lai_y = 38, 55
cbuf_x, cbuf_y = 62, 55
xd_x, xd_y = 86, 55

ros_x, ros_y = 38, 30
stress_x, stress_y = 62, 30
anth_x, anth_y = 86, 30

# Damage mechanisms (bottom)
vuln_x, vuln_y = 30, 8
circadian_x, circadian_y = 46, 8
nonlin_x, nonlin_y = 72, 8

# ============================================================
# DRAW STATE BOXES - Larger boxes for larger fonts
# ============================================================

# Input
draw_state_box(ax, uva_x, uva_y, 18, 12, 'UV-A', '(I_UVA)', colors['input'])

# Growth states
draw_state_box(ax, lai_x, lai_y, 14, 10, 'LAI', '', colors['growth'])
draw_state_box(ax, cbuf_x, cbuf_y, 14, 10, 'C_buf', '', colors['growth'])
draw_state_box(ax, xd_x, xd_y, 16, 10, 'X_d', '(Biomass)', colors['growth'])

# Stress states
draw_state_box(ax, ros_x, ros_y, 14, 10, 'ROS', '', colors['stress'])
draw_state_box(ax, stress_x, stress_y, 14, 10, 'Stress', '', colors['stress'])

# Anthocyanin
draw_state_box(ax, anth_x, anth_y, 14, 10, 'Anth', '', colors['anth'])

# Damage mechanism boxes
draw_state_box(ax, vuln_x, vuln_y, 18, 10, 'LAI', 'Vulnerability', colors['damage'])
draw_state_box(ax, circadian_x + 4, circadian_y, 16, 10, 'Circadian', 'Damage', colors['damage'])
draw_state_box(ax, nonlin_x, nonlin_y, 18, 10, 'Nonlinear', 'Damage', colors['damage'])

# ============================================================
# DRAW ARROWS WITH PATHWAY LABELS
# ============================================================

# UV-A → LAI (Morphological effect) - curved up
draw_arrow(ax, (uva_x + 9, uva_y + 4), (lai_x - 7, lai_y - 2),
           colors['arrow_positive'], linewidth=4, connectionstyle="arc3,rad=-0.3")
add_pathway_label(ax, 22, 54, 'Morphological\nEffect\n(SLA↑, LAI↑)')

# UV-A → ROS (ROS production)
draw_arrow(ax, (uva_x + 9, uva_y - 4), (ros_x - 7, ros_y + 2),
           colors['arrow_negative'], linewidth=4, connectionstyle="arc3,rad=0.2")
add_pathway_label(ax, 22, 34, 'ROS Production\n(k_ros × I_UVA)')

# LAI → C_buf (Light interception)
draw_arrow(ax, (lai_x + 7, lai_y), (cbuf_x - 7, cbuf_y),
           colors['arrow_positive'], linewidth=4)
add_pathway_label(ax, 50, 60, 'Light\nInterception')

# C_buf → X_d (Carbon allocation)
draw_arrow(ax, (cbuf_x + 7, cbuf_y), (xd_x - 8, xd_y),
           colors['arrow_positive'], linewidth=4)
add_pathway_label(ax, 74, 60, 'Carbon\nAllocation')

# ROS → Stress (Damage accumulation)
draw_arrow(ax, (ros_x + 7, ros_y), (stress_x - 7, stress_y),
           colors['arrow_negative'], linewidth=4)
add_pathway_label(ax, 50, 35, 'Damage\nAccumulation')

# Stress → C_buf (Growth inhibition) - curved
draw_arrow(ax, (stress_x, stress_y + 5), (cbuf_x, cbuf_y - 5),
           colors['arrow_negative'], linewidth=4, connectionstyle="arc3,rad=0.3")
add_pathway_label(ax, 56, 44, 'Growth\nInhibition (−)')

# Stress → Anth (Stress-induced synthesis)
draw_arrow(ax, (stress_x + 7, stress_y), (anth_x - 7, anth_y),
           colors['arrow_positive'], linewidth=4)
add_pathway_label(ax, 74, 35, 'Stress-Induced\nSynthesis')

# Anth → X_d (Antioxidant protection) - curved
draw_arrow(ax, (anth_x, anth_y + 5), (xd_x, xd_y - 5),
           colors['arrow_positive'], linewidth=4, connectionstyle="arc3,rad=-0.3")
add_pathway_label(ax, 92, 44, 'Antioxidant\nProtection (−)')

# Stress → LAI (Inhibits morphological effect) - curved
draw_arrow(ax, (stress_x - 4, stress_y + 5), (lai_x + 4, lai_y - 5),
           colors['arrow_negative'], linewidth=4, connectionstyle="arc3,rad=-0.3")
add_pathway_label(ax, 44, 44, 'Inhibits Morph.\nEffect (−)')

# Stress decay (self-loop arrow representation)
draw_arrow(ax, (stress_x + 5, stress_y - 5), (stress_x + 7, stress_y - 3),
           colors['arrow_neutral'], linewidth=3, connectionstyle="arc3,rad=-1.5")
add_pathway_label(ax, 72, 22, 'decay')

# ============================================================
# DAMAGE MECHANISMS CONNECTIONS
# ============================================================

# LAI Vulnerability → ROS
draw_arrow(ax, (vuln_x, vuln_y + 5), (ros_x - 4, ros_y - 5),
           colors['arrow_neutral'], linewidth=3, connectionstyle="arc3,rad=0.2")
add_pathway_label(ax, 28, 18, 'v(LAI) =\nA·exp(−k·LAI)+1')

# ROS → LAI Vulnerability (feedback)
draw_arrow(ax, (ros_x - 4, ros_y - 5), (vuln_x + 4, vuln_y + 5),
           colors['arrow_neutral'], linewidth=3, connectionstyle="arc3,rad=0.2")
add_pathway_label(ax, 38, 18, 'k1·ROS·v(LAI)\n(D12 groups)')

# Circadian → Stress
draw_arrow(ax, (circadian_x + 4, circadian_y + 5), (stress_x - 2, stress_y - 5),
           colors['arrow_negative'], linewidth=3, connectionstyle="arc3,rad=0.2")
add_pathway_label(ax, 52, 18, 'Night UVA\n(L6D6-N)')

# Nonlinear → Stress
draw_arrow(ax, (nonlin_x, nonlin_y + 5), (stress_x + 4, stress_y - 5),
           colors['arrow_negative'], linewidth=3, connectionstyle="arc3,rad=-0.2")
add_pathway_label(ax, 68, 18, 'k2·ROS·nonlin\n(H12D3)')

# Gompertz formula
add_pathway_label(ax, 80, 18, 'Gompertz\n(threshold = 10.5h)')

# ============================================================
# TITLE AND LEGEND
# ============================================================

ax.set_title('Six-State ODE Model: System Dynamics Block Diagram',
             fontsize=title_fontsize, fontweight='bold', pad=30)

# Legend - positioned at bottom left
legend_x, legend_y = 5, 3
legend_items = [
    (colors['growth'], 'Growth States (LAI, C_buf, X_d)'),
    (colors['stress'], 'Stress States (ROS, Stress)'),
    (colors['anth'], 'Anthocyanin (Anth)'),
    (colors['damage'], 'Damage Mechanisms'),
]

for i, (color, text) in enumerate(legend_items):
    rect = patches.Rectangle((legend_x, legend_y + i * 4), 4, 3,
                             facecolor=color, edgecolor='black', linewidth=2)
    ax.add_patch(rect)
    ax.text(legend_x + 6, legend_y + i * 4 + 1.5, text,
            fontsize=legend_fontsize, va='center')

# Save
plt.tight_layout()
plt.savefig('paper_figures/Fig17_system_block_diagram.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig('paper_figures/Fig17_system_block_diagram.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("Fig 17 generated with MUCH larger fonts:")
print("  - Box labels: 36pt")
print("  - Pathway labels: 28pt")
print("  - Title: 48pt")
print("  - Figure size: 36x24 inches")
print("Saved to paper_figures/Fig17_system_block_diagram.png/pdf")
