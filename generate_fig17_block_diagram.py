#!/usr/bin/env python3
"""
Generate Fig 17: System Block Diagram for Six-State ODE Model
- Smaller figure for paper
- Smaller boxes, same text size
- Fix overlap at bottom
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Create figure - paper size (single column width ~3.5in, double ~7in)
fig, ax = plt.subplots(figsize=(14, 10))
ax.set_xlim(0, 140)
ax.set_ylim(0, 105)
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

# Font sizes - keep readable
title_fontsize = 16
box_fontsize = 14          # State box labels
box_detail_fontsize = 11   # Secondary text in boxes
pathway_fontsize = 10      # Pathway labels
legend_fontsize = 10

def draw_state_box(ax, x, y, width, height, label, sublabel, color, fontsize=box_fontsize):
    """Draw a state variable box - SMALLER"""
    box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.5",
                         facecolor=color, edgecolor='black', linewidth=1.5)
    ax.add_patch(box)

    if sublabel:
        ax.text(x, y + 2, label, ha='center', va='center', fontsize=fontsize, fontweight='bold')
        ax.text(x, y - 2.5, sublabel, ha='center', va='center', fontsize=box_detail_fontsize)
    else:
        ax.text(x, y, label, ha='center', va='center', fontsize=fontsize, fontweight='bold')


def draw_arrow(ax, start, end, color='black', style='->', linewidth=2, connectionstyle="arc3,rad=0"):
    """Draw an arrow between two points"""
    arrow = FancyArrowPatch(start, end,
                           arrowstyle=style,
                           mutation_scale=15,
                           color=color,
                           linewidth=linewidth,
                           connectionstyle=connectionstyle)
    ax.add_patch(arrow)


def add_pathway_label(ax, x, y, text, fontsize=pathway_fontsize, color='#333333'):
    """Add pathway label WITHOUT box, with color matching arrow"""
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            color=color, fontweight='bold', style='italic')


# ============================================================
# LAYOUT - Compact but spread out to avoid overlap
# ============================================================

# Input (left side)
uva_x, uva_y = 15, 60

# Main state variables - upper row
lai_x, lai_y = 50, 85
cbuf_x, cbuf_y = 85, 85
xd_x, xd_y = 120, 85

# Middle row - stress related
ros_x, ros_y = 50, 50
stress_x, stress_y = 85, 50
anth_x, anth_y = 120, 50

# Damage mechanisms (bottom) - spread out horizontally
vuln_x, vuln_y = 30, 15
circadian_x, circadian_y = 70, 15
nonlin_x, nonlin_y = 110, 15

# ============================================================
# DRAW STATE BOXES - Smaller boxes
# ============================================================

# Box dimensions - smaller
box_w, box_h = 18, 12
small_box_w, small_box_h = 20, 12

# Input
draw_state_box(ax, uva_x, uva_y, 16, 14, 'UV-A', '(I_UVA)', colors['input'])

# Growth states
draw_state_box(ax, lai_x, lai_y, box_w, box_h, 'LAI', '', colors['growth'])
draw_state_box(ax, cbuf_x, cbuf_y, box_w, box_h, 'C_buf', '', colors['growth'])
draw_state_box(ax, xd_x, xd_y, box_w, box_h, 'X_d', '(Biomass)', colors['growth'])

# Stress states
draw_state_box(ax, ros_x, ros_y, box_w, box_h, 'ROS', '', colors['stress'])
draw_state_box(ax, stress_x, stress_y, box_w, box_h, 'Stress', '', colors['stress'])

# Anthocyanin
draw_state_box(ax, anth_x, anth_y, box_w, box_h, 'Anth', '', colors['anth'])

# Damage mechanism boxes
draw_state_box(ax, vuln_x, vuln_y, small_box_w, box_h, 'LAI', 'Vulnerability', colors['damage'])
draw_state_box(ax, circadian_x, circadian_y, small_box_w, box_h, 'Circadian', 'Damage', colors['damage'])
draw_state_box(ax, nonlin_x, nonlin_y, small_box_w, box_h, 'Nonlinear', 'Damage', colors['damage'])

# ============================================================
# DRAW ARROWS WITH PATHWAY LABELS (colors match arrows)
# ============================================================

# UV-A → LAI (Morphological effect) - POSITIVE
draw_arrow(ax, (uva_x + 8, uva_y + 5), (lai_x - 9, lai_y - 4),
           colors['arrow_positive'], linewidth=2, connectionstyle="arc3,rad=-0.2")
add_pathway_label(ax, 25, 80, 'Morphological\nEffect\n(SLA↑, LAI↑)', color=colors['arrow_positive'])

# UV-A → ROS (ROS production) - NEGATIVE
draw_arrow(ax, (uva_x + 8, uva_y - 5), (ros_x - 9, ros_y + 4),
           colors['arrow_negative'], linewidth=2, connectionstyle="arc3,rad=0.15")
add_pathway_label(ax, 25, 52, 'ROS Production\n(k_ros × I_UVA)', color=colors['arrow_negative'])

# LAI → C_buf (Light interception) - POSITIVE
draw_arrow(ax, (lai_x + 9, lai_y), (cbuf_x - 9, cbuf_y),
           colors['arrow_positive'], linewidth=2)
add_pathway_label(ax, 67.5, 93, 'Light\nInterception', color=colors['arrow_positive'])

# C_buf → X_d (Carbon allocation) - POSITIVE
draw_arrow(ax, (cbuf_x + 9, cbuf_y), (xd_x - 9, xd_y),
           colors['arrow_positive'], linewidth=2)
add_pathway_label(ax, 102.5, 93, 'Carbon\nAllocation', color=colors['arrow_positive'])

# ROS → Stress (Damage accumulation) - NEGATIVE
draw_arrow(ax, (ros_x + 9, ros_y), (stress_x - 9, stress_y),
           colors['arrow_negative'], linewidth=2)
add_pathway_label(ax, 67.5, 57, 'Damage\nAccumulation', color=colors['arrow_negative'])

# Stress → C_buf (Growth inhibition) - curved - NEGATIVE
draw_arrow(ax, (stress_x - 3, stress_y + 6), (cbuf_x - 3, cbuf_y - 6),
           colors['arrow_negative'], linewidth=2, connectionstyle="arc3,rad=0.35")
add_pathway_label(ax, 76, 68, 'Growth\nInhibition (−)', color=colors['arrow_negative'])

# Stress → Anth (Stress-induced synthesis) - POSITIVE
draw_arrow(ax, (stress_x + 9, stress_y), (anth_x - 9, anth_y),
           colors['arrow_positive'], linewidth=2)
add_pathway_label(ax, 102.5, 57, 'Stress-Induced\nSynthesis', color=colors['arrow_positive'])

# Anth → X_d (Antioxidant protection) - curved - POSITIVE
draw_arrow(ax, (anth_x, anth_y + 6), (xd_x, xd_y - 6),
           colors['arrow_positive'], linewidth=2, connectionstyle="arc3,rad=-0.35")
add_pathway_label(ax, 128, 68, 'Antioxidant\nProtection (−)', color=colors['arrow_positive'])

# Stress → LAI (Inhibits morphological effect) - curved - NEGATIVE
draw_arrow(ax, (stress_x - 6, stress_y + 6), (lai_x + 6, lai_y - 6),
           colors['arrow_negative'], linewidth=2, connectionstyle="arc3,rad=-0.25")
add_pathway_label(ax, 60, 68, 'Inhibits Morph.\nEffect (−)', color=colors['arrow_negative'])

# Stress decay - arrow on right side going down and curving back
draw_arrow(ax, (stress_x + 9, stress_y - 2), (stress_x + 9, stress_y + 2),
           colors['arrow_neutral'], linewidth=1.5, connectionstyle="arc3,rad=-1.5")
add_pathway_label(ax, 98, 45, 'decay', color=colors['arrow_neutral'])

# ============================================================
# DAMAGE MECHANISMS CONNECTIONS - spread out to avoid overlap
# ============================================================

# LAI → LAI Vulnerability (feedback line)
draw_arrow(ax, (lai_x - 7, lai_y - 6), (vuln_x + 5, vuln_y + 6),
           colors['arrow_neutral'], linewidth=1.5, connectionstyle="arc3,rad=0.15")
add_pathway_label(ax, 32, 55, 'v(LAI) =\nA·exp(−k·LAI)+1', color=colors['arrow_neutral'])

# LAI Vulnerability → ROS (modulates damage)
draw_arrow(ax, (vuln_x + 8, vuln_y + 6), (ros_x - 7, ros_y - 6),
           colors['arrow_neutral'], linewidth=1.5, connectionstyle="arc3,rad=-0.15")
add_pathway_label(ax, 45, 30, 'k1·ROS·v(LAI)\n(D12 groups)', color=colors['arrow_neutral'])

# Circadian → Stress - NEGATIVE
draw_arrow(ax, (circadian_x, circadian_y + 6), (stress_x - 3, stress_y - 6),
           colors['arrow_negative'], linewidth=1.5, connectionstyle="arc3,rad=0.1")
add_pathway_label(ax, 70, 35, 'Night UVA\n(L6D6-N)', color=colors['arrow_negative'])

# Nonlinear → Stress - NEGATIVE
draw_arrow(ax, (nonlin_x - 5, nonlin_y + 6), (stress_x + 3, stress_y - 6),
           colors['arrow_negative'], linewidth=1.5, connectionstyle="arc3,rad=-0.1")
add_pathway_label(ax, 95, 30, 'k2·ROS·nonlin\n(H12D3)', color=colors['arrow_negative'])

# Gompertz formula - positioned right of Nonlinear box
add_pathway_label(ax, 130, 15, 'Gompertz\n(threshold=10.5h)', color=colors['arrow_neutral'])

# ============================================================
# TITLE AND LEGEND
# ============================================================

ax.set_title('Six-State ODE Model: System Dynamics Block Diagram',
             fontsize=title_fontsize, fontweight='bold', pad=15)

# Legend - horizontal at very bottom
legend_y = 2
legend_items = [
    (5, colors['growth'], 'Growth States (LAI, C_buf, X_d)'),
    (45, colors['stress'], 'Stress States (ROS, Stress)'),
    (80, colors['anth'], 'Anthocyanin (Anth)'),
    (110, colors['damage'], 'Damage Mechanisms'),
]

for lx, color, text in legend_items:
    rect = patches.Rectangle((lx, legend_y), 4, 3,
                             facecolor=color, edgecolor='black', linewidth=1)
    ax.add_patch(rect)
    ax.text(lx + 5, legend_y + 1.5, text, fontsize=legend_fontsize, va='center')

# Save
plt.tight_layout()
plt.savefig('paper_figures/Fig17_system_block_diagram.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig('paper_figures/Fig17_system_block_diagram.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("Fig 17 generated (paper size):")
print("  - Figure: 14x10 inches")
print("  - Box labels: 14pt")
print("  - Pathway labels: 10pt (colors match arrows)")
print("  - DPI: 300")
print("Saved to paper_figures/Fig17_system_block_diagram.png/pdf")
