import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
# Use LaTeX for all text rendering (brings in Computer Modern automatically)
rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
rcParams["font.serif"] = ["Computer Modern"]

# Optional for clean math mode rendering
rcParams["mathtext.fontset"] = "cm"


# === Increase font sizes ===
rcParams["font.size"] = 16         # base font size
rcParams["axes.titlesize"] = 20    # title
rcParams["axes.labelsize"] = 18    # x and y labels
rcParams["xtick.labelsize"] = 14   # x tick labels
rcParams["ytick.labelsize"] = 14   # y tick labels
rcParams["legend.fontsize"] = 14   # legend text
rcParams["figure.titlesize"] = 16  # figure title

def plot_DAM_flex_nonflex_bar_LCOH(scenarios, revenue_categories, expense_categories,
                                            revenue_data_input, expense_data_input, LCOH_data, npv_data, # Added LCOH_data
                                            ylabel=r'Value $(\$MM$)', # Renamed to revenue/expense label
                                            lcoh_ylabel=r'LCOH ($\$/kg$)'): # New label for LCOH axis
    scenarios = np.array(scenarios)
    revenue_data = np.array(revenue_data_input)
    expense_data = np.array(expense_data_input).T # Transpose as per your original example
    
    n = len(scenarios)
    x = np.arange(n)
    width = 0.45

    fig, ax1 = plt.subplots(figsize=(15, 7)) # Renamed ax to ax1 for clarity

    # --- ADD THIS BLOCK FOR SECOND Y-AXIS ---
    ax2 = ax1.twinx() # Create a second Axes that shares the same x-axis
    # --- END ADD BLOCK ---

    # Plot revenue
    bottom = np.zeros(n)
    revenue_hatches = [None]
    for i, cat in enumerate(revenue_categories):
        # Using specific colors to keep it visually appealing
        # color = plt.cm.Greens(0.4 + i * 0.2) # Example: shades of green
        # ax1.bar(x - width/2, revenue_data[:, i], width, bottom=bottom, label=cat, color=color, edgecolor='forestgreen')
        
        color = plt.cm.Blues(0.4 + i * 0.2) # CHANGED to Blues
        hatch = revenue_hatches[i % len(revenue_hatches)] # Get the hatch
        
        ax1.bar(x - width/2, revenue_data[:, i], width, bottom=bottom, label=cat, 
                color=color, edgecolor='darkblue', hatch=hatch) # ADDED hatch=hatch
        
        bottom += revenue_data[:, i]

    # Plot expense (as negative)
    expense_hatches = ['//', 'xx', '..', '+', 'o']
    bottom = np.zeros(n)
    expense_colors = [
        '#FFB300',  # Amber (for 'Plant Capital')
        '#E65100',  # Dark Orange (for 'Plant OPEX')
        '#D48265',
        '#8D6E63',   # Brown (for 'Electrolyzer Stack Replacement')
        '#A1887F',   # A lighter, grayer-brown
    ]
    for i, cat in enumerate(expense_categories):
        # Using specific colors for expenses
        # color = plt.cm.Reds(0.4 + i * 0.2) # Example: shades of red
        # ax1.bar(x - width/2, -expense_data[:, i], width, bottom=-bottom, label=cat, color=color, edgecolor='maroon')
        
        # --- NEW ---
        # color = plt.cm.Oranges(0.4 + i * 0.2) # CHANGED to Oranges
        color = expense_colors[i % len(expense_colors)] # Get color from our new list
        hatch = expense_hatches[i % len(expense_hatches)] # Get the hatch
        
        ax1.bar(x - width/2, -expense_data[:, i], width, bottom=-bottom, label=cat, 
                color=color, edgecolor='black', hatch=hatch) # ADDED hatch=hatch
        
        bottom += expense_data[:, i]

    # Profit bar
    # total_revenue = revenue_data[0] # This was an error in your original code if multiple categories
    # total_revenue = revenue_data.sum(axis=1) # Sums up all revenue categories for each scenario
    # # total_expense = -expense_data.sum(axis=1) # The negative sign is for plotting, sum absolute
    # total_expense = expense_data.sum(axis=1) # Sum of absolute expense categories for each scenario
    # profit = total_revenue - total_expense # Profit is revenue minus absolute expenses

    profit_bars = ax1.bar(x + width/2, npv_data.flatten(), width, label='Net Present Value', color='black', edgecolor='black')

    for bar, val in zip(profit_bars, npv_data.flatten()):
        # Adjusted text position for clarity and LaTeX formatting
        # Using a fixed offset for simplicity as requested, but a dynamic one is often better
        text_offset = 0.35 # Adjust this value if labels overlap
        if val >= 0:
            text_pos = bar.get_height() + text_offset 
            va = 'bottom'
        else:
            text_pos = bar.get_height() - text_offset -0.2
            va = 'top'
        ax1.text(bar.get_x() + bar.get_width()/2, text_pos, f"\\${val:1.2f}MM",
        ha='center', va=va, fontsize=12) # usetex=True is set globally
        
    # Formatting for ax1 (primary y-axis)
    ax1.axhline(0, color='black', linewidth=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(scenarios)
    ax1.set_ylabel(ylabel)
    ax1.set_yticks(np.linspace(-30, 30, num = 13))

    # --- ADD THIS BLOCK FOR LCOH PLOT AND AXIS ---
    lcoh_color = '#007FFF' # A distinct blue for LCOH
    lcoh_marker = 'o'
    lcoh_markersize = 8
    lcoh_linestyle = 'None'

    lcoh_line, = ax2.plot(x, LCOH_data, marker=lcoh_marker, markersize=lcoh_markersize,
                          color=lcoh_color, linestyle=lcoh_linestyle, label='LCOH')
    # --- Second Y-axis Formatting (LCOH Axis) ---
    ax2.set_ylabel(lcoh_ylabel, fontsize=rcParams["axes.labelsize"], color=lcoh_line.get_color())
    ax2.tick_params(axis='y', labelcolor=lcoh_line.get_color())

    # --- Start of lines to change ---
    lcoh_min = np.min(LCOH_data)
    lcoh_max = np.max(LCOH_data)
    
    # Calculate a buffer based on the range of LCOH data
    lcoh_range = lcoh_max - lcoh_min

    lower_bound = max(0.0, lcoh_min - lcoh_range * 0.1) # 10% buffer below min, but not < 0
    upper_bound = lcoh_max + lcoh_range * 0.1 # 10% buffer above max
    
    ax2.set_ylim(0, 6)
    # --- End of lines to change ---

    def lcoh_formatter(x, pos):
        return f"{x:.2f}"
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lcoh_formatter))
    # --- END ADD BLOCK ---
    # --- LCOH PLOT AND LABELS ---
    # lcoh_color = '#007FFF'
    
    # # Flatten the data here to avoid the TypeError
    lcoh_flat = np.array(LCOH_data).flatten()
    
    # lcoh_line, = ax2.plot(x, lcoh_flat, marker='o', markersize=8,
    #                       color=lcoh_color, 
    #                       #linestyle='--', 
    #                       label='LCOH')

    # Add text labels above each dot using the flattened data
    # y_offset = 0.10 
    # for i, val in enumerate(lcoh_flat):
    #     # Now 'val' is a float, so :.2f will work perfectly
    #     ax2.text(x[i], val + y_offset, f'\${val:.2f}', 
    #              ha='center', va='bottom', color=lcoh_color, 
    #              fontsize=11, fontweight='bold')
    # Add text labels with a white background box
    y_offset = 0.35
    for i, val in enumerate(lcoh_flat):
        ax2.text(x[i], val - y_offset, rf'\textbf{{\Large \$}}{val:.2f}', 
                 ha='center', va='bottom', color=lcoh_color, 
                 fontsize=16, fontweight='bold',
                 # The 'bbox' adds the white background
                 bbox=dict(facecolor='white', 
                           edgecolor=lcoh_color, 
                           boxstyle='round,pad=0.3', 
                           alpha=0.9))



    # --- Legend (COMBINED from both axes) ---
    # Collect handles and labels for a single legend
    # Get handles from ax1
    handles1, labels1 = ax1.get_legend_handles_labels()
    
    # Add LCOH line handle
    lcoh_handle = lcoh_line # Use the Line2D object directly
    
    # Combine handles and labels
    all_handles = handles1 + [lcoh_handle]
    all_labels = labels1 + ['LCOH'] # Add LCOH label
    # all_handles = handles1
    # all_labels = labels1  # Add LCOH label

    ax1.legend(handles=all_handles, labels=all_labels, loc='best', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    # --- END LEGEND ---

    plt.subplots_adjust(bottom=0.15, right=0.78, top=0.95, left=0.07)
# Example: Adding a "Market Type" category below the scenarios
    group_labels = ['ERCOT Average', 'Panhandle Node', 'Distribution Pricing']
    group_positions = [0.5, 2.5, 4] # X-coordinates between the bars

    for pos, label in zip(group_positions, group_labels):
        ax1.text(pos, -0.09, label, transform=ax1.get_xaxis_transform(),
                ha='center', va='top', fontsize=18, fontweight='bold',
                bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round,pad=0.3'))



    plt.savefig("plots/figures/stacked_bar_chart_with_lcoh.pdf", format="pdf")
    plt.show()

# FLEX V. NONFLEX:
# -------------------------------
scenarios = ['Flexible', 'Non-Flexible', 'Flexible', 'Non-Flexible', 'Non-Flexible']
revenue_categories = ['Hydrogen Revenue']
expense_categories = ['Plant CAPEX', 'Plant OPEX', 'Stack Replacement', 'Electricity Expenses']
revenue_data = np.array([
[2.553647E7],[2.9829243905890103e7], [2.5619466E7], [2.9829243905890103e7], [2.9829243905890103e7]
])/1000000
expense_data = np.array([
    [3.9952E6,3.9952E6, 3.9952E6, 3.9952E6, 3.9952E6], #capex
    [103479*22,103479*22,103479*22,103479*22,103479*22], #opex
    [550000*2, 550000*3, 550000*2, 550000*3,550000*3], # replacement
    [1.2818e7,	1.9145e7, 1.1959e7, 1.8144e7, 3.059430395994904e7] #electricity
])/1000000
npv_data = np.array([
    [8.6085e6],[5.910e6], # average market (flex, non)
    [9.175e6],[6.485e6], # wind market
    [-1.6301e6]
])/1000000

LCOH_data = np.array([
    [3.44],[3.86], # average market (flex, non)
    [3.34],[3.77], # wind market
    [5.10]
])



# plot_DAM_flex_nonflex_bar_LCOH(scenarios, revenue_categories, expense_categories, revenue_data, expense_data, LCOH_data, npv_data)

def create_comparison_plot(
    profiles,
    categories,
    legend_labels=['Hydrogen Revenue', 'Electricity Expenses'],
    colors=[plt.cm.Blues(0.4), '#A1887F'],
    x_label=r'Value $(\$MM$)',
    figsize=(10, 8)
):
    """
    Creates a multi-subplot comparison plot with grouped horizontal bars,
    baselines, and custom annotations.
    """
    n_profiles = len(profiles)
    fig, axes = plt.subplots(n_profiles, 1, figsize=figsize)
    if n_profiles == 1: axes = [axes]

    bar_width = 0.4
    y_pos = np.arange(len(categories))

    for ax, profile in zip(axes, profiles):
        # 1. Plot Grouped Bars
        ax.barh(y_pos - bar_width/2, profile['values_top'], bar_width,
                label=legend_labels[0], color=colors[0], edgecolor='black')
        ax.barh(y_pos + bar_width/2, profile['values_bottom'], bar_width,
                label=legend_labels[1], color=colors[1], edgecolor='black', hatch='+')

        # 2. Vertical Baselines
        delta = -0.5
        ax.vlines(x=[profile['baseline_top'], profile['baseline_bottom']],
                  ymin=delta+0.1, ymax=len(categories)+delta-0.1, colors='black', linewidth=3)

        # 3. Flexible Annotations
        for ann in profile.get('annotations', []):
            bar_width=0.4
            y_coord = y_pos[ann['y_idx']] + (bar_width/2+0.1 if ann['bar_type'] == 'bottom' else -bar_width/2)
            # x_ref can be the bar end or the baseline
            x_ref = ann.get('x_ref', profile['baseline_top'] if ann['bar_type'] == 'top' else profile['baseline_bottom'])
            ax.text(x_ref + ann.get('offset', 0), y_coord, ann['text'],
                    va=ann.get('va', 'center'), ha=ann.get('ha', 'left'), 
                    fontweight='bold', fontsize=16)

        # 4. Aesthetic Formatting
        ax.set_yticks(y_pos)
        ax.set_yticklabels(categories)
        ax.invert_yaxis()
        ax.set_title(profile['title'], fontweight='bold', pad=15)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.grid(axis='x', linestyle='-', alpha=0.3)

    axes[-1].set_xlabel(x_label)
    
    # Legend - Centered between subplots
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc='center', bbox_to_anchor=(0.55, 0.5),
               ncol=2, frameon=True, edgecolor='black')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if n_profiles > 1: plt.subplots_adjust(hspace=0.6)
    return fig

# Define your y-axis categories
categories = ['Flexible', 'Non-Flexible']

# Define individual profile data
p1_hyd_rev = np.array([2.553647E7, 2.9829e7]) *1e-6
p1_elec_exp = np.array([1.2818e7,	1.9145e7])*1e-6
rel_dif1 = (p1_elec_exp[1]-p1_elec_exp[0])/p1_elec_exp[1]*100
eps = 0.01
p1_data = {
    'title': "ERCOT Hub Average Price Profile",
    'values_top': p1_hyd_rev,       # Bar lengths for 'hydrogen profit'
    'values_bottom': p1_elec_exp,    # Bar lengths for 'Electricity expenses'
    'baseline_top': max(p1_hyd_rev)+eps,                    # Vertical line for top bars
    'baseline_bottom': max(p1_elec_exp),                 # Vertical line for bottom bars
    'annotations': [
        {'text': f'-{rel_dif1:.1f}\%', 'y_idx': 0, 'bar_type': 'bottom', 'x_ref': 15.8, 'ha': 'right', 'va': 'bottom'},
        {'text': '-14.4\%', 'y_idx': 0, 'bar_type': 'top','x_ref': 26, 'ha': 'left', },
        # ... add more labels as needed
    ]
}

p2_hyd_rev = np.array([2.5619466E7, 2.9829e7])*1e-6
p2_elec_exp = np.array([1.1959e7, 1.8144e7])*1e-6
rel_dif2 = (p2_elec_exp[1]-p2_elec_exp[0])/p2_elec_exp[1]*100
p2_data = {
    'title': "ERCOT Panhandle Price Profile",
    'values_top': p2_hyd_rev,       # Bar lengths for 'hydrogen profit'
    'values_bottom': p2_elec_exp,    # Bar lengths for 'Electricity expenses'
    'baseline_top': max(p2_hyd_rev),                    # Vertical line for top bars
    'baseline_bottom': max(p2_elec_exp),                 # Vertical line for bottom bars
    'annotations': [
        {'text': f'-{rel_dif2:.1f}\%', 'y_idx': 0, 'bar_type': 'bottom', 'x_ref': 15.1, 'ha': 'right', 'va': 'bottom'},
        {'text': '-14.1\%', 'y_idx': 0, 'bar_type': 'top','x_ref': 26, 'ha': 'left', },
        # ... add more labels as needed
    ]
}



# Generate the plot
fig = create_comparison_plot([p1_data, p2_data], categories)
# plt.show()

plt.savefig("plots/figures/ercot_case_study.pdf", format="pdf")