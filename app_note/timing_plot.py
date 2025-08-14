import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def time_to_minutes(time_str):
    """Convert time string like '1m23.456s' to minutes as float"""
    # Remove 's' at the end and split by 'm'
    time_str = time_str.rstrip('s')
    if 'm' in time_str:
        parts = time_str.split('m')
        minutes = float(parts[0])
        seconds = float(parts[1]) if parts[1] else 0
    else:
        # Just seconds
        minutes = 0
        seconds = float(time_str)
    return minutes + seconds / 60.0


# Read the data
df = pd.read_csv("run_time.tsv", sep='\t')

# Convert time columns to minutes
virus_columns = df.columns[2:].array.tolist()
for col in virus_columns:
    df[col + '_minutes'] = df[col].apply(time_to_minutes)

# Reshape data for plotting
plot_data = []
for _, row in df.iterrows():
    for virus in virus_columns:
        plot_data.append({
            'platform': row['platform'],
            'virus': virus,
            'run': row['i'],
            'time_minutes': row[virus + '_minutes']
        })

plot_df = pd.DataFrame(plot_data)

# Create the plot
fig, ax = plt.subplots(figsize=(14, 8))

# Define colors for each platform
platform_colors = {
    'linux': '#1f77b4',      # blue
    'macX86': '#ff7f0e',     # orange
    'macArm': '#2ca02c'      # green
}

# Define markers for each platform
platform_markers = {
    'linux': 'o',
    'macX86': 's',
    'macArm': '^'
}

# Create x-axis positions.  Reorder viruses in order of speed
viruses = ['oropouche_m', 'mumps', 'measles', 'RSV-A', 'mpox']
platforms = ['linux', 'macX86', 'macArm']

# Width for grouping platforms within each virus
platform_width = 1 / (len(platforms) + 1)
virus_positions = np.arange(len(viruses))

# Plot points for each platform-virus combination
for i, platform in enumerate(platforms):
    platform_data = plot_df[plot_df['platform'] == platform]

    for j, virus in enumerate(viruses):
        virus_data = platform_data[platform_data['virus'] == virus]

        # X position: virus position + platform offset
        x_pos = virus_positions[j] + (i - 1) * platform_width

        # Plot individual run points with some jitter for visibility
        x_jitter = np.random.normal(0, 0.02, len(virus_data))  # small random offset
        ax.scatter(
            [x_pos] * len(virus_data) + x_jitter,
            virus_data['time_minutes'],
            color=platform_colors[platform],
            marker=platform_markers[platform],
            s=60,
            alpha=0.7,
            label=platform if j == 0 else ""  # Only label once per platform
        )

# Customize the plot
ax.set_xlabel('Viral Species', fontsize=12, fontweight='bold')
ax.set_ylabel('Runtime (minutes)', fontsize=12, fontweight='bold')
ax.set_title('Performance Comparison Across Viral Species and Platforms\n(Individual Run Times)',
             fontsize=14, fontweight='bold', pad=20)

# Set x-axis labels and positions
virus_labels = ['Oropouche', 'Mumps', 'Measles', 'RSV-A', 'Mpox']
ax.set_xticks(virus_positions)
ax.set_xticklabels(virus_labels)

# Add legend
ax.legend(title='Platform', loc='upper left', frameon=True, fancybox=True, shadow=True)

# Add grid for better readability
ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
ax.set_axisbelow(True)

# Set y-axis to start from 0 for better comparison
ax.set_ylim(bottom=0, top=plot_df['time_minutes'].max() * 1.05)

# Add minor ticks for better precision reading
ax.minorticks_on()

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save the figure as a PNG file
plt.savefig("timing_plot.png")

# Show the plot
plt.show()

# Print some summary statistics
print("\nSummary Statistics:")
print("=" * 50)
for platform in platforms:
    print(f"\n{platform.upper()}:")
    platform_data = plot_df[plot_df['platform'] == platform]
    for virus in viruses:
        virus_data = platform_data[platform_data['virus'] == virus]['time_minutes']
        print(f"  {virus:12s}: {virus_data.min():.2f}m, {virus_data.median():.2f}m, {virus_data.max():.2f}m")
