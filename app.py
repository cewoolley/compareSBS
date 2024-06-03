from flask import Flask, render_template, request, jsonify
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Set the backend to 'Agg' else crashes on mac...
import matplotlib.pyplot as plt
import seaborn as sns
import re
import io
import base64
from scipy.spatial.distance import cosine

app = Flask(__name__)

# Load the data from the TSV file
data = pd.read_csv('COSMIC_v3.4_SBS_GRCh38.tsv', sep='\t', index_col=0)

def subtract_signatures(primary_sig, secondary_sig):
    subtracted_sig = data[primary_sig] - data[secondary_sig]
    return subtracted_sig

def plot_signature(signature, primary_sig, secondary_sig, title):
    # Create a DataFrame from the signature Series
    df = pd.DataFrame({'Proportion': signature})
    # Extract mutation types from the index using regular expressions
    df['MutationType'] = df.index.str.extract(r'\[(\w>\w)\]', expand=False)
    # Define colour mapping for mutation types
    mutation_colors = {
        'C>T': '#FF0000',
        'C>G': '#000000',
        'T>A': '#808080',
        'T>G': '#FFC0CB',
        'T>C': '#00FF00',
        'C>A': '#0000FF'
    }
    # Sort the DataFrame by 'MutationType'
    df['MutationType'] = pd.Categorical(df['MutationType'], categories=['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'])
    df = df.sort_values('MutationType')

    # Calculate cosine similarity
    primary_sig_values = data[primary_sig].values
    secondary_sig_values = data[secondary_sig].values
    cosine_sim = 1 - cosine(primary_sig_values, secondary_sig_values)

    # Create a new figure
    fig, ax = plt.subplots(figsize=(20, 6))
    # Create a grouped bar plot
    sns.barplot(x=df.index, y='Proportion', data=df, hue='MutationType', dodge=True, palette=mutation_a, ax=ax)

    ax.set_xlabel('Trinucleotide Context')
    ax.set_ylabel('Proportion of Single Base Substitutions')
    ax.set_title(f"{title}\nCosine Similarity: {cosine_sim:.4f}")
    ax.legend(title='Mutation Types', loc='upper right')
    plt.xticks(rotation=90)
    plt.xticks(fontsize=8)

    # Set background colours above and below y=0 so we can easily see which is which, looks better hopefully. 
    ax.axhspan(0, 1, facecolor='#34a1eb', alpha=0.03)  # Above y=0
    ax.axhspan(-1, 0, facecolor='#fcad03', alpha=0.03)  # Below y=0

    # Add text labels for primary and secondary signatures
    primary_sig_label = "More active in " + primary_sig
    secondary_sig_label = "More active in " + secondary_sig
    ax.text(0.05, 0.95, primary_sig_label, transform=ax.transAxes, fontsize=12, va='top', ha='left')
    ax.text(0.05, 0.05, secondary_sig_label, transform=ax.transAxes, fontsize=12, va='bottom', ha='left')

    # Adjust the bottom margin
    plt.subplots_adjust(bottom=0.2)
    # Colour code the axis labels
    labels = ax.get_xticklabels()
    for label in labels:
        mutation_type = re.search(r'\[(\w>\w)\]', label.get_text()).group(1)
        label.set_color(mutation_colors[mutation_type])

    
    # Find the top 5 points furthest from zero in either direction
    top_5_indices = df['Proportion'].abs().nlargest(5).index
    for idx in top_5_indices:
        mutation_label = idx
        ax.text(idx, df.loc[idx, 'Proportion'], mutation_label, color='black', ha='center', va='bottom' if df.loc[idx, 'Proportion'] > 0 else 'top', rotation=45)
        
    # Save the plot to a bytes buffer
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    # Encode the plot as base64 for embedding in HTML
    plot_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    # Close the figure to free up memory
    plt.close(fig)

    # Plot for the primary signature
    primary_sig_plot = plot_individual_signature(data[primary_sig], primary_sig, primary_sig)

    # Plot for the secondary signature
    secondary_sig_plot = plot_individual_signature(data[secondary_sig], secondary_sig, secondary_sig)

    return plot_base64, primary_sig_plot, secondary_sig_plot

def plot_individual_signature(signature, sig_name, title):
    # Create a DataFrame from the signature Series
    df = pd.DataFrame({'Proportion': signature})
    # Extract mutation types from the index using regular expressions
    df['MutationType'] = df.index.str.extract(r'\[(\w>\w)\]', expand=False)
    # Define colour mapping for mutation types
    mutation_colors = {
        'C>T': '#FF0000',
        'C>G': '#000000',
        'T>A': '#808080',
        'T>G': '#FFC0CB',
        'T>C': '#00FF00',
        'C>A': '#0000FF'
    }
    # Sort the DataFrame by 'MutationType'
    df['MutationType'] = pd.Categorical(df['MutationType'], categories=['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'])
    df = df.sort_values('MutationType')
    # Create a new figure
    fig, ax = plt.subplots(figsize=(10, 3))
    # Create a grouped bar plot
    sns.barplot(x=df.index, y='Proportion', data=df, hue='MutationType', dodge=True, palette=mutation_colors, ax=ax)
    ax.set_xlabel('Trinucleotide Context')
    ax.set_ylabel('Proportion of Single Base Substitutions')
    ax.set_title(title)
    ax.legend().remove()
    plt.xticks(rotation=90)
    plt.xticks(fontsize=6)
    # Adjust the bottom margin
    plt.subplots_adjust(bottom=0.5)
    # Colour code the axis labels
    labels = ax.get_xticklabels()
    for label in labels:
        mutation_type = re.search(r'\[(\w>\w)\]', label.get_text()).group(1)
        label.set_color(mutation_colors[mutation_type])

    # Save the plot to a bytes buffer
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    # Encode the plot as base64 for embedding in HTML
    plot_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    # Close the figure to free up memory
    plt.close(fig)
    return plot_base64

@app.route('/')
def index():
    signatures = {
        'SBS1': 'Spontaneous deamination of 5-methylcytosine (clock-like signature)',
        'SBS2': 'Activity of APOBEC family of cytidine deaminases',
        'SBS3': 'Defective homologous recombination DNA damage repair',
        'SBS4': 'Tobacco smoking',
        'SBS5': 'Unknown (clock-like signature)',
        'SBS6': 'Defective DNA mismatch repair',
        'SBS7a': 'Ultraviolet light exposure',
        'SBS7b': 'Ultraviolet light exposure',
        'SBS7c': 'Ultraviolet light exposure',
        'SBS7d': 'Ultraviolet light exposure',
        'SBS8': 'Unknown',
        'SBS9': 'Polimerase eta somatic hypermutation activity',
        'SBS10a': 'Polymerase epsilon exonuclease domain mutations',
        'SBS10b': 'Polymerase epsilon exonuclease domain mutations',
        'SBS10c': 'Defective POLD1 proofreading',
        'SBS10d': 'Defective POLD1 proofreading',
        'SBS11': 'Temozolomide treatment',
        'SBS12': 'Unknown',
        'SBS13': 'Activity of APOBEC family of cytidine deaminases',
        'SBS14': 'Concurrent polymerase epsilon mutation and defective DNA mismatch repair',
        'SBS15': 'Defective DNA mismatch repair',
        'SBS16': 'Unknown',
        'SBS17a': 'Unknown',
        'SBS17b': 'Unknown',
        'SBS18': 'Damage by reactive oxygen species',
        'SBS19': 'Unknown',
        'SBS20': 'Concurrent POLD1 mutations and defective DNA mismatch repair',
        'SBS21': 'Defective DNA mismatch repair',
        'SBS22a': 'Aristolochic acid exposure',
        'SBS22b': 'Aristolochic acid exposure',
        'SBS23': 'Unknown',
        'SBS24': 'Aflatoxin exposure',
        'SBS25': 'Chemotherapy treatment',
        'SBS26': 'Defective DNA mismatch repair',
        'SBS28': 'Unknown',
        'SBS29': 'Tobacco chewing',
        'SBS30': 'Defective DNA base excision repair due to NTHL1 mutations',
        'SBS31': 'Platinum chemotherapy treatment',
        'SBS32': 'Azathioprine treatment',
        'SBS33': 'Unknown',
        'SBS34': 'Unknown',
        'SBS35': 'Platinum chemotherapy treatment',
        'SBS36': 'Defective DNA base excision repair due to MUTYH mutations',
        'SBS37': 'Unknown',
        'SBS38': 'Indirect effect of ultraviolet light',
        'SBS39': 'Unknown',
        'SBS40a': 'Unknown',
        'SBS40b': 'Unknown',
        'SBS40c': 'Unknown',
        'SBS41': 'Unknown',
        'SBS42': 'Haloalkane exposure',
        'SBS44': 'Defective DNA mismatch repair',
        'SBS84': 'Activity of activation-induced cytidine deaminase (AID)',
        'SBS85': 'Indirect effects of activation-induced cytidine deaminase (AID)',
        'SBS86': 'Unknown chemotherapy treatment',
        'SBS87': 'Thiopurine chemotherapy treatment',
        'SBS88': 'Colibactin exposure (E.coli bacteria carrying pks pathogenicity island)',
        'SBS89': 'Unknown',
        'SBS90': 'Duocarmycin exposure',
        'SBS91': 'Unknown',
        'SBS92': 'Tobacco smoking',
        'SBS93': 'Unknown',
        'SBS94': 'Unknown',
        'SBS96': 'Unknown',
        'SBS97': 'Unknown',
        'SBS98': 'Unknown',
        'SBS99': 'Melphalan exposure'
    }
    return render_template('index.html', signatures=signatures)

@app.route('/plot', methods=['POST'])
def generate_plot():
    primary_signature = request.form['primary_signature']
    secondary_signature = request.form['secondary_signature']
    subtracted_signature = subtract_signatures(primary_signature, secondary_signature)
    plot_base64, primary_sig_plot, secondary_sig_plot = plot_signature(subtracted_signature, primary_signature, secondary_signature, f'{primary_signature} - {secondary_signature}')
    return jsonify({'plot_base64': plot_base64, 'primary_sig_plot': primary_sig_plot, 'secondary_sig_plot': secondary_sig_plot})

if __name__ == '__main__':
    app.run(debug=True)
