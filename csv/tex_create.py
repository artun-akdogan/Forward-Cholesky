import pandas as pd
import os
import glob

def csv_to_latex(
    df,
    caption="Table",
    label="table:label",
    index=False,
    column_format=None,
    escape=True,
):
    """
    Convert DataFrame to Overleaf-compatible LaTeX table.
    """
    # Create default alignment if not provided
    if column_format is None:
        column_format = "l" * len(df.columns)

    latex = df.to_latex(
        index=index,
        caption=caption,
        label=label,
        escape=escape,
        column_format=column_format,
        longtable=False,
        bold_rows=False,
        na_rep="-",
        float_format="%.6g",
    )
    return latex


# ============================================
# Convert ALL CSV FILES in a folder
# ============================================

def convert_all_csv_to_latex(folder="."):
    """
    Converts every CSV in a folder to an Overleaf-compatible LaTeX table.
    Creates a .tex file for each CSV.
    """
    csv_files = glob.glob(os.path.join(folder, "*.csv"))

    if not csv_files:
        print("No CSV files found in the folder:", folder)
        return

    print(f"Found {len(csv_files)} CSV files.\n")

    for path in csv_files:
        filename = os.path.basename(path)
        name_no_ext = os.path.splitext(filename)[0]

        # Load CSV
        df = pd.read_csv(path)

        # Generate LaTeX code
        latex = csv_to_latex(
            df,
            caption=f"Results from {filename}",
            label=f"table:{name_no_ext}",
            index=False,
        )

        # Output .tex file
        tex_path = os.path.join(folder, f"{name_no_ext}.tex")
        with open(tex_path, "w") as f:
            f.write(latex)

        print(f"✔ Converted: {filename} → {name_no_ext}.tex")

    print("\nAll files converted successfully.")


# ------------------------
# Run conversion
# ------------------------

convert_all_csv_to_latex("all_csv")  # current folder
