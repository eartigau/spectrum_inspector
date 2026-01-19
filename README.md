# NIRPS/SPIROU Spectrum Inspector 🔭

A simple tool to visually inspect NIRPS and SPIROU spectra order-by-order, comparing observations with templates and atmospheric transmission.

![Sample plot of Order 35](sample_order35.png)
*Example: Order 35 of TOI-4552 showing oxygen airglow fluorescence contamination in the residuals. This is the kind of artifact you should look for!*

---

## 📥 Installation

### Step 1: Download the code

```bash
git clone https://github.com/eartigau/spectrum_inspector.git
cd spectrum_inspector
```

Or if you don't have git, download the ZIP from: https://github.com/eartigau/spectrum_inspector/archive/refs/heads/main.zip

### Step 2: Install Python dependencies

```bash
pip install numpy matplotlib astropy scipy
```

Or with conda:
```bash
conda install numpy matplotlib astropy scipy
```

### Step 3: Copy the TAPAS file

```bash
# Adjust the path to YOUR lbl installation
cp /path/to/your/lbl/models/tapas_lbl.fits .
```

---

## 📋 What This Tool Does

This script creates PDF plots showing:
- **Page 1: Summary page** with observation info (BERV, systemic velocity, V_tot, files used, etc.)
- **Subsequent pages: One spectral order per page** with:
  - **Your observed spectrum** (black) vs **template spectrum** (red)
  - **Residuals** (observation - template) in orange
  - **TAPAS atmospheric transmission** (water in blue, other gases in green)
  - **OH sky emission lines** (cyan vertical lines)

---

## 🚀 Quick Start

### The Simplest Command

```bash
python inspect_spectrum.py YOUR_SPECTRUM.fits
```

That's it! The script will:
1. Auto-detect the template file based on the object name in your FITS header
2. Generate a multi-page PDF with ALL spectral orders
3. Save it as `NIRPS_OBJECTNAME_DATE_all_orders.pdf`

---

## 📦 Required Files

You need **THREE** files in the same folder as the script:

| File | What it is | Where to get it |
|------|-----------|-----------------|
| `YOUR_SPECTRUM.fits` | Your NIRPS/SPIROU observation | From your data reduction |
| `Template_s1dv_OBJECTNAME_sc1d_v_file_A.fits` | Template spectrum | From LBL reduction |
| `tapas_lbl.fits` | Atmospheric transmission | See Installation Step 3 |

---

## 📖 Usage Examples

```bash
python inspect_spectrum.py NIRPS.2024-09-28T23:54:06.014t.fits
```

**Output:** `NIRPS_TOI4552_2024-09-28_all_orders.pdf` (one page per order)

### Example 2: Specify the template explicitly

```bash
python inspect_spectrum.py NIRPS.2024-09-28T23:54:06.014t.fits Template_s1dv_TOI4552_sc1d_v_file_A.fits
```

### Example 3: Look at ONE specific order

```bash
python inspect_spectrum.py NIRPS.2024-09-28T23:54:06.014t.fits --order 50
```

**Output:** Opens a plot window AND saves `NIRPS_TOI4552_2024-09-28_order50.pdf`

### Example 4: Look at a RANGE of orders (e.g., orders 40 to 60)

```bash
python inspect_spectrum.py NIRPS.2024-09-28T23:54:06.014t.fits --order 40 60
```

**Output:** `NIRPS_TOI4552_2024-09-28_orders40-60.pdf` (21 pages)

### Example 5: Just save PDF, don't show plot window

```bash
python inspect_spectrum.py NIRPS.2024-09-28T23:54:06.014t.fits --order 50 --no-show
```

### Example 6: Hide the OH lines and TAPAS (cleaner plot)

```bash
python inspect_spectrum.py NIRPS.2024-09-28T23:54:06.014t.fits --no-oh --no-tapas
```

### Example 7: SPIROU data instead of NIRPS

```bash
python inspect_spectrum.py SPIROU_spectrum.fits --instrument SPIROU
```

---

## 🎛️ All Options

| Option | What it does | Default |
|--------|--------------|---------|
| `file` | Your spectrum FITS file | **REQUIRED** |
| `template` | Template spectrum file | Auto-detected from FITS header |
| `--tapas FILE` | TAPAS transmission file | `tapas_lbl.fits` |
| `--instrument` | `NIRPS` or `SPIROU` | `NIRPS` |
| `--order N` | Plot only order N | All orders |
| `--order N M` | Plot orders N through M | All orders |
| `--no-oh` | Hide OH emission lines | Show them |
| `--no-tapas` | Hide TAPAS transmission | Show it |
| `--no-save` | Don't save PDF | Save PDF |
| `--no-show` | Don't open plot window | Show plot (single order only) |

---

##  Understanding the Plots

### Page 1: Summary Page
Contains all the key observation metadata:
- **Target name** and **instrument**
- **Date of observation**
- **BERV** (Barycentric Earth Radial Velocity)
- **Systemic velocity** (from `HIERARCH ESO TEL TARG RADVEL`)
- **V_tot** = BERV - Vsys (total velocity shift applied)
- **Files used** (spectrum and template)
- **Wavelength range** and number of orders

### Order Pages (Page 2+)

Each page shows TWO panels:

### Top Panel
- **Black line:** Your normalized observed spectrum
- **Red line:** Template spectrum (shifted to your BERV)
- **Blue line:** TAPAS water vapor transmission
- **Green line:** TAPAS other gases (O₂, CO₂, etc.)
- **Cyan vertical lines:** OH sky emission lines (labeled)

### Bottom Panel
- **Orange line:** Residuals (observed - template)
- Good data should scatter around 0
- Systematic features = something interesting (or a problem!)

---

## ❓ Common Issues

### "Template file not found"

The script couldn't find a template matching your object. Either:
1. Provide the template explicitly as the second argument
2. Make sure your template file is named `Template_s1dv_OBJECTNAME_sc1d_v_file_A.fits`

### "TAPAS file not found"

Copy `tapas_lbl.fits` from your LBL installation:
```bash
cp /your/lbl/path/models/tapas_lbl.fits .
```

### "Order X has no valid data"

That order has NaN values (probably at the edge of the detector). The script automatically skips these orders.

### The PDF has fewer pages than expected

Orders with all-NaN data are automatically skipped. This is normal for edge orders.

---

## 🔧 Using in Python Scripts

You can also import the functions:

```python
from inspect_spectrum import inspect_all_orders, inspect_single_order

# Generate PDF for all orders
inspect_all_orders('spectrum.fits', 'template.fits')

# Generate PDF for orders 30-50
inspect_all_orders('spectrum.fits', 'template.fits', order_min=30, order_max=50)

# Plot a single order (shows plot + saves PDF)
inspect_single_order('spectrum.fits', 'template.fits', order=45)

# Plot without saving
inspect_single_order('spectrum.fits', 'template.fits', order=45, save_pdf=False)
```

---

## 📬 Questions?

Contact: etienne.artigau@umontreal.ca

---
