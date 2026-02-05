#!/usr/bin/env python
"""
Compare two spectra from different pipelines (APERO/ESO).

Usage:
    python compare_spectra.py file1.fits file2.fits template.fits
    python compare_spectra.py file1.fits file2.fits template.fits --order 50
    python compare_spectra.py file1.fits file2.fits template.fits --order 40 60
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import argparse

from inspect_spectrum import (
    load_spectrum, load_template, load_tapas, get_oh_lines,
    doppler, robust_polyfit, get_file_type
)


def find_matching_order(spectrum1, spectrum2, order1):
    """
    Find the order in spectrum2 that best matches the wavelength range of order1 in spectrum1.
    
    Parameters
    ----------
    spectrum1, spectrum2 : dict
        Spectrum data dictionaries
    order1 : int
        Order index in spectrum1
    
    Returns
    -------
    order2 : int or None
        Matching order index in spectrum2, or None if no good match found
    """
    wave1 = spectrum1['wave'][order1]
    valid1 = np.isfinite(wave1)
    if not np.any(valid1):
        return None
    
    wave1_min = np.nanmin(wave1[valid1])
    wave1_max = np.nanmax(wave1[valid1])
    wave1_center = (wave1_min + wave1_max) / 2
    
    best_order = None
    best_overlap = 0
    
    for order2 in range(spectrum2['n_orders']):
        wave2 = spectrum2['wave'][order2]
        valid2 = np.isfinite(wave2)
        if not np.any(valid2):
            continue
        
        wave2_min = np.nanmin(wave2[valid2])
        wave2_max = np.nanmax(wave2[valid2])
        
        # Check if wavelength ranges overlap
        overlap_min = max(wave1_min, wave2_min)
        overlap_max = min(wave1_max, wave2_max)
        
        if overlap_max > overlap_min:
            overlap = overlap_max - overlap_min
            # Prefer orders with better overlap
            if overlap > best_overlap:
                best_overlap = overlap
                best_order = order2
    
    return best_order


def plot_unmatched_order(spectrum_data, template_interp, order, file_type):
    """
    Plot a single spectrum order when there's no match in the other file.
    Shows a warning message on the plot.
    """
    from inspect_spectrum import doppler, robust_polyfit
    
    wave = spectrum_data['wave']
    flux = spectrum_data['flux']
    blaze = spectrum_data['blaze']
    berv = spectrum_data['berv']
    wave_in_rest_frame = spectrum_data.get('wave_in_rest_frame', False)
    needs_blaze_correction = spectrum_data.get('needs_blaze_correction', False)
    
    wave_ord = wave[order]
    
    # For APERO files: divide flux by blaze to remove blaze shape
    # For ESO files: blaze is unity, no effect
    flux_blaze_corrected = flux[order] / blaze[order]
    
    # Check if flux is valid
    if np.all(~np.isfinite(flux_blaze_corrected)):
        return None, None
    flux_median = np.nanmedian(flux_blaze_corrected)
    if not np.isfinite(flux_median) or flux_median == 0:
        return None, None
    
    # Normalize
    plot_tmp = flux_blaze_corrected / flux_median
    
    # Get template
    if wave_in_rest_frame:
        template_tmp = template_interp[0](wave_ord)
        template_mask = template_interp[1](wave_ord) > 0.5
    else:
        template_tmp = template_interp[0](doppler(wave_ord, -berv))
        template_mask = template_interp[1](doppler(wave_ord, -berv)) > 0.5
    template_tmp[~template_mask] = np.nan
    template_tmp /= np.nanmedian(template_tmp)
    
    # Apply blaze correction if needed
    if needs_blaze_correction:
        ratio = plot_tmp / template_tmp
        blaze_fit = robust_polyfit(wave_ord, ratio, degree=4, sigma_clip=3.0, max_iter=10)
        if np.any(np.isfinite(blaze_fit)):
            plot_tmp = plot_tmp / blaze_fit
            plot_tmp /= np.nanmedian(plot_tmp[np.isfinite(plot_tmp)])
    
    # Find valid domain
    flux_valid = np.isfinite(plot_tmp)
    if not np.any(flux_valid):
        return None, None
    
    valid_indices = np.where(flux_valid)[0]
    idx_min, idx_max = valid_indices[0], valid_indices[-1]
    wave_min, wave_max = wave_ord[idx_min], wave_ord[idx_max]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(18, 6), nrows=1)
    
    ax.plot(wave_ord, plot_tmp, alpha=0.7, color='black', label=f'{file_type} Order {order}', linewidth=0.8, rasterized=True)
    ax.plot(wave_ord, template_tmp, alpha=0.5, color='red', label='Template', linewidth=0.8, rasterized=True)
    
    # Add prominent warning
    ax.text(0.5, 0.5, '⚠ NO MATCHING ORDER IN OTHER FILE ⚠', 
            transform=ax.transAxes, fontsize=20, color='red', alpha=0.7,
            ha='center', va='center', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
    
    ax.set_xlim(wave_min, wave_max)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Normalized Flux')
    ax.set_title(f"{spectrum_data['file_type']} | Order {order} | {wave_min:.1f}-{wave_max:.1f} nm | NO MATCH")
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(alpha=0.3)
    
    return fig, ax


def plot_comparison_order(spectrum1, spectrum2, template_interp, order1, order2,
                          tapas=None, wave_oh=None, label_oh=None,
                          show_oh=True, show_tapas=True):
    """
    Plot comparison of two spectra for a single order.
    
    Parameters
    ----------
    spectrum1, spectrum2 : dict
        Spectrum data dictionaries from load_spectrum
    template_interp : tuple
        Template interpolators from load_template
    order1 : int
        Order number in spectrum1
    order2 : int
        Order number in spectrum2 (matched by wavelength)
    
    Returns
    -------
    fig, ax : matplotlib figure and axes, or (None, None) if order is invalid
    """
    
    def process_spectrum(spectrum_data, template_interp, order):
        """Process a single spectrum for plotting."""
        wave = spectrum_data['wave']
        flux = spectrum_data['flux']
        blaze = spectrum_data['blaze']
        berv = spectrum_data['berv']
        wave_in_rest_frame = spectrum_data.get('wave_in_rest_frame', False)
        needs_blaze_correction = spectrum_data.get('needs_blaze_correction', False)
        
        wave_ord = wave[order]
        
        # For APERO files: divide flux by blaze to remove blaze shape
        # For ESO files: blaze is unity, no effect
        flux_blaze_corrected = flux[order] / blaze[order]
        
        # Check if flux is all NaN
        if np.all(~np.isfinite(flux_blaze_corrected)):
            return None, None, None, None
        flux_median = np.nanmedian(flux_blaze_corrected)
        if not np.isfinite(flux_median) or flux_median == 0:
            return None, None, None, None
        
        # Normalize the blaze-corrected flux
        plot_tmp = flux_blaze_corrected / flux_median
        
        # Get template at correct wavelengths
        if wave_in_rest_frame:
            template_tmp = template_interp[0](wave_ord)
            template_mask = template_interp[1](wave_ord) > 0.5
        else:
            template_tmp = template_interp[0](doppler(wave_ord, -berv))
            template_mask = template_interp[1](doppler(wave_ord, -berv)) > 0.5
        template_tmp[~template_mask] = np.nan
        
        # Template is already flat - just normalize by median
        template_tmp /= np.nanmedian(template_tmp)
        
        # Store original for reference
        plot_tmp_original = plot_tmp.copy()
        
        # Apply blaze correction if needed
        if needs_blaze_correction:
            ratio = plot_tmp / template_tmp
            blaze_fit = robust_polyfit(wave_ord, ratio, degree=4, sigma_clip=3.0, max_iter=10)
            if np.any(np.isfinite(blaze_fit)):
                plot_tmp = plot_tmp / blaze_fit
                plot_tmp /= np.nanmedian(plot_tmp[np.isfinite(plot_tmp)])
        
        return wave_ord, plot_tmp, template_tmp, plot_tmp_original
    
    # Process both spectra with their respective order indices
    result1 = process_spectrum(spectrum1, template_interp, order1)
    result2 = process_spectrum(spectrum2, template_interp, order2)
    
    if result1[0] is None or result2[0] is None:
        return None, None
    
    wave1, flux1, template1, orig1 = result1
    wave2, flux2, template2, orig2 = result2
    
    # Find valid domain for both spectra
    valid1 = np.isfinite(flux1)
    valid2 = np.isfinite(flux2)
    
    if not np.any(valid1) or not np.any(valid2):
        return None, None
    
    idx1_min, idx1_max = np.where(valid1)[0][[0, -1]]
    idx2_min, idx2_max = np.where(valid2)[0][[0, -1]]
    
    wave_min = max(wave1[idx1_min], wave2[idx2_min])
    wave_max = min(wave1[idx1_max], wave2[idx2_max])
    
    # Create figure with 3 panels
    fig, ax = plt.subplots(figsize=(18, 10), nrows=3, sharex=True,
                           gridspec_kw={'height_ratios': [2, 1, 1]})
    
    # Labels for the two spectra (include order numbers since they may differ)
    label1 = f"{spectrum1['file_type']} ord{order1} ({spectrum1['fits_file'][:25]}...)"
    label2 = f"{spectrum2['file_type']} ord{order2} ({spectrum2['fits_file'][:25]}...)"
    
    # Panel 0: Both spectra and template
    ax[0].plot(wave1, flux1, alpha=0.7, color='black', label=label1, linewidth=0.8, rasterized=True)
    ax[0].plot(wave2, flux2, alpha=0.7, color='blue', label=label2, linewidth=0.8, rasterized=True)
    ax[0].plot(wave1, template1, alpha=0.5, color='red', label='Template', linewidth=0.8, rasterized=True)
    
    # Panel 1: Individual residuals (spectrum - template)
    residuals1 = flux1 - template1
    residuals2 = flux2 - template2
    ax[1].plot(wave1, residuals1, alpha=0.7, color='black', label=f'{label1} - Template', linewidth=0.8, rasterized=True)
    ax[1].plot(wave2, residuals2, alpha=0.7, color='blue', label=f'{label2} - Template', linewidth=0.8, rasterized=True)
    
    # Panel 2: Difference between the two spectra
    # Interpolate spectrum2 onto spectrum1's wavelength grid for direct comparison
    from scipy.interpolate import InterpolatedUnivariateSpline as IUS
    
    # Only interpolate where both are valid
    valid2_interp = np.isfinite(flux2)
    if np.sum(valid2_interp) > 10:
        flux2_interp = IUS(wave2[valid2_interp], flux2[valid2_interp], k=1, ext=1)(wave1)
        diff = flux1 - flux2_interp
        ax[2].plot(wave1, diff, alpha=0.7, color='purple', 
                   label=f'{spectrum1["file_type"]} - {spectrum2["file_type"]}', linewidth=0.8, rasterized=True)
        ax[2].axhline(0, color='grey', linestyle='--', alpha=0.5)
    
    # Add OH lines
    if show_oh and wave_oh is not None:
        keep = (wave_oh >= wave_min) & (wave_oh <= wave_max)
        for w in wave_oh[keep]:
            for a in ax:
                a.axvline(w, color='cyan', alpha=1.0, linewidth=1.0)
    
    # Add TAPAS
    if show_tapas and tapas is not None:
        keep = (tapas['WAVELENGTH'] >= wave_min) & (tapas['WAVELENGTH'] <= wave_max)
        tapas_trim = tapas[keep]
        if len(tapas_trim) > 0:
            ax[0].plot(tapas_trim['WAVELENGTH'], tapas_trim['ABSO_OTHERS'], 
                      color='green', alpha=0.9, label='TAPAS Other', linewidth=1.0, rasterized=True)
            ax[0].plot(tapas_trim['WAVELENGTH'], tapas_trim['ABSO_WATER'],
                      color='blue', alpha=0.9, linestyle='--', label='TAPAS Water', linewidth=1.0, rasterized=True)
    
    # Set axis properties
    ax[0].set_xlim(wave_min, wave_max)
    
    # Y-limits for flux plot
    ymax_flux = 1.2 * max(np.nanpercentile(orig1, 90), np.nanpercentile(orig2, 90))
    ax[0].set_ylim(0, ymax_flux)
    
    # Y-limits for residuals - symmetric and capped
    for i in [1, 2]:
        data = residuals1 if i == 1 else diff if 'diff' in dir() else residuals1
        valid_data = data[np.isfinite(data)]
        if len(valid_data) > 0:
            ymin, ymax = np.min(valid_data), np.max(valid_data)
            ymin, ymax = max(ymin, -0.5), min(ymax, 0.5)
            ylim = (np.abs(ymin) + np.abs(ymax)) / 2
            ax[i].set_ylim(-ylim, ylim)
        else:
            ax[i].set_ylim(-0.5, 0.5)
    
    # Labels and legends
    ax[0].legend(loc='upper right', fontsize=7)
    ax[1].legend(loc='upper right', fontsize=7)
    ax[2].legend(loc='upper right', fontsize=7)
    
    ax[0].set_ylabel('Normalized Flux')
    ax[1].set_ylabel('Residuals (vs Template)')
    ax[2].set_ylabel('Difference')
    ax[2].set_xlabel('Wavelength (nm)')
    
    # Title - show both order numbers if they differ
    if order1 == order2:
        title = f"Comparison | Order {order1} | {wave_min:.1f}-{wave_max:.1f} nm"
    else:
        title = f"Comparison | {spectrum1['file_type']} ord{order1} vs {spectrum2['file_type']} ord{order2} | {wave_min:.1f}-{wave_max:.1f} nm"
    ax[0].set_title(title)
    
    for a in ax:
        a.grid(alpha=0.3)
    
    return fig, ax


def compare_all_orders(file1, file2, template_file, tapas_file='reference_data/tapas_lbl.fits',
                       instrument='NIRPS', order_min=None, order_max=None,
                       show_oh=True, show_tapas=True):
    """
    Create multi-page PDF comparing two spectra across all orders.
    
    Loops through orders of file1 and finds matching orders in file2 by wavelength.
    """
    # Load both spectra
    spectrum1 = load_spectrum(file1, instrument)
    spectrum2 = load_spectrum(file2, instrument)
    template_interp = load_template(template_file)
    
    # Use file1's orders as reference
    n_orders1 = spectrum1['n_orders']
    
    if order_min is None:
        order_min = 0
    if order_max is None:
        order_max = n_orders1 - 1
    
    order_min = max(0, order_min)
    order_max = min(n_orders1 - 1, order_max)
    
    # Load OH lines and TAPAS
    wave_oh, label_oh = None, None
    tapas = None
    
    if show_oh:
        full_wave_min = min(np.nanmin(spectrum1['wave']), np.nanmin(spectrum2['wave']))
        full_wave_max = max(np.nanmax(spectrum1['wave']), np.nanmax(spectrum2['wave']))
        wave_oh, label_oh = get_oh_lines(full_wave_min, full_wave_max)
    
    if show_tapas:
        tapas = load_tapas(tapas_file)
    
    # Generate PDF name
    timestamp1 = spectrum1['date_obs'].replace(':', '-')
    timestamp2 = spectrum2['date_obs'].replace(':', '-')
    pdf_path = f"comparison_{spectrum1['file_type']}_{spectrum2['file_type']}_{spectrum1['target'].replace(' ', '_')}_{timestamp1}_vs_{timestamp2}_orders{order_min}-{order_max}.pdf"
    
    print(f"Comparing {spectrum1['file_type']} ({spectrum1['n_orders']} orders) vs {spectrum2['file_type']} ({spectrum2['n_orders']} orders)")
    print(f"Looping through file1 orders {order_min} to {order_max}, matching by wavelength...")
    
    n_matched = 0
    n_unmatched = 0
    with PdfPages(pdf_path) as pdf:
        for order1 in range(order_min, order_max + 1):
            # Find matching order in spectrum2 by wavelength
            order2 = find_matching_order(spectrum1, spectrum2, order1)
            
            if order2 is None:
                print(f"  Order {order1}: no matching order in file2, adding warning page...")
                # Create a warning page for unmatched order
                fig, ax = plot_unmatched_order(spectrum1, template_interp, order1, spectrum1['file_type'])
                if fig is not None:
                    n_unmatched += 1
                    plt.tight_layout()
                    pdf.savefig(fig, dpi=150)
                    plt.close(fig)
                continue
            
            print(f"  Processing file1 order {order1} <-> file2 order {order2}...", end='\r')
            
            fig, ax = plot_comparison_order(
                spectrum1, spectrum2, template_interp, order1, order2,
                tapas=tapas, wave_oh=wave_oh, label_oh=label_oh,
                show_oh=show_oh, show_tapas=show_tapas
            )
            
            if fig is None:
                continue
            
            n_matched += 1
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close(fig)
    
    print(f"\nMatched {n_matched} orders, {n_unmatched} unmatched")
    print(f"Saved: {pdf_path}")
    return pdf_path


def main():
    parser = argparse.ArgumentParser(
        description='Compare two spectra from different pipelines',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python compare_spectra.py file1.fits file2.fits template.fits
  python compare_spectra.py file1.fits file2.fits template.fits --order 50
  python compare_spectra.py file1.fits file2.fits template.fits --order 40 60
"""
    )
    parser.add_argument('file1', help='First FITS file')
    parser.add_argument('file2', help='Second FITS file')
    parser.add_argument('template', help='Template spectrum file')
    parser.add_argument('--tapas', default='reference_data/tapas_lbl.fits',
                        help='TAPAS transmission file')
    parser.add_argument('--instrument', default='NIRPS', choices=['NIRPS', 'SPIROU'],
                        help='Instrument name')
    parser.add_argument('--order', type=int, nargs='+',
                        help='Order(s) to plot. Single value for one order, two values for range.')
    parser.add_argument('--no-oh', action='store_true', help='Hide OH lines')
    parser.add_argument('--no-tapas', action='store_true', help='Hide TAPAS')
    
    args = parser.parse_args()
    
    # Check files exist
    for f in [args.file1, args.file2, args.template]:
        if not os.path.exists(f):
            print(f"Error: File '{f}' not found.")
            exit(1)
    
    # Validate file types
    for f in [args.file1, args.file2]:
        file_type = get_file_type(f)
        if file_type is None:
            print(f"Error: '{f}' is not a valid APERO (t.fits) or ESO (r.) file.")
            exit(1)
        print(f"  {os.path.basename(f)}: {file_type}")
    
    if args.order is None:
        compare_all_orders(
            args.file1, args.file2, args.template, args.tapas,
            instrument=args.instrument,
            show_oh=not args.no_oh,
            show_tapas=not args.no_tapas
        )
    elif len(args.order) == 1:
        compare_all_orders(
            args.file1, args.file2, args.template, args.tapas,
            instrument=args.instrument,
            order_min=args.order[0],
            order_max=args.order[0],
            show_oh=not args.no_oh,
            show_tapas=not args.no_tapas
        )
    else:
        compare_all_orders(
            args.file1, args.file2, args.template, args.tapas,
            instrument=args.instrument,
            order_min=args.order[0],
            order_max=args.order[1],
            show_oh=not args.no_oh,
            show_tapas=not args.no_tapas
        )


if __name__ == '__main__':
    main()
